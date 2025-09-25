# heat_transfer/solve/continuous_pass.py
from __future__ import annotations
from dataclasses import dataclass
import math

# optional dense solution
try:
    from scipy.integrate import solve_ivp
    _HAVE_SCIPY = True
except Exception:
    _HAVE_SCIPY = False

# --- use your existing modules ---
from heat_transfer.flue_gas.dimensionless import GasReynoldsNumber, GasNusselt
from heat_transfer.flue_gas.convection import GasConvectiveCoefficient
from heat_transfer.flue_gas.flow import GasVelocityCalculator
from heat_transfer.flue_gas.pressure_drop import DarcyFriction, TubePressureDrop
from heat_transfer.flue_gas.properties import PrandtlCalculator
from heat_transfer.htc.gas_h_rad import GasRadiationCoefficient
from heat_transfer.htc.water_h_boil import PoolBoilingHTC
from heat_transfer.flux import WaterResistance, GasResistance, FoulingResistance, TotalResistance, Flux

# --- states and config ---
@dataclass
class GasState:
    T: float          # K
    P: float          # Pa
    m_dot: float      # kg/s
    X: dict[str, float]

@dataclass
class WaterState:
    P: float          # Pa
    m_dot: float      # kg/s
    h: float          # J/kg

@dataclass
class MarchConfig:
    L: float                  # pass length [m]
    eps_wall: float           # wall emissivity [-]
    eps_gas: float            # effective gas emissivity [-]
    roughness: float          # inner roughness [m]
    fouling_i: tuple[float, float]  # (delta_i [m], k_i [W/m/K])
    fouling_o: tuple[float, float]  # (delta_o [m], k_o [W/m/K])
    sigma: float = 5.670374419e-8
    nusselt: str = "db"       # "db" (Dittus–Boelter) for now
    N_eval: int = 401         # fallback RK4 points

@dataclass
class ContinuousProfile:
    x: list[float]
    Tg: list[float]
    Pg: list[float]
    hw: list[float]
    qpp: list[float]
    Tw: list[float]
    # simple interpolator
    def __call__(self, xq: float) -> dict[str, float]:
        xs = self.x
        if xq <= xs[0]: i = 0
        elif xq >= xs[-1]: i = len(xs)-2
        else:
            i = next(j for j in range(len(xs)-1) if xs[j] <= xq <= xs[j+1])
        t = (xq - xs[i]) / max(xs[i+1]-xs[i], 1e-30)
        lerp = lambda a,b: a + t*(b-a)
        return dict(
            Tg=lerp(self.Tg[i], self.Tg[i+1]),
            Pg=lerp(self.Pg[i], self.Pg[i+1]),
            hw=lerp(self.hw[i], self.hw[i+1]),
            qpp=lerp(self.qpp[i], self.qpp[i+1]),
            Tw=lerp(self.Tw[i], self.Tw[i+1]),
        )

# --- helpers binding to your geometry and property providers ---
def _geom_totals(geom: "PassGeomCalc") -> tuple[float,float,float,float,int,float]:
    # uses your PassGeomCalc.p (schemas.Pass)
    g = geom.p.geometry
    D = float(g.inner_diameter)
    ri = 0.5*D
    ro = float(g.wall.thickness) + ri
    Nt = int(g.number_of_tubes)
    A_flow = Nt * math.pi * (0.5*D)**2
    P_w = Nt * math.pi * D
    k_wall = float(g.wall.conductivity)
    return D, ri, ro, A_flow, Nt, k_wall

def _Rpp_wall_inner_basis(ri: float, ro: float, k: float) -> float:
    # correct cylindrical wall resistance per *inner area* basis
    return (ri / k) * math.log(ro/ri)

def _friction_factor(Re: float, rel_rough: float) -> float:
    if Re < 2300:
        return 64.0 / max(Re, 1e-12)
    # Newton step Colebrook using your class
    f0 = 0.02
    return DarcyFriction.colebrook(Re, rel_rough, f0=f0, tol=1e-8, max_iter=50)

# --- main marcher ---
def continuous_pass(
    gas0: GasState,
    water0: WaterState,
    geom: "PassGeomCalc",
    props_g: "GasProps",    # must expose k(T,P,X), mu(T,P,X), rho(T,P,X), cp(T,P,X)
    props_w: "WaterProps",  # must expose Tsat(P), rho_l(T), rho_v(T), mu_l(T), k_l(T), cp_l(T), sigma(T), h_fg(P)
    cfg: MarchConfig,
) -> tuple[GasState, WaterState, ContinuousProfile]:

    D, ri, ro, A_flow, _Nt, k_wall = _geom_totals(geom)
    L = cfg.L
    F_gw = 1.0  # tube core to wall view factor approximation

    delta_i, k_i = cfg.fouling_i
    delta_o, k_o = cfg.fouling_o

    vel = lambda rho: GasVelocityCalculator(gas0.m_dot, rho, A_flow).velocity()
    tube_dp = TubePressureDrop()

    diag: list[tuple[float,float,float,float,float,float]] = []  # x, Tg, Pg, hw, qpp, Tw

    def _hc(Tg, Pg, Tw):
        rho = props_g.rho(Tg, Pg, gas0.X)
        mu  = props_g.mu(Tg, Pg, gas0.X)
        k   = props_g.k(Tg, Pg, gas0.X)
        cp  = props_g.cp(Tg, Pg, gas0.X)
        v   = vel(rho)
        Re  = GasReynoldsNumber(rho, v, D, mu).calculate()
        Pr  = PrandtlCalculator(mu, cp, k).compute()
        if cfg.nusselt == "db":
            Nu = GasNusselt(Re, Pr, n=0.4).calculate()
        else:
            Nu = GasNusselt(Re, Pr, n=0.4).calculate()  # placeholder; extend when you add more correlations
        h_c = GasConvectiveCoefficient(k, D, Nu).h()
        return h_c, Re, v, rho, mu, Pr

    def _water_htc(Tw, P_wtr):
        Tsat = props_w.Tsat(P_wtr)
        dT = max(Tw - Tsat, 1e-9)
        mu_l = props_w.mu_l(Tsat)
        k_l  = props_w.k_l(Tsat)
        cp_l = props_w.cp_l(Tsat)
        rho_l = props_w.rho_l(Tsat)
        rho_v = props_w.rho_v(Tsat)
        sigma = props_w.sigma(Tsat)
        h_fg  = props_w.h_fg(P_wtr)
        Pr_l = mu_l*cp_l/k_l
        h_bo = PoolBoilingHTC.rohsenow_h(mu_l, h_fg, 9.80665, rho_l, rho_v, sigma, cp_l, dT, C_sf=0.013, Pr_l=Pr_l, n=1.7)
        return h_bo, Tsat

    def _qpp_and_Tw(Tg, Pg, hw):
        # fixed-point on Tw using your resistances
        Tw = Tg - 0.3*(Tg - props_w.Tsat(water0.P))
        hr_model = GasRadiationCoefficient()
        for _ in range(6):
            h_c, *_ = _hc(Tg, Pg, Tw)
            hr = hr_model.h(mean_temperature=Tg, emissivity=min(max(cfg.eps_gas*cfg.eps_wall, 0.0), 1.0), sigma=cfg.sigma)
            Rg = GasResistance.per_area(hr, h_c)
            Rw = _Rpp_wall_inner_basis(ri, ro, k_wall)
            Rfi = FoulingResistance.per_area(delta_i / max(k_i, 1e-30))
            Rfo = FoulingResistance.per_area(delta_o / max(k_o, 1e-30))
            h_bo, Tsat = _water_htc(Tw, water0.P)
            Rwtr = WaterResistance.per_area(h_bo)
            Rtot = TotalResistance.per_area(Rwtr, Rw, Rg, Rfi, Rfo)
            qpp = Flux.per_area(Tg - Tsat, Rtot)
            # update Tw from gas-side drop: q'' = (Tg - Tw)/R''_gas
            Tw = Tg - qpp * Rg
        return qpp, Tw

    P_wet = _Nt * math.pi * D  # total wetted perimeter

    def rhs(x, y):
        Tg, Pg, hw = y
        # heat flux and wall temperature
        qpp, Tw = _qpp_and_Tw(Tg, Pg, hw)
        # gas properties for cp and momentum
        cp_g = props_g.cp(Tg, Pg, gas0.X)
        rho = props_g.rho(Tg, Pg, gas0.X)
        mu  = props_g.mu(Tg, Pg, gas0.X)
        v   = vel(rho)
        Re  = GasReynoldsNumber(rho, v, D, mu).calculate()
        f = _friction_factor(Re, cfg.roughness/max(D,1e-30))
        dTgdx = - qpp * P_wet / (gas0.m_dot * cp_g)
        dhwdx = + qpp * P_wet /  max(water0.m_dot, 1e-30)
        dPdx  = - f * (rho * v * v / 2.0) / D
        diag.append((x, Tg, Pg, hw, qpp, Tw))
        return (dTgdx, dPdx, dhwdx)

    # integrate
    x0, x1 = 0.0, L
    y0 = (gas0.T, gas0.P, water0.h)

    if _HAVE_SCIPY:
        sol = solve_ivp(rhs, (x0, x1), y0, method="RK45", atol=1e-6, rtol=1e-6)
        xs = sol.t.tolist()
        Tg, Pg, hw = sol.y[0].tolist(), sol.y[1].tolist(), sol.y[2].tolist()
    else:
        # simple RK4 fallback
        N = max(5, cfg.N_eval)
        xs = [x0 + i*(L/(N-1)) for i in range(N)]
        Tg, Pg, hw = [y0[0]], [y0[1]], [y0[2]]
        for i in range(N-1):
            h = xs[i+1]-xs[i]
            y = (Tg[-1], Pg[-1], hw[-1])
            k1 = rhs(xs[i], y)
            k2 = rhs(xs[i]+0.5*h, tuple(y[j]+0.5*h*k1[j] for j in range(3)))
            k3 = rhs(xs[i]+0.5*h, tuple(y[j]+0.5*h*k2[j] for j in range(3)))
            k4 = rhs(xs[i]+h,     tuple(y[j]+h*k3[j]      for j in range(3)))
            yn = tuple(y[j] + (h/6.0)*(k1[j]+2*k2[j]+2*k3[j]+k4[j]) for j in range(3))
            Tg.append(yn[0]); Pg.append(yn[1]); hw.append(yn[2])

    # diagnostics aligned to xs length (simple padding if needed)
    if len(diag) < len(xs):
        last = diag[-1] if diag else (0, y0[0], y0[1], y0[2], 0.0, gas0.T)
        diag += [last]*(len(xs)-len(diag))
    qpp = [d[4] for d in diag[:len(xs)]]
    Tw  = [d[5] for d in diag[:len(xs)]]

    prof = ContinuousProfile(x=xs, Tg=Tg, Pg=Pg, hw=hw, qpp=qpp, Tw=Tw)

    gas_out = GasState(T=Tg[-1], P=Pg[-1], m_dot=gas0.m_dot, X=gas0.X)
    water_out = WaterState(P=water0.P, m_dot=water0.m_dot, h=hw[-1])
    return gas_out, water_out, prof
