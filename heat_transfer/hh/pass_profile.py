from heat_transfer.heat_rate import WallResistance, GasResistance, WaterResistance, FoulingResistance, TotalResistance
from heat_transfer.fluid_props.GasProps import HotFlueGas
from heat_transfer.hh.afad import HeatLossCalculator, GasTemperature
from heat_transfer.convection.pressure_drop import DarcyFriction, TubePressureDrop
from heat_transfer.convection.flow import GasVelocityCalculator, GasNusselt, GasReynoldsNumber, PrandtlCalculator
from heat_transfer.convection.convection import GasConvectiveCoefficient

class PassRunner:
    def __init__(self,
                 gas: HotFlueGas,
                 heat_loss_calc: HeatLossCalculator,
                 gas_temp: GasTemperature,
                 darcy: DarcyFriction,
                 tube_dp: TubePressureDrop,
                 r_water: WaterResistance,           # instance
                 r_wall: WallResistance,             # instance
                 r_foul_inner: FoulingResistance,    # instance
                 r_foul_outer: FoulingResistance,    # instance
                 h_rad_func,                         # fn(T,P,X)-> h_rad [W/m^2-K]
                 cp_func,                            # fn(T,P,X)-> cp [J/kg-K]
                 mass_func,                          # fn(T,P,rho,dx)-> mass [kg]
                 h_func, inv_h_func,                 # h(T)->J/kg, inv_h(h)->K
                 mdot, A, D, X,
                 rel_roughness, f0, tol, max_iter,
                 nusselt_n):
        self.gas = gas
        self.hl = heat_loss_calc
        self.gt = gas_temp
        self.darcy = darcy
        self.tube_dp = tube_dp
        self.r_water = r_water
        self.r_wall = r_wall
        self.r_fi = r_foul_inner
        self.r_fo = r_foul_outer
        self.h_rad = h_rad_func
        self.cp = cp_func
        self.mass = mass_func
        self.h = h_func
        self.inv_h = inv_h_func
        self.mdot = mdot
        self.A = A
        self.D = D
        self.X = X
        self.rr = rel_roughness
        self.f0, self.tol, self.max_iter = f0, tol, max_iter
        self.nusselt_n = nusselt_n

    def run(self, L, dx, T0, P0):
        x, T, P = [0.0], [T0], [P0]
        v_list, Re_list, Pr_list, Nu_list = [], [], [], []
        h_conv_list, f_list, dp_list = [], [], []
        qflux_list, qrate_list = [], []

        n = int(L / dx)
        Ti, Pi = T0, P0

        for i in range(n):
            xi = (i + 1) * dx

            k  = self.gas.thermal_conductivity(Ti, Pi, self.X)
            mu = self.gas.viscosity(Ti, Pi, self.X)
            rho= self.gas.density(Ti, Pi, self.X)
            cp = self.cp(Ti, Pi, self.X)

            vel_calc = GasVelocityCalculator(self.mdot, rho, self.A)
            v = vel_calc.velocity()

            Pr = PrandtlCalculator(mu, cp, k).compute()
            Re = GasReynoldsNumber(rho, v, self.D, mu).calculate()
            Nu = GasNusselt(Re, Pr, self.nusselt_n).calculate()
            h_conv = GasConvectiveCoefficient(k, self.D, Nu).h()
            h_rad = self.h_rad(Ti, Pi, self.X)

            r_water_pa = self.r_water.resistance_per_area()
            r_wall_pi  = self.r_wall.cyl_per_inner_area()
            r_gas_pa   = GasResistance(h_rad, h_conv).resistance_per_area()
            r_fi_pa    = self.r_fi.resistance_per_area()
            r_fo_pa    = self.r_fo.resistance_per_area()
            R_tot = TotalResistance(r_water_pa, r_wall_pi, r_gas_pa, r_fi_pa, r_fo_pa).resistance_per_area()

            # user supplies ΔT(x,T,P) externally if needed; here compute directly from coefficients
            # if you already have a deltaT_func, call it here and replace dT below
            dT = None  # placeholder when an external ΔT is injected elsewhere

            # if an external ΔT is not injected, compute q from existing classes requires ΔT
            # keep interface minimal: expect caller to set dT before run or inject a function
            # here we assume dT is provided per step via a callable attribute deltaT(x,T,P)
            q_flux = self.deltaT(xi, Ti, Pi) / R_tot
            q_rate = q_flux * self.A

            heat_loss = self.hl.heat_loss(q_rate, v, dx)
            m = self.mass(Ti, Pi, rho, dx)
            Ti = self.gt.new_temperature(Ti, heat_loss, m, self.h, self.inv_h)

            f = self.darcy.colebrook(Re, self.rr, self.f0, self.tol, self.max_iter)
            dp = self.tube_dp(dx, f, self.D, rho, v)
            Pi = Pi - dp

            x.append(xi); T.append(Ti); P.append(Pi)
            v_list.append(v); Re_list.append(Re); Pr_list.append(Pr); Nu_list.append(Nu)
            h_conv_list.append(h_conv); f_list.append(f); dp_list.append(dp)
            qflux_list.append(q_flux); qrate_list.append(q_rate)

        return {
            "x": x, "T": T, "P": P,
            "v": v_list, "Re": Re_list, "Pr": Pr_list, "Nu": Nu_list,
            "h_conv": h_conv_list, "f": f_list, "dp": dp_list,
            "q_flux": qflux_list, "q_rate": qrate_list
        }
