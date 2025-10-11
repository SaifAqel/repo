# heat_transfer/functions/stage_solver.py
from dataclasses import dataclass
from math import pi, log
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from common.units import Q_, ureg
from heat_transfer.config.models import GasStream, WaterStream, FirePass, SmokePass, Reversal
from heat_transfer.functions.htc_water import HTCFunctions
from iapws import IAPWS97
import numpy as np

@dataclass
class HeatStageSolver:
    stage: FirePass | SmokePass | Reversal
    gas: GasStream
    water: WaterStream

    def _wall_flux(self, Tg: Q_, Tc: Q_) -> tuple[Q_, Q_, Q_]:
        # Geometry and wall
        ri = (self.stage.hot_side.inner_diameter/2).to("m")
        ro = (self.stage.hot_side.outer_diameter/2).to("m")
        kw = self.stage.hot_side.wall.conductivity.to("W/(m*K)")
        L  = self.stage.hot_side.inner_length

        # Hot-side convection (bulk based). Radiation depends on Twi.
        h_g = self.gas.convective_coefficient.to("W/(m^2*K)")

        Dh_shell = self.stage.cold_side.hydraulic_diameter
        # Root function in Twi
        def F(Twi_K: float) -> float:
            Twi = Q_(Twi_K, "K")
            h_r = self.gas.radiation_coefficient(Twi).to("W/(m^2*K)")
            qprime_hot = 2*pi*ri * (h_g + h_r) * (Tg - Twi)  # W/m
            # Conduction to outer wall
            Two = Twi - qprime_hot*log(ro/ri)/(2*pi*kw)
            # Cold-side HTC at Two
            qpp_hot = (qprime_hot / (2*pi*ro)).to("W/m^2")
            h_w = HTCFunctions.shell(self.water, Dh_shell, L, Two, qpp_hot, x=0.0)
            qprime_cold = 2*pi*ro * h_w * (Two - Tc)
            return (qprime_hot - qprime_cold).to("W/m").magnitude

        # Bracket: Tc < Twi < Tg

        # Bracket: Tc < Twi < Tg
        a = Tc.to("K").magnitude + 1e-6
        b = Tg.to("K").magnitude - 1e-6
        Fa = F(a); Fb = F(b)
        if Fa*Fb > 0.0:
            # Coarse scan to recover a valid bracket without moving near Tg
            xs = np.linspace(a, b, 33)
            Fs = [F(x) for x in xs]
            brkt = None
            for i in range(len(xs)-1):
                if Fs[i]*Fs[i+1] <= 0.0:
                    brkt = (xs[i], xs[i+1]); break
            if brkt is None:
                raise ValueError(f"No sign change in F(Twi) on [{a:.3f},{b:.3f}]. "
                                 f"F(a)={Fa:.3e}, F(b)={Fb:.3e}")
            a, b = brkt
        Twi_sol = Q_(brentq(F, a, b), "K")
        # Recover flux and Two
        h_r = self.gas.radiation_coefficient(Twi_sol).to("W/(m^2*K)")
        qprime = 2*pi*ri * (h_g + h_r) * (Tg - Twi_sol)
        Two = Twi_sol - qprime*log(ro/ri)/(2*pi*kw)
        return qprime.to("W/m")

    def rhs(self, x, y):
        Tg, pg, hw = y
        TgQ = Q_(Tg, "K")
        pgQ = Q_(pg, "Pa")
        hwQ = Q_(hw, "J/kg")
     
        # sync gas state
        self.gas.temperature = TgQ
        self.gas.pressure = pgQ
        self.water.enthalpy = hwQ

        # water bulk temperature from (P, h)
        PwQ = self.water.pressure
        w = IAPWS97(P=PwQ.to("megapascal").magnitude, h=hwQ.to("kJ/kg").magnitude)
        TcQ = Q_(w.T, "K")

        # local wall solve → heat flux per unit length
        qprime = self._wall_flux(TgQ, TcQ)

        # balances
        rho_g = self.gas.density.to("kg/m^3").magnitude
        cp_g  = self.gas.specific_heat.to("J/(kg*K)").magnitude
        D     = self.stage.hot_side.inner_diameter.to("m").magnitude
        A     = (pi*D**2/4)
        m_g   = self.gas.mass_flow_rate.to("kg/s").magnitude
        P     = pi*D
        v_g   = m_g/(rho_g*A)
        fD    = self.gas.friction_factor

        m_w   = self.water.mass_flow_rate.to("kg/s").magnitude

        dTgdx = - qprime.to("W/m").magnitude / (m_g*cp_g)
        dpdx  = - fD * rho_g * v_g**2 / (2*D)
        dhwdx =   qprime.to("W/m").magnitude /  m_w
        return [dTgdx, dpdx, dhwdx]

    def solve(self):
        L = self.stage.hot_side.inner_length.to("m").magnitude
        y0 = [
            self.gas.temperature.to("K").magnitude,
            self.gas.pressure.to("Pa").magnitude,
            self.water._h.to("J/kg").magnitude,
        ]
        return solve_ivp(self.rhs, (0,L), y0=y0)