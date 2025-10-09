from dataclasses import dataclass
from math import pi
from scipy.integrate import solve_ivp
from common.units import Q_, ureg
from heat_transfer.config.models import GasStream, WaterStream
from heat_transfer.functions.UA import UA
from iapws import IAPWS97

@dataclass
class HeatStageSolver:
    geom: object
    gas: GasStream
    water: WaterStream

    def rhs(self, x, y):
        Tg, pg, hw = y
        TgQ, pgQ, hwQ = Q_(Tg, "K"), Q_(pg, "Pa"), Q_(hw, "J/kg")

        # sync stream state for property-based model
        self.gas.temperature = TgQ
        self.gas.pressure = pgQ

        PwQ = self.water.pressure
        w = IAPWS97(P=PwQ.to("megapascal").magnitude, h=hwQ.to("kJ/kg").magnitude)
        Tc_star = Q_(w.T, "K").to("K").magnitude

        rho_g = self.gas.density.to("kg/m^3").magnitude        # property, no call
        cp_g  = self.gas.specific_heat.to("J/(kg*K)").magnitude # property, no call
        fD    = self.gas.friction_factor                        # property

        D = self.geom.geometry.inner_diameter.to("m").magnitude
        A = (pi * D**2 / 4)
        P = (pi * D)
        m_g = self.gas.mass_flow_rate.to("kg/s").magnitude
        m_w = self.water.mass_flow_rate.to("kg/s").magnitude
        v_g = m_g / (rho_g * A)

        U = UA(self).UA

        dTgdx = -(U * P / (m_g * cp_g)) * (Tg - Tc_star)
        dpdx  = - fD * rho_g * v_g**2 / (2 * D)
        dhwdx =  (U * P /  m_w)        * (Tg - Tc_star)
        return [dTgdx, dpdx, dhwdx]

    def solve(self):
        L = self.geom.geometry.inner_length.to("m").magnitude
        y0 = [
            self.gas.temperature.to("K").magnitude,
            self.gas.pressure.to("Pa").magnitude,
            self.water.enthalpy.to("J/kg").magnitude,
        ]
        return solve_ivp(self.rhs, [0, L], y0)
