from dataclasses import dataclass
from math import pi
from scipy.integrate import solve_ivp
from common.units import Q_, ureg
from heat_transfer.models.streams import GasStream, WaterStream
from heat_transfer.functions.UA import UA

@dataclass
class HeatStageSolver:
    geom: object
    gas: GasStream
    water: WaterStream

    def rhs(self, x, y):
        Tg, pg, hw = y
        TgQ, pgQ, hwQ = Q_(Tg, "K"), Q_(pg, "Pa"), Q_(hw, "J/kg")

        Tc_starQ, _ = self.water.map_state(hwQ)
        Tc_star = Tc_starQ.to("K").magnitude

        rho_g = self.gas.density(TgQ, pgQ).to("kg/m^3").magnitude
        cp_g  = self.gas.specific_heat(TgQ, pgQ).to("J/(kg*K)").magnitude
        fD    = self.gas.friction_factor(TgQ, pgQ)

        D = self.geom.geometry.inner_diameter.to("m").magnitude
        A = (pi * D**2 / 4)
        P = (pi * D)
        m_g = self.gas.mass_flow_rate.to("kg/s").magnitude
        m_w = self.water.mass_flow_rate.to("kg/s").magnitude
        v_g = m_g / (rho_g * A)

        U = UA.UA(self)

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
        sol = solve_ivp(self.rhs, [0, L], y0, dense_output=True)
        return sol
