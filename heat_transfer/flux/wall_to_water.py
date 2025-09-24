import math
from typing import Optional
from water_props.water_props import WaterProps

class PoolBoilingHT:
    def __init__(self, fluid: WaterProps, C_sf: float, n: float):
        self.f = fluid
        self.C_sf = C_sf
        self.n = n

    def qpp_rohsenow(self, T_wall: float, T_sat: Optional[float] = None) -> float:
        T_sat = self.f.T_sat if T_sat is None else T_sat
        dT = max(T_wall - T_sat, 0.0)
        if dT <= 0.0:
            return 0.0
        f = self.f
        g_term = math.sqrt(9.80665 * (f.rho_l - f.rho_v) / f.sigma)
        bracket = (f.cp_l * dT) / (self.C_sf * f.h_fg * (f.Pr_l ** self.n))
        return f.mu_l * f.h_fg * g_term * (bracket ** 3)

    def chf_zuber(self) -> float:
        f = self.f
        return 0.131 * f.h_fg * math.sqrt(f.rho_v) * (f.sigma * 9.80665 * (f.rho_l - f.rho_v))**0.25

    def check_boiling_limits(model, T_wall):
        qpp = model.qpp_rohsenow(T_wall)
        chf = model.chf_zuber()
        if qpp >= chf:
            print(f"Warning: q'' = {qpp:.3e} W/m² exceeds CHF = {chf:.3e} W/m². Nucleate boiling model invalid.")
        else:
            print(f"q'' = {qpp:.3e} W/m² (below CHF = {chf:.3e} W/m²).")
        return qpp, chf

    def h_boiling(self, T_wall: float, T_sat: Optional[float] = None) -> float:
        T_sat = self.f.T_sat if T_sat is None else T_sat
        dT = T_wall - T_sat
        if dT <= 0.0:
            return 0.0
        return self.qpp_rohsenow(T_wall, T_sat) / dT