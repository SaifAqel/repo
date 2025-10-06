from scipy.integrate import solve_ivp
from heat_transfer.thermal.heat import HeatSystem
class FireTubeGasODE:

    def __init__(self, pass_with_calc, gas_stream_with_calc, water_with_calc):
        self.geom = pass_with_calc
        self.stream = gas_stream_with_calc
        self.heat_system = HeatSystem(geom = pass_with_calc,
                                       gas = gas_stream_with_calc,
                                       water = water_with_calc
                                       )

    def rhs(self):
        dTdz = - self.heat_system.q_ / (self.m_dot * self.cp)
        dpdz = - (self.fD / (2.0 * self.D_i)) * (self.G*self.G) / self.rho
        return [dTdz, dpdz]

    def solve(self, z_span, y0, method="BDF", **kwargs):
        return solve_ivp(self.rhs, z_span, y0, method=method, **kwargs)