from scipy.integrate import solve_ivp
from heat_transfer.thermal.heat import HeatSystem
from common.units import ureg, Q_

class FireTubeGasODE:
    def __init__(self, pass_with_calc, gas_stream_with_calc, water_with_calc):
        self.geom = pass_with_calc
        self.gas = gas_stream_with_calc
        self.water = water_with_calc
        self.heat_system = HeatSystem(geom=pass_with_calc,
                                      gas=gas_stream_with_calc,
                                      water=water_with_calc)

    def rhs(self, z, y):
        T_bulk, P = y
        self.gas.gas_stream.temperature = T_bulk  # update bulk temp in stream
        q_dot = self.heat_system.q_
        dTdz = - q_dot / (self.gas.gas_stream.mass_flow_rate * self.gas.specific_heat)
        # dpdz can be added if needed
        dpdz = 0
        return [dTdz.magnitude, dpdz]

    def solve(self, z_span, y0, method="BDF", **kwargs):
        return solve_ivp(self.rhs, z_span, y0, method=method, **kwargs)
