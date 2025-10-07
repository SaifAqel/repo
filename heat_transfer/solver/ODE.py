from scipy.integrate import solve_ivp
from heat_transfer.thermal.heat import HeatSystem
from common.units import ureg, Q_
# heat_transfer/solver/ODE.py  (replace FireTubeGasODE implementation)

from scipy.integrate import solve_ivp
from heat_transfer.thermal.heat import HeatSystem
from common.units import ureg, Q_
from heat_transfer.solver.water_balance import WaterStateConverter, WaterEnergyRHS

class FireTubeGasODE:
    """
    Coupled ODE for gas T, gas P and water specific enthalpy h_w.
    State y = [T_g (K), P (Pa), h_w (J/kg)].
    """

    def __init__(self, pass_with_calc, gas_stream_with_calc, water_with_calc, water_converter: WaterStateConverter = None, water_rhs: WaterEnergyRHS = None):
        self.geom = pass_with_calc
        self.gas = gas_stream_with_calc
        self.water = water_with_calc
        self.heat_system = HeatSystem(geom=pass_with_calc,
                                      gas=gas_stream_with_calc,
                                      water=water_with_calc)
        self.water_converter = water_converter or WaterStateConverter(self.water)
        self.water_rhs = water_rhs or WaterEnergyRHS(self.water, pass_with_calc)

    def rhs(self, z, y):
        T_bulk = Q_(y[0], 'kelvin')
        P = Q_(y[1], 'pascal')
        h_w = Q_(y[2], 'J/kg')

        # update local gas state
        self.gas.gas_stream.temperature = T_bulk
        self.gas.gas_stream.pressure = P

        # convert water enthalpy -> local temperature and expose to HeatSystem
        T_w_local = self.water_converter.T_from_h(h_w, self.water.water_stream.pressure)
        # place local value into the Water dataclass so HeatSystem reads it if needed
        self.water.water_stream.temperature = T_w_local

        # solve wall temperature and get local heat flux (W/m^2)
        T_wall = self.heat_system.wall_temperature()
        q_pp = self.heat_system.heat_flux()  # W / m^2

        # heat per unit length (W/m)
        q_per_length = (q_pp * self.geom.tube_inner_perimeter).to("W / m")

        # gas energy equation: dTg/dz = - q' / (m_dot_gas * cp)
        dTdz = (- q_per_length) / (self.gas.gas_stream.mass_flow_rate * self.gas.specific_heat)

        # pressure drop (same form as before)
        f = self.gas.friction_factor
        rho = self.gas.density.magnitude
        v = self.gas.velocity.magnitude
        di = self.geom.geometry.inner_diameter.magnitude
        dpdz = - (f / di) * (0.5 * rho * v**2)

        # water energy: dh/dz = q' / m_dot_per_circuit
        dhw_dz = self.water_rhs.dhdz(h_w, q_per_length)

        return [dTdz.to('kelvin / meter').magnitude, dpdz, dhw_dz.to('J / (kg * m)').magnitude]

    def solve(self, z_span, y0, method="BDF", **kwargs):
        return solve_ivp(self.rhs, z_span, y0, method=method, **kwargs)
