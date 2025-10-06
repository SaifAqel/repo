from math import log, pi
from dataclasses import dataclass
from common.units import ureg, Q_
from scipy.optimize import fsolve

@dataclass
class HeatSystem:
    geom: object
    gas: object
    water: object

    def water_resistance(self, T_wall: Q_) -> Q_:
        h_w = self.water.rohsenow_h(T_wall)
        di = self.geom.geometry.inner_diameter
        r_i = di / 2
        r_o = self.geom.outer_diameter / 2
        return (r_i / r_o) * (1 / h_w).to("m**2 * K / W")

    def gas_resistance(self, T_wall: Q_) -> Q_:
        """Convective + radiative resistance from gas to wall"""
        h_g = self.gas.convective_coefficient + self.gas.radiation_coefficient(T_wall)
        return (1 / h_g).to("K*m**2/W")

    @property
    def wall_resistance(self) -> Q_:
        di = self.geom.geometry.inner_diameter
        thickness = self.geom.geometry.wall.thickness
        r_i = di / 2
        r_o = r_i + thickness
        ln_term = log(r_o / r_i)
        return (r_i * ln_term / self.geom.geometry.wall.conductivity)

    @property
    def fouling_resistance_inner(self) -> Q_:
        return (self.geom.surfaces.inner.fouling_thickness / self.geom.surfaces.inner.fouling_conductivity).to("K*m**2/W")

    @property
    def fouling_resistance_outer(self) -> Q_:
        thickness = self.geom.surfaces.outer.fouling_thickness
        k_foul = self.geom.surfaces.outer.fouling_conductivity
        di = self.geom.geometry.inner_diameter
        r_i = di / 2
        r_o = self.geom.outer_diameter / 2
        return (r_i / r_o) * (thickness / k_foul).to("m**2 * K / W")

    def total_resistance(self, T_wall: Q_) -> Q_:
        return (
            self.water_resistance(T_wall)
            + self.wall_resistance
            + self.gas_resistance(T_wall)
            + self.fouling_resistance_inner
            + self.fouling_resistance_outer
        ).to("K*m**2/W")
    
    @property
    def fixed_resistance(self) -> Q_:
        return self.fouling_resistance_inner + self.wall_resistance + self.fouling_resistance_outer
    def wall_temperature(self) -> Q_:
        T_bulk = self.gas.gas_stream.temperature
        T_sat = self.water.saturation_temperature
        R_fixed = self.fixed_resistance

        def residual(T_w_mag: float) -> float:
            T_w = Q_(T_w_mag, 'kelvin')  # Gas-side wall temperature
            h_g = self.gas.convective_coefficient + self.gas.radiation_coefficient(T_w)
            q = h_g * (T_bulk - T_w)  # Heat flux from gas
            T_w_water = T_w - q * R_fixed  # Water-side wall temperature, accounting for fixed resistances
            h_w = self.water.rohsenow_h(T_w_water)
            q_water = h_w * (T_w_water - T_sat)  # Heat flux to water
            return (q - q_water).magnitude  # Residual in W/m²

        T_guess = ((T_bulk + T_sat) / 2).magnitude
        T_w_mag = fsolve(residual, T_guess)[0]
        return Q_(T_w_mag, 'kelvin')
    def heat_flux(self) -> Q_:
        T_w = self.wall_temperature()
        h_g = self.gas.convective_coefficient + self.gas.radiation_coefficient(T_w)
        return h_g * (self.gas.gas_stream.temperature - T_w)
    @property
    def q_(self) -> Q_:
        """Total heat rate (W)"""
        return self.heat_flux() * self.geom.tube_inner_perimeter
