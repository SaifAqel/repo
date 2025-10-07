from math import log, pi
from dataclasses import dataclass
from common.units import ureg, Q_
from scipy.optimize import fsolve
from scipy.optimize import root_scalar


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
        """Per-area thermal resistance for cylindrical wall (m² K / W)"""
        r_i = self.geom.geometry.inner_diameter / 2
        r_o = r_i + self.geom.geometry.wall.thickness
        R = log(r_o / r_i)
        return (r_i * R / self.geom.geometry.wall.conductivity).to("m**2 * K / W")

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
        T_g = self.gas.gas_stream.temperature.to('kelvin').magnitude
        T_sat = self.water.saturation_temperature.to('kelvin').magnitude

        def residual(T_w):
            # Gas-side flux (W/m²), including fouling inner
            R_gas = self.gas_resistance(Q_(T_w, 'kelvin')) + self.fouling_resistance_inner
            q_g = (T_g - T_w) / R_gas.to('m**2 * K / W').magnitude  # magnitude for scalar

            # Water-side flux (W/m²)
            delta_T = T_w - T_sat
            if delta_T <= 0:  # No boiling if T_wall <= T_sat; assume small h_w (convection only)
                return q_g - 0  # Or implement single-phase h here
            q_w = self.water.rohsenow_h(Q_(T_w, 'kelvin')).to('W / (m**2 * K)').magnitude * delta_T

            # Include wall + outer fouling in balance if needed, but approx here
            return q_g - q_w

        # Solve with bracket; guess T_sat + 10K
        sol = root_scalar(residual, bracket=[T_sat + 1e-6, T_g], x0=T_sat + 10)
        if not sol.converged:
            raise ValueError("T_wall solver failed")
        return Q_(sol.root, 'kelvin')

    def heat_flux(self) -> Q_:
        T_w = self.wall_temperature()
        R_gas = self.gas_resistance(T_w) + self.fouling_resistance_inner
        return ((self.gas.gas_stream.temperature - T_w) / R_gas).to("W / m**2")

    @property
    def q_(self) -> Q_:
        return self.heat_flux() * self.geom.tube_inner_perimeter * self.geom.geometry.number_of_tubes
