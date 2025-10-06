from math import log, pi
from dataclasses import dataclass
from common.units import ureg, Q_

@dataclass
class HeatSystem:
    geom: object
    gas: object
    water: object

    def water_resistance(self, T_wall: Q_) -> Q_:
        """Convective resistance from wall to water (boiling/condensing)"""
        h_w = self.water.rohsenow_h(T_wall)
        return (1 / h_w).to("K*m**2/W")

    def gas_resistance(self, T_wall: Q_) -> Q_:
        """Convective + radiative resistance from gas to wall"""
        h_g = self.gas.convective_coefficient + self.gas.radiation_coefficient(T_wall)
        return (1 / h_g).to("K*m**2/W")

    @property
    def wall_resistance(self) -> Q_:
        """Conduction through tube wall"""
        R = log((self.geom.r_o / self.geom.r_i).to_base_units().magnitude)
        return Q_(R, "") / (2 * pi * self.geom.k * self.geom.r_i).to("W/K")

    @property
    def fouling_resistance_inner(self) -> Q_:
        return (self.geom.thickness_inner / self.geom.k_foul_inner).to("K*m**2/W")

    @property
    def fouling_resistance_outer(self) -> Q_:
        return (self.geom.thickness_outer / self.geom.k_foul_outer).to("K*m**2/W")

    def total_resistance(self, T_wall: Q_) -> Q_:
        return (
            self.water_resistance(T_wall)
            + self.wall_resistance
            + self.gas_resistance(T_wall)
            + self.fouling_resistance_inner
            + self.fouling_resistance_outer
        ).to("K*m**2/W")

    def wall_temperature(self) -> Q_:
        """
        Solve for T_wall implicitly from energy balance:
        (T_bulk - T_wall)/R_gas = (T_wall - T_water)/R_water
        """
        R_g = self.gas_resistance  # function
        R_w = self.water_resistance  # function
        T_bulk = self.gas.gas_stream.temperature
        T_water = self.water.saturation_temperature

        # Solve explicitly:
        Rg = R_g(T_bulk)  # approximate at current step
        Rw = R_w(T_water)  # approximate at current step
        T_wall = (T_bulk * Rw + T_water * Rg) / (Rg + Rw)
        return T_wall

    def heat_flux(self) -> Q_:
        """Compute heat flux per unit area"""
        T_w = self.wall_temperature()
        return (self.gas.gas_stream.temperature - self.wall_temperature()) / self.gas_resistance(T_w)

    @property
    def q_(self) -> Q_:
        """Total heat rate (W)"""
        return self.heat_flux() * self.geom.tube_inner_perimeter
