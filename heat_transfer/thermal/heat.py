from math import log, pi
from common.units import ureg, Q_
from dataclasses import dataclass

@dataclass
class HeatSystem:

    geom: object    
    gas: object
    water: object

    @property
    def water_resistance(self):
        return (1 / self.water.rohsenow_h).to("kelvin*meter**2/watt")
    
    @property
    def gas_resistance(self):
        return (1 / (self.gas.radiation_coefficient + self.gas.Convective_coefficient)).to("kelvin*meter**2/watt")
    
    @property
    def wall_resistance(self):
        value = log((self.geom.r_o / self.geom.r_i).to_base_units().magnitude)
        return Q_(value, "") / (2 * pi * self.geom.k * self.geom.r_i).to("watt/kelvin")
    
    @property
    def fouling_resistance_inner(self):
        return (self.geom.thickness_inner / self.geom.k_foul_inner).to("kelvin*meter**2/watt")
    
    @property
    def fouling_resistance_outer(self):
        return (self.geom.thickness_outer / self.geom.k_foul_outer).to("kelvin*meter**2/watt")
    
    @property
    def total_resistance(self):
        return (
            self.water_resistance
            + self.wall_resistance
            + self.gas_resistance
            + self.fouling_resistance_inner
            + self.fouling_resistance_outer
        ).to("kelvin*meter**2/watt")
    
    @property
    def delta_T(self):
        return self.gas.gas_stream.temperature - self.water.saturation_temperature

    @property
    def heat_flux(self):
        return self.delta_T / self.total_resistance
    @property
    def q_(self):
        return self.heat_flux * self.geom.tube_inner_perimeter
