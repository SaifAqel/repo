from heat_transfer.config.models import GasStream, WaterStream, FirePass, SmokePass, Reversal, Economiser
from heat_transfer.functions.htc_water import WaterHTC
from common.units import Q_
from math import pi, log

class HeatRate:
    def __init__(self, stage: FirePass | SmokePass | Reversal | Economiser, gas: GasStream, water: WaterStream):
        self.stage = stage
        self.gas = gas
        self.water = water

    def gas_resistance_per_length(self) -> Q_:
        h = self.gas.radiation_coefficient + self.gas.convective_coefficient
        return 1 / ( h * self.stage.hot_side.inner_perimeter )
    
    def inner_fouling_resistance_per_length(self) -> Q_:
        return self.stage.hot_side.wall.surfaces.inner.fouling_thickness / ( self.stage.hot_side.wall.surfaces.inner.fouling_conductivity * self.stage.hot_side.inner_perimeter )

    def wall_resistance_per_length(self) -> Q_:
        return log ( self.stage.hot_side.outer_diameter / self.stage.hot_side.inner_diameter ) / ( 2 * pi * self.stage.hot_side.wall.conductivity )

    def outer_fouling_resistance_per_length(self) -> Q_:
        return self.stage.hot_side.wall.surfaces.outer.fouling_thickness / ( self.stage.hot_side.wall.surfaces.outer.fouling_conductivity * self.stage.hot_side.outer_perimeter )
    
    def water_resistance_per_length(self) -> Q_:
        return 1 / ( self.stage.hot_side.outer_perimeter * WaterHTC(self).calc_htc() )

    def total_resistance_per_length(self) -> Q_:
        return ( 
            self.gas_resistance_per_length() 
            + self.inner_fouling_resistance_per_length() 
            + self.wall_resistance_per_length() 
            + self.outer_fouling_resistance_per_length() 
            + self.water_resistance_per_length() 
        )
    
    def heat_rate_per_length(self) -> Q_:
        return ( self.gas.temperature - self.water.temperature ) / self.total_resistance_per_length