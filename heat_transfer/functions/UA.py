from math import log
from common.units import ureg, Q_
from heat_transfer.models.streams import GasStream, WaterStream
from dataclasses import dataclass

@dataclass
class UA:
    stage = object
    gas = GasStream
    water = WaterStream
    
    @property
    def R_gas(self) -> Q_:
        h = self.gas.convective_coefficient + self.gas.radiation_coefficient
        return 1/h
    
    @property
    def R_foul_gas(self) -> Q_:
        return self.geom.wall.surfaces.inner.fouling_thickness / self.geom.wall.surfaces.inner.fouling_conductivity

    @property
    def R_wall(self) -> Q_:
        r_i = self.stage.geom.inner_diameter / 2
        r_o = self.stage.geom.outer_diameter / 2
        k = self.stage.wall.conductivity
        return r_i * log(r_o/r_i) / k

    @property
    def R_foul_water(self) -> Q_:
        return 0

    @property
    def R_water(self) -> Q_:
        h = self.water.convective_coefficient
        return 1/h

    @property
    def R_total(self) -> Q_:
        return (
                self.R_gas(self)
                + self.R_foul_gas(self)
                + self.R_wall(self)
                + self.R_foul_water(self)
                + self.R_water(self)
        )

    @property
    def U_overall(self) -> Q_:
        return 1 / self.R_total(self)

    @property
    def UA(self) -> Q_:
        A = self.stage.geom.HX_area
        return self.U_overall(self) * A
