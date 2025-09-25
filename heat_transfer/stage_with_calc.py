from heat_transfer.config.schemas import Reversal
from heat_transfer.config.schemas import Pass
from heat_transfer.config.schemas import Drum
from dataclasses import dataclass
from common.units import ureg, Q_
from math import pi

@dataclass
class PassWithCalc(Pass):

    @property
    def tube_inner_flow_area(self) -> Q_:
        di = self.geometry.inner_diameter
        return pi * (di/2)**2

    @property
    def tube_inner_perimeter(self) -> Q_:
        di = self.geometry.inner_diameter
        return pi * di
    
    @property
    def tube_outer_perimeter(self) -> Q_:
        do = self.outer_diameter
        return pi * do

    @property
    def tube_inner_heat_transfer_area(self) -> Q_:
        Pi = self.tube_inner_perimeter
        L = self.geometry.inner_length
        return Pi * L

    @property
    def outer_diameter(self) -> Q_:
        di = self.geometry.inner_diameter
        k = self.geometry.wall.thickness
        return di + (2 * k)

    @property
    def cross_section_outer_area(self) -> Q_:
        do = self.outer_diameter
        return pi * (do/2)**2

@dataclass
class DrumWithCalc(Drum):

    @property
    def cross_section_inner_area(self) -> Q_:
        di = self.geometry.inner_diameter
        return pi * (di/2)**2
    
    @property
    def outer_diameter(self) -> Q_:
        di = self.geometry.inner_diameter
        k = self.geometry.wall.thickness
        return di + (2 * k)

    @property
    def cross_section_outer_area(self) -> Q_:
        do = self.outer_diameter
        return pi * (do/2)**2