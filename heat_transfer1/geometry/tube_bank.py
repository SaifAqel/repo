# geometry/tube_bank.py
from dataclasses import dataclass
from math import pi
from .base import IGeometryStage
from common.units import Q_

@dataclass
class TubeBankGeom(IGeometryStage):
    L: Q_
    Di: Q_
    Do: Q_
    Ntubes: int
    pitch: Q_
    arrangement: str
    roughness: Q_
    thickness: Q_
    frontal_width: Q_
    frontal_height: Q_

    def length(self) -> Q_:
        return self.L.to("m")

    def hydraulic_diameter(self) -> Q_:
        return self.Di.to("m")

    def flow_area(self) -> Q_:
        return (self.Ntubes * (pi * self.Di**2 / 4)).to("m^2")

    def heat_transfer_area(self) -> Q_:
        return (self.Ntubes * (pi * self.Di * self.L)).to("m^2")

    # external helpers
    def external_frontal_area(self) -> Q_:
        return (self.frontal_width * self.frontal_height).to("m^2")

    def external_blockage_area(self) -> Q_:
        return (self.Ntubes * (pi * self.Do**2 / 4)).to("m^2")

    def external_crossflow_min_area(self) -> Q_:
        return (self.external_frontal_area() - self.external_blockage_area()).to("m^2")
