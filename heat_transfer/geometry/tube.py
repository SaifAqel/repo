# geometry/tube.py
from dataclasses import dataclass
from .base import IGeometryStage
@dataclass
class SingleTubeGeom(IGeometryStage):
    L: float
    Di: float
    Do: float
    roughness: float
    thickness: float
    R_fg: float
    R_fw: float
    def length(self) -> float:
        return self.L
    def hydraulic_diameter(self) -> float:
        return self.Di
    def flow_area(self) -> float:
        from math import pi
        return pi*(self.Di**2)/4.0
    def heat_transfer_area(self) -> float:
        from math import pi
        return pi*self.Di*self.L
