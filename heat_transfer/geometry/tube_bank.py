# geometry/tube_bank.py
from dataclasses import dataclass
from .base import IGeometryStage
from math import pi
@dataclass
class TubeBankGeom(IGeometryStage):
    L: float
    Di: float
    Do: float
    Ntubes: int
    pitch: float
    arrangement: str
    roughness: float
    thickness: float
    frontal_width: float
    frontal_height: float
    def length(self) -> float:
        return self.L
    def hydraulic_diameter(self) -> float:
        return self.Di
    def flow_area(self) -> float:
        return self.frontal_width*self.frontal_height - self.Ntubes*pi*(self.Do**2)/4.0
    def heat_transfer_area(self) -> float:
        return self.Ntubes*pi*self.Do*self.L
