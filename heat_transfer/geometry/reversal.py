# geometry/reversal.py
from dataclasses import dataclass
from .base import IGeometryStage
@dataclass
class ReversalChamberGeom(IGeometryStage):
    L_eq: float
    Dh: float
    area_flow: float
    area_ht: float
    def length(self) -> float:
        return self.L_eq
    def hydraulic_diameter(self) -> float:
        return self.Dh
    def flow_area(self) -> float:
        return self.area_flow
    def heat_transfer_area(self) -> float:
        return self.area_ht*self.L_eq
