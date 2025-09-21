# geometry/reversal.py
from dataclasses import dataclass
from .base import IGeometryStage

@dataclass
class ReversalChamberGeom(IGeometryStage):
    L_eq_m: float                  # equivalent straight length
    area_flow_m2: float             # cross-sectional flow area A
    wetted_perimeter_m: float      # wetted perimeter P for Dh = 4A/P
    area_ht_per_length_m: float    # local heat-transfer perimeter (P_ht)

    # Note: chamber is not circular; avoid Di/Do fields

    def length(self) -> float:
        return self.L_eq_m

    def hydraulic_diameter(self) -> float:
        return 4.0 * self.area_flow_m2/ self.wetted_perimeter_m

    def flow_area(self) -> float:
        return self.area_flow_m2

    def heat_transfer_area(self) -> float:
        # heat-transfer area = perimeter in contact with fluid × length
        return self.area_ht_per_length_m * self.L_eq_m
