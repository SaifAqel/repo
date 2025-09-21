# geometry/tube.py
from dataclasses import dataclass
from .base import IGeometryStage
from math import pi

@dataclass
class SingleTubeGeom(IGeometryStage):
    L_m: float
    Di_m: float
    Do_m: float
    roughness_m: float
    thickness_m: float
    R_fg_m2K_W: float
    R_fw_m2K_W: float
    # gas inside; water outside

    def length(self) -> float:
        return self.L_m

    def hydraulic_diameter(self) -> float:
        return self.Di_m

    def flow_area(self) -> float:
        return pi * (self.Di_m**2) / 4.0

    def heat_transfer_area(self) -> float:
        # inner area (consistent with gas-side h_i and with TubeBank)
        return pi * self.Di_m * self.L_m
