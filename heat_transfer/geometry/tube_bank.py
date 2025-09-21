# geometry/tube_bank.py
from dataclasses import dataclass
from .base import IGeometryStage
from math import pi

@dataclass
class TubeBankGeom(IGeometryStage):
    L_m: float           # tube length
    Di_m: float          # tube inner diameter
    Do_m: float          # tube outer diameter
    Ntubes: int
    pitch_m: float
    arrangement: str
    roughness_m: float
    thickness_m: float
    frontal_width_m: float
    frontal_height_m: float
    # gas flows INSIDE; water is shell-side

    def length(self) -> float:
        return self.L_m

    def hydraulic_diameter(self) -> float:
        # single circular conduit -> Dh = Di
        return self.Di_m

    def flow_area(self) -> float:
        # total internal flow area of all tubes
        return self.Ntubes * (pi * (self.Di_m**2) / 4.0)

    def heat_transfer_area(self) -> float:
        # choose the basis you use for U consistently.
        # For gas-side h_i correlations: inner area
        return self.Ntubes * (pi * self.Di_m * self.L_m)

    # Optional: external shell-side (water) geometry, if you need it later
    def external_frontal_area(self) -> float:
        return self.frontal_width_m * self.frontal_height_m

    def external_blockage_area(self) -> float:
        return self.Ntubes * (pi * (self.Do_m**2) / 4.0)

    def external_crossflow_min_area(self) -> float:
        # only meaningful for crossflow outside tubes
        return self.external_frontal_area() - self.external_blockage_area()
