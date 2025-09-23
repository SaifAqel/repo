# geometry/tube.py
from dataclasses import dataclass
from math import pi
from .base import IGeometryStage
from common.units import ureg, Q_

@dataclass
class SingleTubeGeom(IGeometryStage):
    L: Q_            # tube length [m]
    Di: Q_           # inner diameter [m]
    Do: Q_           # outer diameter [m]
    roughness: Q_    # surface roughness [m]
    thickness: Q_    # wall thickness [m]
    R_fg: Q_         # gas-side fouling resistance [m²*K/W]
    R_fw: Q_         # water-side fouling resistance [m²*K/W]
    # gas inside; water outside

    def length(self) -> Q_:
        return self.L.to("m")

    def hydraulic_diameter(self) -> Q_:
        return self.Di.to("m")

    def flow_area(self) -> Q_:
        return (pi * self.Di**2 / 4).to("m^2")

    def heat_transfer_area(self) -> Q_:
        # inner area, consistent with gas-side h_i and TubeBank convention
        return (pi * self.Di * self.L).to("m^2")
