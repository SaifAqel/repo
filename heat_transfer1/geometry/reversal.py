# geometry/reversal.py
from dataclasses import dataclass
from .base import IGeometryStage
from common.units import ureg, Q_

@dataclass
class ReversalChamberGeom(IGeometryStage):
    L_eq: Q_                 # equivalent straight length [m]
    area_flow: Q_            # flow cross-sectional area [m²]
    wetted_perimeter: Q_     # wetted perimeter [m]
    ht_perimeter: Q_         # local heat-transfer perimeter per length [m]

    def length(self) -> Q_:
        return self.L_eq.to("m")

    def hydraulic_diameter(self) -> Q_:
        return (4.0 * self.area_flow / self.wetted_perimeter).to("m")

    def flow_area(self) -> Q_:
        return self.area_flow.to("m^2")

    def heat_transfer_area(self) -> Q_:
        # perimeter in contact with fluid × length
        return (self.ht_perimeter * self.L_eq).to("m^2")
