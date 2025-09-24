from dataclasses import dataclass
from heat_transfer.config.schemas import Pass
from common.units import ureg, Q_
import math

@dataclass
class PassGeomCalc:
    p: Pass

    @property
    def flow_area(self) -> Q_:
        d = self.p.geometry.inner_diameter
        return math.pi * (d/2)**2
    
