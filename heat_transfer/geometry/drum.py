from dataclasses import dataclass
from config.schemas import Drum
from common.units import ureg, Q_
import math

@dataclass
class DrumGeomCalc:
    """Container for a Drum to extend with calculations later."""
    drum: Drum

   