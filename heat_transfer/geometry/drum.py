from dataclasses import dataclass

from config.schemas import Drum

@dataclass
class DrumGeomCalc:
    """Container for a Drum to extend with calculations later."""
    drum: Drum