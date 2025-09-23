from dataclasses import dataclass

from config.schemas import Reversal

@dataclass
class ReversalChamberGeomCalc:
    """Container for a ReversalChamber to extend with calculations later."""
    rc: Reversal