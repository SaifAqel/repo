from dataclasses import dataclass

from config.schemas import Pass

@dataclass
class PassGeomCalc:
    """Container for a Pass object to extend with geometry calculations."""
    p: Pass