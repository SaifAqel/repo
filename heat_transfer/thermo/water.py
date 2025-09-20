# thermo/water.py
from dataclasses import dataclass
from typing import Dict
@dataclass
class SaturationProps:
    Tsat: float
    rho_l: float
    rho_v: float
    h_fg: float
def drum_inventory_update(m_w: float, dm: float) -> float:
    return m_w + dm
def to_saturation(props: Dict[str,float]) -> SaturationProps:
    return SaturationProps(props["Tsat"],props["rho_l"],props["rho_v"],props["h_fg"])
