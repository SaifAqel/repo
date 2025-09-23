# thermo/water.py
from dataclasses import dataclass
from typing import Dict
from common.units import ureg, Q_

@dataclass
class SaturationProps:
    Tsat: float   # K
    rho_l: float  # kg/m^3
    rho_v: float  # kg/m^3
    h_fg: float   # J/kg

def drum_inventory_update(m_w: float, dm: float) -> float:
    m_w = Q_(m_w, ureg.kg)
    dm  = Q_(dm,  ureg.kg)
    return (m_w + dm).to(ureg.kg)

def to_saturation(props: Dict[str, float]) -> SaturationProps:
    Tsat  = Q_(props["Tsat"], ureg.kelvin)
    rho_l = Q_(props["rho_l"], ureg.kg / ureg.meter**3)
    rho_v = Q_(props["rho_v"], ureg.kg / ureg.meter**3)
    h_fg  = Q_(props["h_fg"], ureg.joule / ureg.kg)
    return SaturationProps(Tsat, rho_l, rho_v, h_fg)
