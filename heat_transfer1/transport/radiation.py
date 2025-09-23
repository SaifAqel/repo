# transport/radiation.py
from math import pi
from common.units import ureg, Q_

def linearized_h_rad(eps_g: float, sigma: float, Tg: float, Tw: float) -> float:
    eps_g = Q_(eps_g, ureg.dimensionless)
    sigma = Q_(sigma, ureg.watt / (ureg.meter**2 * ureg.kelvin**4))
    Tg    = Q_(Tg, ureg.kelvin)
    Tw    = Q_(Tw, ureg.kelvin)

    h_rad = eps_g * Q_(4.0, ureg.dimensionless) * sigma * (((Tg + Tw) / Q_(2.0, ureg.dimensionless)) ** 3)
    return h_rad.to(ureg.watt / (ureg.meter**2 * ureg.kelvin))

def mean_beam_length(volume: float, area: float) -> float:
    volume = Q_(volume, ureg.meter**3)
    area   = Q_(area, ureg.meter**2)
    mbl = Q_(4.0, ureg.dimensionless) * volume / area
    return mbl.to(ureg.meter)
