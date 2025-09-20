# transport/radiation.py
from math import pi
def linearized_h_rad(eps_g: float, sigma: float, Tg: float, Tw: float) -> float:
    return eps_g*4.0*sigma*((Tg+Tw)/2.0)**3
def mean_beam_length(volume: float, area: float) -> float:
    return 4.0*volume/area
