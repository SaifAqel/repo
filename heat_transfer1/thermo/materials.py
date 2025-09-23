# thermo/materials.py
from typing import Callable
from common.units import ureg, Q_

def k_wall_linear(a: float, b: float) -> Callable:
    """
    Returns a thermal conductivity function k(T) = a + b*T
    with pint quantities.
    a : W/(m*K)
    b : W/(m*K^2)
    """
    a = Q_(a, ureg.watt / (ureg.meter * ureg.kelvin))
    b = Q_(b, ureg.watt / (ureg.meter * ureg.kelvin**2))

    def k(T):
        # T must be Quantity in Kelvin
        T = Q_(T, ureg.kelvin)
        return a + b * T

    return k
