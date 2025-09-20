# core/units.py
from pint import UnitRegistry
ureg = UnitRegistry()
Q_ = ureg.Quantity
def ensure_quantity(x, units):
    q = Q_(x)
    return q.to(units)
