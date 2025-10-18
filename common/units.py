from __future__ import annotations  # at top of every module
import pint

ureg = pint.UnitRegistry()
Q_ = ureg.Quantity

class Converter:

    @staticmethod
    def _dim(q: Q_):  return q.to("dimensionless")

    @staticmethod
    def _MPa(q: Q_):  return q.to("MPa")

    @staticmethod
    def _K(q: Q_):    return q.to("K")

    @staticmethod
    def _kJkg(q: Q_): return q.to("kJ/kg")