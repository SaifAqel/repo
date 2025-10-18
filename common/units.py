import pint
ureg = pint.UnitRegistry()
Q_ = ureg.Quantity

class Converter:
    def _dim(q: Q_):  return q.to("dimensionless")
    def _MPa(q: Q_):  return q.to("MPa")
    def _K(q: Q_):    return q.to("K")
    def _kJkg(q: Q_): return q.to("kJ/kg")