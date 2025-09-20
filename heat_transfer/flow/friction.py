# flow/friction.py
from abc import ABC, abstractmethod
from math import log10
class IFrictionModel(ABC):
    @abstractmethod
    def f(self, Re: float, rel_rough: float) -> float:
        ...
class ColebrookWhite(IFrictionModel):
    def f(self, Re, rel_rough):
        x = 1.0
        for _ in range(50):
            g = -2.0*log10(rel_rough/3.7 + 2.51/(Re*x)) - x
            dg = 2.0/(2.302585*(rel_rough/3.7 + 2.51/(Re*x))) * (2.51/(Re*x*x)) -1.0
            x = x - g/dg
        return x*x
class Churchill(IFrictionModel):
    def f(self, Re, rel_rough):
        A = (2.457*log10((7.0/Re)**0.9 + 0.27*rel_rough))**16
        B = (37530.0/Re)**16
        return 8.0*((8.0/Re)**12 + 1.0/(A+B)**1.5)**(1.0/12.0)
class Haaland(IFrictionModel):
    def f(self, Re, rel_rough):
        return 1.0/(-1.8*log10((rel_rough/3.7)**1.11 + 6.9/Re))**2
