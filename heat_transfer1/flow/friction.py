# flow/friction.py
from abc import ABC, abstractmethod
import numpy as np
from common.units import Q_

class IFrictionModel(ABC):
    @abstractmethod
    def f(self, Re: Q_, rel_rough: Q_) -> Q_: ...

class ColebrookWhite(IFrictionModel):
    def f(self, Re: Q_, rel_rough: Q_) -> Q_:
        Re_m = Re.to("dimensionless").magnitude
        rr_m = rel_rough.to("dimensionless").magnitude
        if Re_m <= 0:
            raise ValueError("Re must be > 0")
        # Start guess (Haaland)
        f = 1.0/(-1.8*np.log10((rr_m/3.7)**1.11 + 6.9/Re_m))**2
        for _ in range(50):
            g = 1/np.sqrt(f) + 2.0*np.log10(rr_m/3.7 + 2.51/(Re_m*np.sqrt(f)))
            df = 1e-8
            g2 = 1/np.sqrt(f+df) + 2.0*np.log10(rr_m/3.7 + 2.51/(Re_m*np.sqrt(f+df)))
            dg = (g2 - g)/df
            step = g/dg
            f_new = f - step
            if not np.isfinite(f_new) or f_new <= 1e-8:
                break
            if abs(step) < 1e-12:
                f = f_new
                break
            f = f_new
        return Q_(float(f), "dimensionless")

class Haaland(IFrictionModel):
    def f(self, Re: Q_, rel_rough: Q_) -> Q_:
        Re_m = Re.to("dimensionless").magnitude
        rr_m = rel_rough.to("dimensionless").magnitude
        return Q_(1.0/(-1.8*np.log10((rr_m/3.7)**1.11 + 6.9/Re_m))**2, "dimensionless")

class Churchill(IFrictionModel):
    """
    Churchill (1977) explicit correlation for Darcy friction factor.
    Valid across laminar, transition, and turbulent regimes.
    """
    def f(self, Re: Q_, rel_rough: Q_) -> Q_:
        Re_m = Re.to("dimensionless").magnitude
        rr_m = rel_rough.to("dimensionless").magnitude
        if Re_m <= 0:
            raise ValueError("Re must be > 0")
        # Natural logs per original paper
        A = (2.457 * np.log(1.0 / ((7.0 / Re_m)**0.9 + 0.27 * rr_m)))**16
        B = (37530.0 / Re_m)**16
        f = 8.0 * ((8.0 / Re_m)**12 + 1.0 / (A + B)**1.5)**(1.0 / 12.0)
        return Q_(float(f), "dimensionless")
