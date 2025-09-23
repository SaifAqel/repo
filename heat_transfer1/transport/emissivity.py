# transport/emissivity.py
from abc import ABC, abstractmethod
from typing import Any, Dict, Sequence, Tuple
import numpy as np
from common.units import Q_, ureg

SPECIES_ORDER: Tuple[str, ...] = ("CO2", "H2O", "O2", "N2")

def _yfraction(y: Any):
    if isinstance(y, dict):
        def get(name: str) -> float:
            v = y.get(name, 0.0)
            return float(v.to("dimensionless").magnitude) if hasattr(v, "to") else float(v)
        return get
    if isinstance(y, tuple) and len(y)==2 and isinstance(y[0], np.ndarray) and isinstance(y[1], Sequence):
        arr, names = y
        lookup = {str(n): float(arr[i]) for i, n in enumerate(names)}
        return lambda name: float(lookup.get(name, 0.0))
    if isinstance(y, np.ndarray):
        arr = np.asarray(y, dtype=float).ravel()
        lookup = {n: float(arr[i]) if i < arr.size else 0.0 for i, n in enumerate(SPECIES_ORDER)}
        return lambda name: float(lookup.get(name, 0.0))
    raise TypeError("Unsupported y")

class IEpsilonModel(ABC):
    @abstractmethod
    def epsilon(self, T: Q_, P: Q_, y: Dict[str, Q_], Lb: Q_) -> Q_: ...

class SimpleGrayGas(IEpsilonModel):
    def epsilon(self, T, P, y, Lb):
        yget = _yfraction(y)
        # crude weighted gray gas for stability
        eps = 0.25*yget("CO2") + 0.75*yget("H2O")
        eps = max(0.0, min(1.0, eps))
        return Q_(eps, "dimensionless")

class LecknerModel(IEpsilonModel):
    """Leckner correlation for CO2–H2O gas emissivity."""

    def epsilon(self, T, P, y, Lb):
        yget = _yfraction(y)
        p_co2 = (P.to("atm").magnitude) * yget("CO2")
        p_h2o = (P.to("atm").magnitude) * yget("H2O")
        # Leckner’s chart is piecewise; below is a simple curve-fit form
        eps = 1 - np.exp(-0.65 * (p_co2 + p_h2o) ** 0.45 * (Lb.to("m").magnitude) ** 0.65)
        return Q_(max(0.0, min(1.0, eps)), "dimensionless")


class SmithModel(IEpsilonModel):
    """Smith, Shen, and Friedman (SSF) model for gas emissivity."""

    def epsilon(self, T, P, y, Lb):
        yget = _yfraction(y)
        p_co2 = (P.to("atm").magnitude) * yget("CO2")
        p_h2o = (P.to("atm").magnitude) * yget("H2O")
        L = Lb.to("m").magnitude

        # Typical SSF effective pressure term
        p_eff = (p_co2 ** 0.5 + p_h2o ** 0.5) ** 2
        eps = 1 - np.exp(-1.5 * (p_eff * L) ** 0.65)
        return Q_(max(0.0, min(1.0, eps)), "dimensionless")


class WSGGM_SimpleModel(IEpsilonModel):
    """Weighted-Sum-of-Gray-Gases (WSGGM) simplified model."""

    def epsilon(self, T, P, y, Lb):
        yget = _yfraction(y)
        Tc = T.to("K").magnitude
        L = Lb.to("m").magnitude
        p_co2 = (P.to("atm").magnitude) * yget("CO2")
        p_h2o = (P.to("atm").magnitude) * yget("H2O")

        # Simple 3-gray-gas WSGGM
        a = [0.2, 0.3, 0.5]
        k = [0.1, 1.0, 10.0]  # absorption coefficients [1/m] (example values)
        p_eff = p_co2 + p_h2o

        eps = 0.0
        for ai, ki in zip(a, k):
            eps += ai * (1 - np.exp(-ki * p_eff * L * (Tc / 1000) ** 0.5))

        return Q_(max(0.0, min(1.0, eps)), "dimensionless")
