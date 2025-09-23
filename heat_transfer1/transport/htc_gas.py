# transport/htc_gas.py
from abc import ABC, abstractmethod
import numpy as np
from common.units import Q_

class IHTCGasModel(ABC):
    @abstractmethod
    def h(self, Re: Q_, Pr: Q_, k: Q_, Dh: Q_, T_bulk: Q_, T_wall: Q_) -> Q_: ...

class GnielinskiModel(IHTCGasModel):
    def h(self, Re, Pr, k, Dh, T_bulk, T_wall):
        Re_m = Re.to("dimensionless").magnitude
        Pr_m = Pr.to("dimensionless").magnitude
        k_m  = k.to("W/m/K").magnitude
        Dh_m = Dh.to("m").magnitude
        if Re_m < 2300:
            Nu = 3.66  # simple laminar placeholder
        else:
            f = (0.79*np.log(Re_m) - 1.64)**-2
            Nu = (f/8)*(Re_m-1000)*Pr_m / (1 + 12.7*np.sqrt(f/8)*(Pr_m**(2/3)-1))
        h = Nu * k_m / Dh_m
        return Q_(h, "W/m^2/K")

class ChurchillBernsteinModel(IHTCGasModel):
    def h(self, Re, Pr, k, Dh, T_bulk, T_wall):
        Re_m = Re.to("dimensionless").magnitude
        Pr_m = Pr.to("dimensionless").magnitude
        k_m  = k.to("W/m/K").magnitude
        Dh_m = Dh.to("m").magnitude
        term = 0.3 + (0.62*Re_m**0.5*Pr_m**(1/3)) / (1+(0.4/Pr_m)**(2/3))**0.25 * (1+(Re_m/282000)**(5/8))**(4/5)
        h = term * k_m / Dh_m
        return Q_(h, "W/m^2/K")

class DittusBoelterModel(IHTCGasModel):
    """Turbulent internal flow. Valid roughly for Re > 1e4, 0.7 < Pr < 160."""
    def h(self, Re, Pr, k, Dh, T_bulk, T_wall):
        Re_m = Re.to("dimensionless").magnitude
        Pr_m = Pr.to("dimensionless").magnitude
        k_m  = k.to("W/m/K").magnitude
        Dh_m = Dh.to("m").magnitude

        # Choose exponent based on heating/cooling of the fluid
        Tb = T_bulk.to("K").magnitude
        Tw = T_wall.to("K").magnitude
        n = 0.4 if Tw > Tb else 0.3

        # Simple laminar fallback (keeps parity with Gnielinski placeholder)
        if Re_m < 2300:
            Nu = 3.66
        else:
            Nu = 0.023 * (Re_m**0.8) * (Pr_m**n)

        h = Nu * k_m / Dh_m
        return Q_(h, "W/m^2/K")
