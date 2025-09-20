# transport/htc_gas.py
from abc import ABC, abstractmethod
from typing import Dict
from math import log10
import numpy as np
class IHTCGasModel(ABC):
    @abstractmethod
    def h(self, Re: float, Pr: float, k: float, Dh: float, T_bulk: float, T_wall: float) -> float:
        ...
class GnielinskiModel(IHTCGasModel):
    def h(self, Re, Pr, k, Dh, T_bulk, T_wall):
        f = (0.79*np.log10(Re)-1.64)**-2
        Nu = (f/8.0)*(Re-1000.0)*Pr/(1.0+12.7*(f/8.0)**0.5*(Pr**(2.0/3.0)-1.0))
        return Nu*k/Dh
class DittusBoelterModel(IHTCGasModel):
    def h(self, Re, Pr, k, Dh, T_bulk, T_wall):
        n = 0.4 if T_wall<T_bulk else 0.3
        Nu = 0.023*(Re**0.8)*(Pr**n)
        return Nu*k/Dh
class ChurchillBernsteinModel(IHTCGasModel):
    def h(self, Re, Pr, k, Dh, T_bulk, T_wall):
        Nu = 0.3+0.62*(Re**0.5)*(Pr**(1.0/3.0))/((1.0+(0.4/Pr)**(2.0/3.0))**0.25)*(1.0+(Re/282000.0)**(5.0/8.0))**(4.0/5.0)
        return Nu*k/Dh
