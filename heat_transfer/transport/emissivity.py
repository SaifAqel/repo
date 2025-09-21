# transport/emissivity.py
from abc import ABC, abstractmethod
from typing import Dict
from math import exp
class IEpsilonModel(ABC):
    @abstractmethod
    def epsilon(self, T: float, P: float, y: Dict[str,float], Lb: float) -> float:
        ...
class LecknerModel(IEpsilonModel):
    def epsilon(self, T,P,y,Lb):
        pco2 = y.get("CO2",0.0)*P
        ph2o = y.get("H2O",0.0)*P
        eb_co2 = 1.0-exp(-0.06*abs((pco2*Lb)**0.65))
        eb_h2o = 1.0-exp(-0.19*abs((ph2o*Lb)**0.45))
        return 1.0-(1.0-eb_co2)*(1.0-eb_h2o)
class SmithModel(IEpsilonModel):
    def epsilon(self,T,P,y,Lb):
        pco2 = y.get("CO2",0.0)*P
        ph2o = y.get("H2O",0.0)*P
        eb = 1.0-exp(-0.1*abs(((pco2+ph2o)*Lb)**0.75))
        return eb
class WSGGM_SimpleModel(IEpsilonModel):
    def epsilon(self,T,P,y,Lb):
        p = (y.get("CO2",0.0)+y.get("H2O",0.0))*P
        a = 1.0-exp(-0.3*abs((p*Lb)**0.5))
        return a
