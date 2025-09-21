# transport/htc_water_boil.py
from abc import ABC, abstractmethod
class IHTCWaterBoilingModel(ABC):
    def h(self, G: float, x: float, qpp: float, P: float, Dh: float) -> float:
        ...
class ChenBoilingModel(IHTCWaterBoilingModel):
    def h(self,G,x,qpp,P,Dh):
        Csf = 0.013
        F = 1.0
        hfc = 0.023*(G*Dh)**0.8*(1.0)**0.4/(Dh)
        hnb = 0.00122*(abs(qpp)**0.79)*(P**0.45)
        return F*hfc + S(hfc,qpp,P)*hnb
def S(hfc,qpp,P):
    return 1.0/(1.0+2.53e-6*hfc**1.17*abs(qpp)**0.7*P**-1.02)
class ThomModel(IHTCWaterBoilingModel):
    def h(self,G,x,qpp,P,Dh):
        return 55.0*(P/1e6)**0.12*(abs(qpp)/1e4)**0.67
class GungorWintertonModel(IHTCWaterBoilingModel):
    def h(self,G,x,qpp,P,Dh):
        Re_lo = G*Dh/1e-3
        Pr_lo = 1.0
        h_lo = 0.023*(Re_lo**0.8)*(Pr_lo**0.4)/Dh
        Xtt = ((1.0-x)/max(x,1e-6))**0.9
        F = 1.0+2400.0/x**1.16 if x>0 else 1.0
        return h_lo*(1.0+3000.0*Xtt**-1.16)+F*(abs(qpp)**0.79)*(P**0.45)*1e-6
