# flow/pressure_drop.py
from typing import List, Tuple
from common.units import Q_
from .friction import IFrictionModel
from .minor_losses import K_total

def darcy_dp(f: Q_, L: Q_, Dh: Q_, rho: Q_, u: Q_) -> Q_:
    f  = f.to("dimensionless")
    L  = L.to("m")
    Dh = Dh.to("m")
    rho= rho.to("kg/m^3")
    u  = u.to("m/s")
    return f * (L/Dh) * (rho * u**2 / 2)

class StagePressureDrop:
    def __init__(self, fm: IFrictionModel):
        self.fm = fm

    def total(self,
              Re: Q_, rel_rough: Q_,
              L: Q_, Dh: Q_,
              rho: Q_, u: Q_,
              elements: List[Tuple[str, Q_]]) -> Q_:
        f = self.fm.f(Re, rel_rough).to("dimensionless")
        dp_fric = darcy_dp(f, L, Dh, rho, u)
        K = K_total(elements).to("dimensionless")
        dp_minor = (rho.to("kg/m^3") * u.to("m/s")**2 / 2) * K
        return dp_fric + dp_minor
