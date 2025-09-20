# flow/pressure_drop.py
from .friction import IFrictionModel
from .minor_losses import K_total
def darcy_dp(f: float, L: float, Dh: float, rho: float, u: float) -> float:
    return f*(L/Dh)*0.5*rho*u*u
class PressureDropCalculator:
    def __init__(self, friction_model: IFrictionModel):
        self.fm = friction_model
    def stage_dp(self, Re: float, rel_rough: float, L: float, Dh: float, rho: float, u: float, elements) -> float:
        f = self.fm.f(Re,rel_rough)
        dp_fric = darcy_dp(f,L,Dh,rho,u)
        dp_minor = 0.5*rho*u*u*K_total(elements)
        return dp_fric + dp_minor
