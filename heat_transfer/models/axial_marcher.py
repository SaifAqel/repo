# models/axial_marcher.py
from ..geometry.base import IGeometryStage
from ..thermo.provider import IThermoProvider
from ..transport.htc_gas import IHTCGasModel
from ..transport.emissivity import IEpsilonModel
from ..transport.htc_water_boil import IHTCWaterBoilingModel
from ..flow.friction import IFrictionModel
from typing import List, Tuple
from dataclasses import dataclass
from typing import List, Dict
from .heat_balance import LocalHeatBalance
from ..core.types import CellResult, GasState
from ..transport.radiation import linearized_h_rad
from math import pi
sigma = 5.670374419e-8
@dataclass
class AxialMarcher:
    geom: IGeometryStage
    thermo: IThermoProvider
    htc_gas: IHTCGasModel
    eps_model: IEpsilonModel
    htc_water: IHTCWaterBoilingModel
    friction: IFrictionModel
    L: float
    N: int
    flow: str
    R_fg: float
    R_fw: float
    losses: List[Tuple[str, float]]
    def march(self, gas_in: GasState, Tsat: float, G_water: float) -> List[CellResult]:
        dx = self.L/self.N
        A_flow = self.geom.flow_area()
        A_ht = self.geom.heat_transfer_area()/self.N
        Dh = self.geom.hydraulic_diameter()
        x = 0.0
        Tg = gas_in.T
        Pg = gas_in.P
        y = gas_in.y
        props = self.thermo.gas_props(Tg,Pg,y)
        rho = props["rho"]
        mu = props["mu"]
        k = props["k"]
        cp = props["cp"]
        Pr = props["Pr"]
        u = gas_in.m_dot/(rho*A_flow)
        cells = []
        hb = LocalHeatBalance(lambda T: 1.0)
        Twi = Tsat+1.0
        Two = Tg-1.0
        for i in range(self.N):
            Re = rho*u*Dh/mu
            h_g = self.htc_gas.h(Re,Pr,k,Dh,Tg,Two)
            eps = self.eps_model.epsilon(Tg,Pg,y, Dh)
            h_rad = linearized_h_rad(eps,sigma,Tg,Two)
            h_w = self.htc_water.h(G_water,0.0,1e4,Pg,Dh)
            qpp,Twi,Two = hb.solve(h_g,h_rad,h_w,self.R_fg,self.R_fw,self.geom.Di if hasattr(self.geom,"Di") else Dh,self.geom.Do if hasattr(self.geom,"Do") else Dh,Tg,Tsat,Twi,Two)
            Qcell = qpp*A_ht
            Tg = Tg - Qcell/(gas_in.m_dot*cp)
            dp = self.friction.f(Re,0.0)*dx/Dh*0.5*rho*u*u
            Pg = Pg - dp
            props = self.thermo.gas_props(Tg,Pg,y)
            rho = props["rho"]
            mu = props["mu"]
            k = props["k"]
            cp = props["cp"]
            Pr = props["Pr"]
            u = gas_in.m_dot/(rho*A_flow)
            cells.append(CellResult(x,Tg,Pg,rho,u,cp,mu,k,Pr,h_g,eps,h_rad,h_w,qpp,Twi,Two))
            x += dx
        return cells
