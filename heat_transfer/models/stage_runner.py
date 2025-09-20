# models/stage_runner.py
from dataclasses import dataclass
from typing import Dict
from .axial_marcher import AxialMarcher
from ..core.types import StageResult
from ..geometry.base import IGeometryStage

@dataclass
class StageRunner:
    geom: IGeometryStage
    marcher: AxialMarcher
    def run(self, gas_in, Tsat, G_water):
        cells = self.marcher.march(gas_in,Tsat,G_water)
        Q = sum(c.qpp*(self.geom.heat_transfer_area()/self.marcher.N) for c in cells)
        dP = gas_in.P - cells[-1].Pg
        T_out = cells[-1].Tg
        P_out = cells[-1].Pg
        Tw_min = min(c.Twi for c in cells)
        Tw_max = max(c.Two for c in cells)
        h_g_avg = sum(c.h_g for c in cells)/len(cells)
        eps_g_avg = sum(c.eps_g for c in cells)/len(cells)
        return StageResult(cells,Q,dP,T_out,P_out,Tw_min,Tw_max,h_g_avg,eps_g_avg)
