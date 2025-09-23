# models/stage_runner.py
from dataclasses import dataclass
from common.units import Q_
from ..geometry.base import IGeometryStage
from ..core.types import StageResult
from .axial_marcher import AxialMarcher

@dataclass
class StageRunner:
    geom: IGeometryStage
    marcher: AxialMarcher

    def run(self, gas_in, Tsat, G_water, Dh_water, x_water_profile, qpp_init_profile):
        cells = self.marcher.march(gas_in, Tsat, G_water, Dh_water, x_water_profile, qpp_init_profile)
        Q = sum(c.qpp for c in cells) * (self.geom.heat_transfer_area()/self.geom.length()/Q_(self.marcher.N,"1"))
        dP = gas_in.P - cells[-1].Pg
        T_out = cells[-1].Tg
        P_out = cells[-1].Pg
        Tw_min = min(c.Twi for c in cells)
        Tw_max = max(c.Two for c in cells)
        h_g_avg = sum(c.h_g for c in cells) / Q_(len(cells), "1")
        eps_g_avg = sum(c.eps_g for c in cells) / Q_(len(cells), "1")
        return StageResult(cells, Q.to("W"), dP.to("Pa"), T_out.to("K"), P_out.to("Pa"), Tw_min.to("K"), Tw_max.to("K"), h_g_avg.to("W/m^2/K"), eps_g_avg.to("1"))
