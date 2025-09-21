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
from math import pi, isfinite
import logging
logging.basicConfig(level=logging.DEBUG)
import numpy as np

def _limit_dp(dx, Dh, rho, u, f, p, max_frac=0.2):
    """Limit Darcy drop per cell to a fraction of local pressure."""
    dp = f * (dx / Dh) * 0.5 * rho * u * u
    while dp > max_frac * p and dx > 1e-6:
        dx *= 0.5
        dp = f * (dx / Dh) * 0.5 * rho * u * u
    return dp, dx


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
        dx_base = self.L / self.N

        # per-path area from geometry; total paths from geometry (default 1)
        A_single = self.geom.flow_area()
        N_paths  = max(getattr(self.geom, "N_paths", 1), 1)
        A_eff    = A_single * N_paths             # total flow area used for u and Re
        Dh       = self.geom.hydraulic_diameter()
        A_ht     = self.geom.heat_transfer_area() / self.N

        x_ax = 0.0
        Tg, Pg, y = gas_in.T, gas_in.P, gas_in.y
        props = self.thermo.gas_props(Tg, Pg, y)
        rho, mu, k, cp, Pr = props["rho"], props["mu"], props["k"], props["cp"], props["Pr"]

        print(f"stage start: A_single={A_single}, N_paths={N_paths}, A_eff={A_eff}, Dh={Dh}, Pg_in={Pg}, Tg_in={Tg}")

        cells = []
        hb = LocalHeatBalance(lambda T: 1.0)
        Twi, Two = Tsat + 1.0, Tg - 1.0

        for i in range(self.N):
            dx = dx_base

            # use TOTAL area with TOTAL m_dot
            u  = gas_in.m_dot / (rho * A_eff)
            Re = (gas_in.m_dot / A_eff) * Dh / mu  # == rho*u*Dh/mu

            h_g  = self.htc_gas.h(Re, Pr, k, Dh, Tg, Two)
            eps  = self.eps_model.epsilon(Tg, Pg, y, Dh)
            h_rad = linearized_h_rad(eps, sigma, Tg, Two)
            h_w   = self.htc_water.h(G_water, 0.0, 1e4, Pg, Dh)

            qpp, Twi, Two = hb.solve(
                h_g, h_rad, h_w, self.R_fg, self.R_fw,
                self.geom.Di if hasattr(self.geom, "Di") else Dh,
                self.geom.Do if hasattr(self.geom, "Do") else Dh,
                Tg, Tsat, Twi, Two
            )

            Qcell = qpp * A_ht
            Tg    = Tg - Qcell / (gas_in.m_dot * cp)

            f = self.friction.f(Re, 0.0)

            # compressible Darcy step integrated over dx:
            mdot  = gas_in.m_dot
            T_avg = Tg
            C = (f / Dh) * 0.5 * (mdot**2) * self.R_fg / (A_eff**2)

            Pg_old = Pg
            Pg = (Pg*Pg - 2.0 * C * T_avg * dx)**0.5

            print(f"dp_terms: f={f}, dx={dx}, Dh={Dh}, rho={rho}, u2={u*u}, "
                f"P_old^2={Pg_old*Pg_old}, term={2.0*C*T_avg*dx}, P_new={Pg}")


            props = self.thermo.gas_props(Tg, Pg, y)
            rho, mu, k, cp, Pr = props["rho"], props["mu"], props["k"], props["cp"], props["Pr"]
            print(f"after dp: Pg={Pg}")
            print(f"props_ok? rho>0={rho>0 if np.isfinite(rho) else None}, mu={mu}, k={k}, cp={cp}, Pr={Pr}")

            cells.append(CellResult(x_ax, Tg, Pg, rho, u, cp, mu, k, Pr, h_g, eps, h_rad, h_w, qpp, Twi, Two))
            x_ax += dx_base

        return cells