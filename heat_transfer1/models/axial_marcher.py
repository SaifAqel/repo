# models/axial_marcher.py
from dataclasses import dataclass
from typing import List, Tuple, Optional, Dict
from common.units import Q_, ureg
from ..geometry.base import IGeometryStage
from ..thermo.provider import IThermoProvider
from ..transport.htc_gas import IHTCGasModel
from ..transport.emissivity import IEpsilonModel
from ..transport.htc_water_boil import IHTCWaterBoilingModel
from ..flow.friction import IFrictionModel
from ..core.types import GasState, CellResult
from ..transport.radiation import linearized_h_rad
from common.numcheck import assert_real_finite_q
SIGMA = Q_(5.670374419e-8, "W/m^2/K^4")

@dataclass
class AxialMarcher:
    geom: IGeometryStage
    thermo: IThermoProvider
    htc_gas: IHTCGasModel
    eps_model: IEpsilonModel
    htc_water: IHTCWaterBoilingModel
    friction: IFrictionModel
    L: Q_
    N: int
    flow: str
    R_fg: Q_
    R_fw: Q_
    losses: List[Tuple[str, Q_]]

    def march(self,
            gas_in: GasState,
            Tsat: Q_,
            G_water: Q_,
            Dh_water: Optional[Q_],
            x_water_profile: Optional[List[Q_]],
            qpp_init_profile: Optional[List[Q_]]) -> List[CellResult]:

        L = self.geom.length().to("m")
        Dh = self.geom.hydraulic_diameter().to("m")
        A_flow = self.geom.flow_area().to("m^2")
        A_ht = self.geom.heat_transfer_area().to("m^2")
        dx = L / Q_(self.N, "1")
        A_cell = A_ht * (dx / L)

        # Optional: fail early if provided profiles are too short
        if x_water_profile is not None and len(x_water_profile) < self.N:
            raise ValueError(f"x_water_profile length {len(x_water_profile)} < N {self.N}")
        if qpp_init_profile is not None and len(qpp_init_profile) < self.N:
            raise ValueError(f"qpp_init_profile length {len(qpp_init_profile)} < N {self.N}")

        cells: List[CellResult] = []
        Tg = gas_in.T.to("K")
        Pg = gas_in.P.to("Pa")
        m_dot = gas_in.m_dot.to("kg/s")
        y = gas_in.y

        for i in range(self.N):
            x_pos = dx * Q_(i + 1, "1")

            # properties
            props = self.thermo.gas_props(Tg, Pg, y)
            rho = props["rho"]; cp = props["cp"]; mu = props["mu"]; k = props["k"]; Pr = props["Pr"]
            u = (m_dot / (rho * A_flow)).to("m/s")
            Re = (rho * u * Dh / mu).to("dimensionless")

            # film coefficients
            h_g = self.htc_gas.h(Re, Pr, k, Dh, Tg, Tg)
            eps_g = self.eps_model.epsilon(Tg, Pg, y, (A_flow * L / A_ht))
            h_rad = linearized_h_rad(eps_g, SIGMA, Tg, Tsat)
            h_g_total = (h_g + h_rad).to("W/m^2/K")

            # water-side inputs with safe indexing
            if x_water_profile is not None:
                xq = x_water_profile[i]  # will not IndexError due to early check
            else:
                xq = Q_(0.0, "1")

            Dh_w = Dh_water if Dh_water is not None else Dh

            if qpp_init_profile is not None:
                qpp_guess = qpp_init_profile[i]
            else:
                qpp_guess = Q_(1e4, "W/m^2")

            h_w = self.htc_water.h(G_water, xq, qpp_guess, Pg, Dh_w)

            # overall U
            R_tot = (1 / h_g_total + self.R_fg + 1 / h_w + self.R_fw).to("m^2*K/W")
            U = (1 / R_tot).to("W/m^2/K")

            # heat flux and state update
            DT = (Tg - Tsat).to("K")
            qpp = (U * DT).to("W/m^2")
            Qcell = (qpp * A_cell).to("W")
            dT = (Qcell / (m_dot * cp)).to("K")
            Tg_new = (Tg - dT).to("K")

            # pressure drop per cell
            rel_rough = Q_(1e-5, "1")
            dp = (self.friction.f(Re, rel_rough) * (dx / Dh) * (rho * u**2 / 2)).to("Pa")
            Pg_new = (Pg - dp).to("Pa")

            cells.append(CellResult(
                x=x_pos.to("m"),
                Tg=Tg_new, Pg=Pg_new,
                rhog=rho, cp=cp, mu=mu, k=k, Pr=Pr, ug=u,
                Twi=(Tsat + Q_(0.1, "K")),
                Two=(Tsat + Q_(0.0, "K")),
                qpp=qpp, h_g=h_g_total, eps_g=eps_g
            ))

            # advance marching variables
            Tg = Tg_new
            Pg = Pg_new

            # sanity checks on current state (do NOT reset to inlet)
            assert_real_finite_q(Tg, "Tg")
            assert_real_finite_q(Pg, "Pg")
            for k_sp, v_sp in y.items():
                assert_real_finite_q(v_sp, f"y[{k_sp}]")

        return cells

