# services/plant_runner.py
import yaml
from dataclasses import dataclass
from typing import List
from ..io.case_loader import load_case
from ..core.types import PlantResult, GasState

from common.units import ureg, Q_
from common.numcheck import assert_real_finite_q
@dataclass
class MultiStageHEX:
    settings_path: str

    def run(self) -> PlantResult:
        cfg, builders = load_case(self.settings_path)
        gas = builders["gas_inlet"]()
        assert_real_finite_q(gas.T, "gas_in.T")
        assert_real_finite_q(gas.P, "gas_in.P")
        for k, v in gas.y.items():
            assert_real_finite_q(v, f"gas_in.y[{k}]")
        stages = []
        Q_total = Q_(0.0, ureg.watt)
        dP_total = Q_(0.0, ureg.pascal)

        for s in builders["stages"]:
            res = s["runner"].run(
                gas,
                s["Tsat"],
                s["G_water"],
                s["Dh_water"],
                s["x_water_profile"],
                s["qpp_init_profile"]
            )
            stages.append(res)
            Q_total += res.Q
            dP_total += res.dP

            gas = GasState(
                res.T_out,                  # Kelvin (Quantity)
                res.P_out,                  # Pascal (Quantity)
                gas.m_dot,                  # kg/s
                gas.y,                      # composition (dimensionless)
                res.cells[-1].rhog,         # kg/m^3 (Quantity)
                res.cells[-1].cp,           # J/(kg*K) (Quantity)
                res.cells[-1].mu,           # Pa*s (Quantity)
                res.cells[-1].k,            # W/(m*K) (Quantity)
                res.cells[-1].Pr,           # dimensionless (Quantity)
                res.cells[-1].ug            # m/s (Quantity)
            )

        return PlantResult(stages, Q_total, dP_total, gas.T, gas.P)
