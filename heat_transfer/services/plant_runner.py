# services/plant_runner.py
import yaml
from dataclasses import dataclass
from typing import List
from ..io.case_loader import load_case
from ..core.types import PlantResult, GasState
@dataclass
class MultiStageHEX:
    settings_path: str
    case_path: str
    def run(self) -> PlantResult:
        settings, case, builders = load_case(self.settings_path,self.case_path)
        gas = builders["gas_inlet"]()
        stages = []
        Q_total = 0.0
        dP_total = 0.0
        for s in builders["stages"]:
            res = s["runner"].run(gas,s["Tsat"],s["G_water"])
            stages.append(res)
            Q_total += res.Q
            dP_total += res.dP
            gas = GasState(res.T_out,res.P_out,gas.m_dot,gas.y,res.cells[-1].rhog,res.cells[-1].cp,res.cells[-1].mu,res.cells[-1].k,res.cells[-1].Pr,res.cells[-1].ug)
        return PlantResult(stages,Q_total,dP_total,gas.T,gas.P)
