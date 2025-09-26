from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, Any, Mapping, Optional
import fnmatch
import tomllib
from common.units import ureg, Q_
from heat_transfer.fluid_props.GasProps import HotFlueGas
from heat_transfer.calc_ops.stage_with_calc import PassWithCalc
class Stream:
    def __init__(self, mass_flow_rate, HotFlueGas, PassWithCalc):
        self.m_dot = mass_flow_rate
        self.rho = HotFlueGas.density
        self.PassWithCalc = PassWithCalc
        inlet_temperature
        inlet_pressure
        composition

    def velocity(self) -> Q_:
        return self.m_dot / (self.rho * self.PassWithCalc.tube_inner_flow_area)
