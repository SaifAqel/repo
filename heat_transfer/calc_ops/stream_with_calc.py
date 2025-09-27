from __future__ import annotations
from dataclasses import dataclass
from typing import Dict
from common.units import ureg, Q_

@dataclass
class GasStream:
    mass_flow_rate: Q_
    inlet_temperature: Q_
    inlet_pressure: Q_
    composition: Dict[str, Q_]


@dataclass
class GasStreamWithCalc(GasStream):
    def __init__(self, mass_flow_rate, HotFlueGas, PassWithCalc, correlation_exponent):
        self.m_dot = mass_flow_rate
        self.rho = HotFlueGas.density
        self.mu = HotFlueGas.viscosity
        self.cp = HotFlueGas.cp
        self.k = HotFlueGas.thermal_conductivity
        self.di = PassWithCalc.geometry.inner_diameter
        self.A = PassWithCalc.tube_inner_flow_area
        self.nt = PassWithCalc.geometry.number_of_tubes
        self.n = correlation_exponent

    @property
    def velocity(self) -> Q_:
        return self.m_dot / (self.rho * self.PassWithCalc.tube_inner_flow_area * self.nt)

    @property
    def reynolds_number(self) -> Q_:
        return (self.rho * self.velocity * self.di) / self.mu
    
    @property
    def prandt_number(self) -> Q_:
        return self.mu * self.cp / self.k
    
    @property
    def nusselt_number(self) -> Q_:
        return 0.023 * (self.reynolds_number ** 0.8) * (self.prandt_number ** self.n)
    
    @property
    def Convective_coefficient(self) -> Q_:
        return self.nusselt_number * self.k / self.di