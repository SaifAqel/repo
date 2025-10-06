from __future__ import annotations
from dataclasses import dataclass
from typing import Dict
from common.units import ureg, Q_
from math import exp, log, sqrt, log10
from heat_transfer.fluid_props.WaterProps import WaterProps
from heat_transfer.fluid_props.GasProps import GasProps

@dataclass
class GasStream:
    mass_flow_rate: Q_
    temperature: Q_
    pressure: Q_
    composition: Dict[str, Q_]
    spectroscopic_data: Dict[str, Q_]

@dataclass
class Water:
    mass_flow_rate: Q_
    inlet_temperature: Q_
    outlet_temperature: Q_
    pressure: Q_
    superheat: Q_
    mu_exp: Q_
    C_sf: Q_
    n: Q_
    composition: Dict[str, Q_]



@dataclass
class GasStreamWithCalc:

    gas_props: object
    gas_stream: object
    geometry: object

    @property
    def density(self) -> Q_:
        return self.gas_props.density(self.gas_stream.temperature, self.gas_stream.pressure, self.gas_stream.composition)
    
    @property
    def specific_heat(self) -> Q_:
        return self.gas_props.cp(self.gas_stream.temperature, self.gas_stream.pressure, self.gas_stream.composition)
    
    @property
    def dynamic_viscosity(self) -> Q_:
        return self.gas_props.viscosity(self.gas_stream.temperature, self.gas_stream.pressure, self.gas_stream.composition)
    
    @property
    def thermal_conductivity(self) -> Q_:
        return self.gas_props.thermal_conductivity(self.gas_stream.temperature, self.gas_stream.pressure, self.gas_stream.composition)
    
    @property
    def mass_flux(self) -> Q_:
        return self.gas_stream.mass_flow_rate / self.geometry.tube_inner_flow_area

    @property
    def velocity(self) -> Q_:
        return self.gas_stream.mass_flow_rate / (self.density * self.geometry.tube_inner_flow_area * self.geometry.geometry.number_of_tubes)

    @property
    def reynolds_number(self) -> Q_:
        return (self.density * self.velocity * self.geometry.geometry.inner_diameter) / self.dynamic_viscosity
    
    @property
    def prandt_number(self) -> Q_:
        return self.dynamic_viscosity * self.specific_heat / self.thermal_conductivity
    
    @property
    def nusselt_number(self) -> Q_:
        n = 0.3
        return 0.023 * (self.reynolds_number ** 0.8) * (self.prandt_number ** n)
    
    @property
    def convective_coefficient(self) -> Q_:
        return self.nusselt_number * self.thermal_conductivity / self.geometry.geometry.inner_diameter
    
    @property
    def absorption_coefficient (self) -> Q_:
        return sum(self.gas_stream.composition[s] * self.gas_stream.spectroscopic_data[s] for s in self.gas_stream.composition)

    @property
    def emmissivity(self) -> Q_:
        return 1.0 - exp(-self.absorption_coefficient * self.geometry.path_length)
        
    def radiation_coefficient(self, T_wall) -> Q_:
        sigma = 5.670374419e-8 * ureg.watt / (ureg.meter**2 * ureg.kelvin**4)
        mean_temperature = (self.gas_stream.temperature + T_wall) / 2
        return 4.0 * sigma * (mean_temperature**3) * self.emmissivity

    @property
    def friction_factor(self) -> float:
        f = 0.02 
        Re = self.reynolds_number.magnitude
        rel_rough = self.geometry.rel_roughness
        tol = 1e-6
        max_iter = 50

        for _ in range(max_iter):
            inner = rel_rough/3.7 + 2.51/(Re * sqrt(f))
            g = 1.0/sqrt(f) + 2.0*log10(inner)
            dg = (-0.5)*f**(-1.5) + (2.0/log(10))*(1.0/inner)*(2.51/Re)*(-0.5)*f**(-1.5)
            f = f - g/dg
            if abs(g) < tol:
                break
        else:
            raise RuntimeError("Colebrook did not converge")

        return f

@dataclass
class WaterWithCalc:
    Water: object
    WaterProps: object

    @property
    def saturation_temperature(self) -> Q_:
        return self.WaterProps.Tsat(self.Water.pressure)
    
    @property
    def density_liquid(self) -> Q_:
        return self.WaterProps.rho_l(self.saturation_temperature)
    
    @property
    def density_vapor(self) -> Q_:
        return self.WaterProps.rho_v(self.saturation_temperature)
    
    @property
    def dynamic_viscosity_liquid(self) -> Q_:
        return self.WaterProps.mu_l(self.saturation_temperature)
    
    @property
    def latent_heat_of_vaporization(self) -> Q_:
        return self.WaterProps.h_fg(self.Water.pressure)
    
    @property
    def surface_tension(self) -> Q_:
        return self.WaterProps.sigma(self.saturation_temperature)
    
    @property
    def specific_capacity_liquid(self) -> Q_:
        return self.WaterProps.cp_l(self.saturation_temperature)
    
    @property
    def thermal_conductivity(self) -> Q_:
        return self.WaterProps.k_l(self.saturation_temperature)

    @property
    def prandt_number(self) -> Q_:
        return self.dynamic_viscosity_liquid * self.specific_capacity_liquid / self.thermal_conductivity
    
    def rohsenow_h(self, T_wall: Q_) -> Q_:
        g = 9.81 * ureg.m / (ureg.s**2)
        delta_T = T_wall - self.saturation_temperature # No boiling if ΔT <= 0; handle as single-phase if needed

        q_dot = (
            self.dynamic_viscosity_liquid
            * self.latent_heat_of_vaporization
            * ((g * (self.density_liquid - self.density_vapor) / self.surface_tension)**0.5)
            * ((self.specific_capacity_liquid * delta_T) / (self.Water.C_sf * self.latent_heat_of_vaporization * self.prandt_number**self.Water.n))**3
        )
        h = q_dot / delta_T
        return h.to("W / (m^2 * K)")
