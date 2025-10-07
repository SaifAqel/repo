from __future__ import annotations
from dataclasses import dataclass
from typing import Dict
from common.units import ureg, Q_
from math import exp, log, sqrt, log10
from heat_transfer.fluid_props.WaterProps import WaterProps
from heat_transfer.fluid_props.GasProps import GasProps
from heat_transfer.config.models import GasStream, Water

@dataclass
class GasStreamWithCalc:

    gas_props: GasProps
    gas_stream: GasStream
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
        return self.gas_stream.mass_flow_rate / (self.density * self.geometry.tube_inner_flow_area)

    @property
    def reynolds_number(self) -> Q_:
        return (self.density * self.velocity * self.geometry.geometry.inner_diameter) / self.dynamic_viscosity
    
    @property
    def prandt_number(self) -> Q_:
        return self.dynamic_viscosity * self.specific_heat / self.thermal_conductivity
    
    @property
    def nusselt_number(self) -> Q_:
        n = 0.3
        Re_mag = self.reynolds_number.magnitude
        Pr_mag = self.prandt_number.magnitude
        nu_mag = 0.023 * (Re_mag ** 0.8) * (Pr_mag ** n)
        return Q_(nu_mag, ureg.dimensionless)
    
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
    water_stream: Water
    water_props: WaterProps

    @property
    def saturation_temperature(self) -> Q_:
        return self.water_props.Tsat(self.water_stream.pressure)

    @property
    def phase(self) -> str:
        if self.water_stream.temperature < self.saturation_temperature:
            return "liquid"
        elif self.water_stream.temperature > self.saturation_temperature:
            return "vapor"
        else:
            return "saturated"
        
    @property
    def quality(self) -> Q_:
        if self.phase != "saturated":
            return Q_(0.0 if self.phase == "liquid" else 1.0, "dimensionless")
        
        h = self.water_props.h_l(self.water_stream.temperature)
        h_l = self.water_props.h_l(self.water_stream.temperature)
        h_v = self.water_props.h_v(self.water_stream.temperature)
        x = (h - h_l) / (h_v - h_l)
        return Q_(x, "dimensionless")

    @property
    def enthalpy(self) -> Q_:
        if self.phase == "liquid":
            return self.water_props.h_l(self.water_stream.temperature)
        elif self.phase == "vapor":
            return self.water_props.h_v(self.water_stream.temperature)
        else:  # saturated/mixed
            x = self.quality.magnitude
            h_l = self.water_props.h_l(self.saturation_temperature)
            h_v = self.water_props.h_v(self.saturation_temperature)
            return Q_(h_l.magnitude * (1-x) + h_v.magnitude * x, "J/kg")

    @property
    def density(self) -> Q_:
        if self.phase == "liquid":
            return self.water_props.rho_l(self.water_stream.temperature)
        elif self.phase == "vapor":
            return self.water_props.rho_v(self.water_stream.temperature)
        else:
            x = self.quality.magnitude
            rho_l = self.water_props.rho_l(self.saturation_temperature)
            rho_v = self.water_props.rho_v(self.saturation_temperature)
            return Q_(rho_l.magnitude * (1-x) + rho_v.magnitude * x, "kg/m^3")
        
    @property
    def specific_heat(self) -> Q_:
        """Specific heat capacity of the water at current temperature and phase."""
        if self.phase == "liquid":
            return self.water_props.cp_l(self.water_stream.temperature)
        elif self.phase == "vapor":
            return self.water_props.cp_v(self.water_stream.temperature)
        else:  # saturated/mixed
            cp_l = self.water_props.cp_l(self.saturation_temperature)
            cp_v = self.water_props.cp_v(self.saturation_temperature)
            x = self.quality.magnitude
            return Q_(cp_l.magnitude * (1 - x) + cp_v.magnitude * x, "J/(kg*K)")

    @property
    def thermal_conductivity(self) -> Q_:
        """Thermal conductivity of the water at current temperature and phase."""
        if self.phase == "liquid":
            return self.water_props.k_l(self.water_stream.temperature)
        elif self.phase == "vapor":
            return self.water_props.k_v(self.water_stream.temperature)
        else:  # saturated/mixed
            k_l = self.water_props.k_l(self.saturation_temperature)
            k_v = self.water_props.k_v(self.saturation_temperature)
            x = self.quality.magnitude
            return Q_(k_l.magnitude * (1 - x) + k_v.magnitude * x, "W/(m*K)")

    @property
    def dynamic_viscosity(self) -> Q_:
        """Dynamic viscosity of the water at current temperature and phase."""
        if self.phase == "liquid":
            return self.water_props.mu_l(self.water_stream.temperature)
        elif self.phase == "vapor":
            return self.water_props.mu_v(self.water_stream.temperature)
        else:  # saturated/mixed
            mu_l = self.water_props.mu_l(self.saturation_temperature)
            mu_v = self.water_props.mu_v(self.saturation_temperature)
            x = self.quality.magnitude
            return Q_(mu_l.magnitude * (1 - x) + mu_v.magnitude * x, "Pa*s")


    @property
    def surface_tension(self) -> Q_:
        """Surface tension of water at current temperature."""
        # Usually surface tension is for liquid; for saturated mixtures, use saturation temperature
        T = self.water_stream.temperature if self.phase == "liquid" else self.saturation_temperature
        return self.water_props.sigma(T)  # Units: N/m or kg/s^2

    @property
    def latent_heat_of_vaporization(self) -> Q_:
        """Latent heat of vaporization at current temperature."""
        # For liquid or vapor or saturated, use saturation temperature
        return self.water_props.h_v(self.saturation_temperature) - self.water_props.h_l(self.saturation_temperature)  # Units: J/kg
    

    
    
    @property
    def prandt_number(self) -> Q_:
        return self.dynamic_viscosity * self.specific_heat / self.thermal_conductivity
    
    def rohsenow_h(self, T_wall: Q_) -> Q_:
        g = 9.81 * ureg.m / (ureg.s**2)
        delta_T = T_wall - self.saturation_temperature
        n = 1
        Csf = 1

        q_dot = (
            self.dynamic_viscosity
            * self.latent_heat_of_vaporization
            * ((g * (self.density - self.density) / self.surface_tension)**0.5)
            * ((self.specific_heat * delta_T) / (Csf * self.latent_heat_of_vaporization * self.prandt_number**n))**3
        )
        h = q_dot / delta_T
        return h.to("W / (m^2 * K)")
