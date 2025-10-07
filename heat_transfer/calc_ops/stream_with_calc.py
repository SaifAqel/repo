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
    water_stream: object  # expects attributes: temperature (Q_), pressure (Q_), mass_flow_rate (Q_)
    water_props: WaterProps

    @property
    def saturation_temperature(self) -> Q_:
        return self.water_props.Tsat(self.water_stream.pressure)

    @property
    def phase(self) -> str:
        Tsat = self.saturation_temperature
        if self.water_stream.temperature < Tsat: return "liquid"
        if self.water_stream.temperature > Tsat: return "vapor"
        return "saturated"

    @property
    def quality(self) -> Q_:
        P = self.water_stream.pressure
        if self.phase != "saturated":
            return Q_(0.0 if self.phase == "liquid" else 1.0, "dimensionless")

        h_l = self.water_props.h_l_sat(P)
        h_v = self.water_props.h_v_sat(P)
        # avoid division by zero
        x = 0.0
        if (h_v - h_l).magnitude != 0:
            x = 0.0
        return Q_(x, "dimensionless")

    @property
    def enthalpy(self) -> Q_:
        P = self.water_stream.pressure
        T = self.water_stream.temperature
        if self.phase == "liquid":
            return self.water_props.h_single(P, T)
        if self.phase == "vapor":
            return self.water_props.h_single(P, T)
        # saturated mixture: compute from saturated h_l and h_v with quality if available
        h_l = self.water_props.h_l_sat(P)
        h_v = self.water_props.h_v_sat(P)
        # fallback assume saturated liquid
        return Q_(h_l.magnitude, "J/kg")

    @property
    def density(self) -> Q_:
        P = self.water_stream.pressure
        T = self.water_stream.temperature
        if self.phase == "liquid":
            return self.water_props.rho_l(P, T)
        if self.phase == "vapor":
            return self.water_props.rho_v(P, T)
        # saturated mixture: linear interpolation by quality using saturation properties
        rho_l = self.water_props.rho_l_sat(P)
        rho_v = self.water_props.rho_v_sat(P)
        x = 0.0
        return Q_(rho_l.magnitude * (1 - x) + rho_v.magnitude * x, "kg/m^3")

    @property
    def specific_heat(self) -> Q_:
        P = self.water_stream.pressure
        T = self.water_stream.temperature
        if self.phase == "liquid":
            return self.water_props.cp_l(P, T)
        if self.phase == "vapor":
            return self.water_props.cp_v(P, T)
        cp_l = self.water_props.cp_l_sat(P)
        cp_v = self.water_props.cp_v_sat(P)
        x = 0.0
        return Q_(cp_l.magnitude * (1 - x) + cp_v.magnitude * x, "J/(kg*K)")

    @property
    def thermal_conductivity(self) -> Q_:
        P = self.water_stream.pressure
        T = self.water_stream.temperature
        if self.phase == "liquid":
            return self.water_props.k_l(P, T)
        if self.phase == "vapor":
            return self.water_props.k_v(P, T)
        k_l = self.water_props.k_l_sat(P)
        k_v = self.water_props.k_v_sat(P)
        x = 0.0
        return Q_(k_l.magnitude * (1 - x) + k_v.magnitude * x, "W/(m*K)")

    @property
    def dynamic_viscosity(self) -> Q_:
        P = self.water_stream.pressure
        T = self.water_stream.temperature
        if self.phase == "liquid":
            return self.water_props.mu_l(P, T)
        if self.phase == "vapor":
            return self.water_props.mu_v(P, T)
        mu_l = self.water_props.mu_l_sat(P)
        mu_v = self.water_props.mu_v_sat(P)
        x = 0.0
        return Q_(mu_l.magnitude * (1 - x) + mu_v.magnitude * x, "Pa*s")

    @property
    def surface_tension(self) -> Q_:
        P = self.water_stream.pressure
        # always use saturation T for surface tension when near boiling
        T_use = self.saturation_temperature if self.phase != "liquid" else self.water_stream.temperature
        return self.water_props.sigma(P, T_use)

    @property
    def latent_heat_of_vaporization(self) -> Q_:
        P = self.water_stream.pressure
        return self.water_props.h_fg(P)

    @property
    def prandt_number(self) -> Q_:
        return self.dynamic_viscosity * self.specific_heat / self.thermal_conductivity

    def rohsenow_h(self, T_wall: Q_) -> Q_:
        """Use saturated liquid/vapor properties at the fixed pressure for Rohsenow correlation."""
        P = self.water_stream.pressure
        Tsat = self.water_props.Tsat(P)
        delta_T = T_wall - Tsat
        if delta_T.magnitude <= 0:
            return Q_(0.0, "W/(m^2*K)")

        g = 9.81 * ureg.m / (ureg.s**2)
        n = 1
        Csf = 0.013

        # use saturated properties at P
        mu_l = self.water_props.mu_l_sat(P)
        h_fg = self.water_props.h_fg(P)
        rho_l = self.water_props.rho_l_sat(P)
        rho_v = self.water_props.rho_v_sat(P)
        sigma = self.water_props.sigma_sat(P)
        cp_l = self.water_props.cp_l_sat(P)
        k_l = self.water_props.k_l_sat(P)

        q_dot = (
            mu_l
            * h_fg
            * ((g * (rho_l - rho_v) / sigma) ** 0.5)
            * ((cp_l * delta_T) / (Csf * h_fg * (mu_l * cp_l / k_l) ** n)) ** 3
        )
        h = q_dot / delta_T
        return h.to("W / (m^2 * K)")
