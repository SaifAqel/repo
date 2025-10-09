from __future__ import annotations
from dataclasses import dataclass
from typing import Dict
from common.units import ureg, Q_
from math import exp, log, sqrt, log10
from heat_transfer.functions.fluid_props import WaterProps, GasProps
from heat_transfer.functions.htc_water import HTCFunctions

@dataclass(frozen=True)
class Surface:
    roughness: Q_
    emissivity: Q_
    fouling_thickness: Q_
    fouling_conductivity: Q_

@dataclass(frozen=True)
class Surfaces:
    inner: Surface
    outer: Surface

@dataclass(frozen=True)
class Wall:
    thickness: Q_
    conductivity: Q_
    surfaces: Surfaces

@dataclass(frozen=True)
class Nozzle:
    pass
@dataclass(frozen=True)
class Nozzles:
    inner: Nozzle
    outer: Nozzle

@dataclass(frozen=True)
class TubeGeometry:
    inner_diameter: Q_
    inner_length: Q_

@dataclass(frozen=True)
class ReversalChamberGeometry:
    pass

@dataclass(frozen=True)
class BankLayout:
    tubes_number: Q_
    shape: str
    pitch: Q_

@dataclass(frozen=True)
class ShellGeometry:
    flow_area: Q_
    wetted_perimeter: Q_

@dataclass(frozen=True)
class Shell:
    geometry: ShellGeometry
    wall: Wall



@dataclass(frozen=True)
class FirePass:
    geometry: TubeGeometry
    wall: Wall
    shell: Shell

@dataclass(frozen=True)
class SmokePass:
    geometry: TubeGeometry
    layout: BankLayout
    wall: Wall
    shell: Shell

@dataclass(frozen=True)
class Reversal:
    geometry: ReversalChamberGeometry
    nozzles: Nozzles
    wall: Wall
    shell: Shell

@dataclass(frozen=True)
class HotSide:
    area: Q_
    
@dataclass(frozen=True)
class ColdSide:
    area: Q_

@dataclass(frozen=True)
class Economiser:
    hot_side: HotSide
    cold_side: ColdSide


@dataclass(frozen=True)
class Stages:
    HX_1: FirePass
    HX_2: Reversal
    HX_3: SmokePass
    HX_4: Reversal
    HX_5: SmokePass
    HX_6: Economiser

    def __iter__(self):
        return iter((self.HX_1, self.HX_2, self.HX_3, self.HX_4, self.HX_5, self.HX_6))




@dataclass
class GasStream:
        mass_flow_rate: Q_
        temperature: Q_
        pressure: Q_
        composition: Dict[str, Q_]
        spectroscopic_data: Dict[str, Q_]

        stage: object
        gas_props: GasProps

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
class WaterStream:
        mass_flow_rate: Q_
        temperature: Q_
        pressure: Q_
        composition: Dict[str, Q_]

        stage: object
        water_props: WaterProps

        @property
        def saturation_temperature(self) -> Q_:
            return self.water_props.Tsat(self.pressure)

        @property
        def phase(self) -> str:
            Tsat = self.saturation_temperature
            if self.temperature < Tsat: return "liquid"
            if self.temperature > Tsat: return "vapor"
            return "saturated"

        @property
        def quality(self) -> Q_:
            P = self.pressure
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
            P = self.pressure
            T = self.temperature
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
            P = self.pressure
            T = self.temperature
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
            P = self.pressure
            T = self.temperature
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
            P = self.pressure
            T = self.temperature
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
            P = self.pressure
            T = self.temperature
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
            P = self.pressure
            # always use saturation T for surface tension when near boiling
            T_use = self.saturation_temperature if self.phase != "liquid" else self.temperature
            return self.water_props.sigma(P, T_use)

        @property
        def latent_heat_of_vaporization(self) -> Q_:
            P = self.pressure
            return self.water_props.h_fg(P)

        @property
        def prandt_number(self) -> Q_:
            return self.dynamic_viscosity * self.specific_heat / self.thermal_conductivity
        
        @property
        def velocity(self) -> Q_:
            return self.mass_flow_rate / self.stage.shell.flow_area
        
        @property
        def reynlods_number(self) -> Q_:
            return self.density * self.velocity * self.stage.hydraulic_diamtere / self.dynamic_viscosity

        @property
        def htc(self) -> Q_:
            return HTCFunctions.htc_shell(self)