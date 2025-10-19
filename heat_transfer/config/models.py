from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, Callable
from common.units import ureg, Q_
from math import exp, log10, pi
from heat_transfer.functions.fluid_props import WaterProps, GasProps

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
    k: Q_
    
@dataclass(frozen=True)
class Nozzles:
    inlet: Nozzle
    outlet: Nozzle

@dataclass(frozen=True)
class TubeGeometry:
    inner_diameter: Q_
    inner_length: Q_
    wall: Wall

    @property
    def outer_diameter(self) -> Q_:
        return self.inner_diameter + ( 2 * self.wall.thickness )

    @property
    def flow_area(self) -> Q_:
        return pi * (self.inner_diameter / 2 )**2
    
    @property
    def rel_roughness(self) -> Q_:
        return self.wall.surfaces.inner.roughness / self.inner_diameter
    
    @property
    def inner_perimeter(self) -> Q_:
        return pi * self.inner_diameter
    
    @property
    def outer_perimeter(self) -> Q_:
        return pi * self.outer_diameter
    
    @property
    def path_length(self) -> Q_:
        return self.inner_diameter * 0.9
    
    @property
    def hydraulic_diameter(self) -> Q_:
        return self.inner_diameter

@dataclass(frozen=True)
class ReversalGeometry:
    inner_length: Q_
    inner_diameter: Q_
    curvature_radius: Q_
    nozzles: Nozzles
    wall: Wall

    @property
    def outer_diameter(self) -> Q_:
        return self.inner_diameter + ( 2 * self.wall.thickness)
    
    @property
    def flow_area(self) -> Q_:
        return pi * (self.inner_diameter / 2 )**2

    @property
    def rel_roughness(self) -> Q_:
        return self.wall.surfaces.inner.roughness / self.inner_diameter
    
    @property
    def inner_perimeter(self) -> Q_:
        return pi * self.inner_diameter
    
    @property
    def outer_perimeter(self) -> Q_:
        return pi * self.outer_diameter
    
    @property
    def path_length(self) -> Q_:
        return self.inner_diameter * 0.9
    
    @property
    def hydraulic_diameter(self) -> Q_:
        return self.inner_diameter
    
@dataclass(frozen=True)
class BankGeometry:
    inner_diameter: Q_
    inner_length: Q_
    tubes_number: Q_
    layout: str
    pitch: Q_
    wall: Wall

    @property
    def flow_area(self) -> Q_:
        return self.tubes_number * pi * (self.inner_diameter / 2 )**2
    
    @property
    def rel_roughness(self) -> Q_:
        return self.wall.surfaces.inner.roughness / self.inner_diameter
    
    @property
    def outer_diameter(self) -> Q_: 
        return self.inner_diameter + 2*self.wall.thickness
    
    @property
    def inner_perimeter(self) -> Q_:
        return pi * self.inner_diameter
    
    @property
    def outer_perimeter(self) -> Q_:
        return pi * self.outer_diameter

    @property
    def hydraulic_diameter(self) -> Q_:
        return self.inner_diameter
    
    @property
    def path_length(self) -> Q_:
        return self.inner_diameter * 0.9

@dataclass(frozen=True)
class EconomiserHot:
    inner_length: Q_
    inner_diameter: Q_
    wall: Wall

    @property
    def flow_area(self) -> Q_:
        return  pi * (self.inner_diameter / 2 )**2
    
    @property
    def rel_roughness(self) -> Q_:
        return self.wall.surfaces.inner.roughness / self.inner_diameter
    
    @property
    def outer_diameter(self): return self.inner_diameter + 2*self.wall.thickness

    @property
    def hydraulic_diameter(self) -> Q_:
        return self.inner_diameter
    
    @property
    def path_length(self) -> Q_:
        return self.inner_diameter * 0.9


@dataclass(frozen=True)
class EconomiserCold:
    inner_length: Q_
    inner_diameter: Q_
    wall: Wall

    @property
    def flow_area(self) -> Q_:
        return pi * (self.inner_diameter / 2 )**2
    
    @property
    def rel_roughness(self) -> Q_:
        return self.wall.surfaces.inner.roughness / self.inner_diameter
    
    @property
    def outer_diameter(self): return self.inner_diameter + 2*self.wall.thickness

    @property
    def hydraulic_diameter(self) -> Q_:
        return self.inner_diameter
    
    @property
    def path_length(self) -> Q_:
        return self.inner_diameter * 0.9

@dataclass(frozen=True)
class ShellGeometry:
    flow_area: Q_
    wetted_perimeter: Q_
    wall: Wall

    @property
    def hydraulic_diameter(self) -> Q_:
        return 4 * self.flow_area / self.wetted_perimeter
    
@dataclass(frozen=True)
class Drum:
    flow_area: Q_

######################### Stages #########################

@dataclass(frozen=True)
class FirePass:
    hot_side: TubeGeometry
    cold_side: ShellGeometry

@dataclass(frozen=True)
class SmokePass:
    hot_side: BankGeometry
    cold_side: ShellGeometry

@dataclass(frozen=True)
class Reversal:
    hot_side: ReversalGeometry
    cold_side: ShellGeometry

@dataclass(frozen=True)
class Economiser:
    hot_side: EconomiserHot
    cold_side: EconomiserCold

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
    stage: FirePass | SmokePass | Reversal | Economiser
    wall_temperature: Q_ | None = None  

    ######################### Properties #########################
    @property
    def density(self) -> Q_:
        return GasProps.density(self)
    
    @property
    def specific_heat(self) -> Q_:
        return GasProps.specific_heat(self)
    
    @property
    def dynamic_viscosity(self) -> Q_:
        return GasProps.viscosity(self)
    
    @property
    def thermal_conductivity(self) -> Q_:
        return GasProps.thermal_conductivity(self)
    
    ######################### Flow #########################
    @property
    def mass_flux(self) -> Q_:
        return self.mass_flow_rate / self.stage.hot_side.flow_area

    @property
    def velocity(self) -> Q_:
        return self.mass_flow_rate / (self.density * self.stage.hot_side.flow_area)

    @property
    def reynolds_number(self) -> Q_:
        return (self.density * self.velocity * self.stage.hot_side.hydraulic_diameter) / self.dynamic_viscosity
    
    ######################### Convection #########################
    @property
    def prandtl_number(self) -> Q_:
        return self.dynamic_viscosity * self.specific_heat / self.thermal_conductivity

    @property
    def nusselt_number(self) -> Q_:
        n = 0.3
        Re_mag = self.reynolds_number.magnitude
        Pr_mag = self.prandtl_number.magnitude
        nu_mag = 0.023 * (Re_mag ** 0.8) * (Pr_mag ** n)
        return Q_(nu_mag, ureg.dimensionless)
    
    @property
    def convective_coefficient(self) -> Q_:
        return self.nusselt_number * self.thermal_conductivity / self.stage.hot_side.hydraulic_diameter
    
    ######################### Radiation #########################
    @property
    def absorption_coefficient (self) -> Q_:
        return sum(self.composition[s] * self.spectroscopic_data[s] for s in self.composition)

    @property
    def emissivity(self) -> Q_:
        return Q_(1.0 - exp((-self.absorption_coefficient * self.stage.hot_side.path_length).magnitude), "dimensionless")
    
    @property
    def film_temperature(self) -> Q_:
        return (self.wall_temperature + self.temperature) / 2
    
    @property
    def radiation_coefficient(self) -> Q_:
        sigma = 5.670374419e-8 * ureg.watt / (ureg.meter**2 * ureg.kelvin**4)
        return 4.0 * sigma * (self.film_temperature**3) * self.emissivity
    
    ######################### HTC #########################
    @property
    def htc(self) -> Q_:
        return self.radiation_coefficient + self.convective_coefficient
    
    ######################### Friction #########################
    @property
    def friction_factor(self) -> float:
        f = ( -1.8*log10( (self.stage.hot_side.rel_roughness.magnitude/3.7)**1.11 + 6.9/self.reynolds_number.magnitude ) )**-2
        return f


@dataclass
class WaterStream:
    mass_flow_rate: Q_
    enthalpy: Q_
    pressure: Q_
    composition: Dict[str, Q_]
    drum: Drum
    stage: FirePass | SmokePass | Reversal | Economiser
    wall_temperature: Q_ | None = None  # Two
    q_flux: Q_ | None = None

    # --- properties ---
    @property
    def quality(self) -> Q_:
        return WaterProps.quality_from_h(self)
    
    @property
    def temperature(self) -> Q_:
        return WaterProps.temperature(self)

    @property
    def density(self) -> Q_:
        return WaterProps.density(self)

    @property
    def specific_heat(self) -> Q_:
        return WaterProps.specific_heat_cp(self)

    @property
    def thermal_conductivity(self) -> Q_:
        return WaterProps.thermal_conductivity(self)

    @property
    def dynamic_viscosity(self) -> Q_:
        return WaterProps.dynamic_viscosity(self)
    
    @property
    def molecular_weight(self) -> Q_:
        return Q_(18.0, "kg/kmol")
    
    # --- Film ---
    @property
    def film(self) -> Film:
        return Film(bulk=self)

    # --- Saturation Properties ---
    @property
    def saturation_temperature(self) -> Q_:
        return WaterProps.saturation_temperature(self)
    
    @property
    def liquid_saturation_enthalpy(self) -> Q_:
        return WaterProps.saturation_enthalpy_liquid(self)
    
    @property
    def latent_heat_of_vaporization(self) -> Q_:
        return WaterProps.latent_heat(self)

    @property
    def surface_tension(self) -> Q_:
        return WaterProps.surface_tension(self)

    # --- Flow ---
    @property
    def prandtl_number(self) -> Q_:
        return self.dynamic_viscosity * self.specific_heat / self.thermal_conductivity

    @property
    def velocity(self) -> Q_:
        return self.mass_flow_rate / (self.stage.cold_side.flow_area * self.density)

    @property
    def reynolds_number(self) -> Q_:
        return self.density * self.velocity * self.stage.cold_side.hydraulic_diameter / self.dynamic_viscosity

    @property
    def reynolds_gap(self) -> Q_:
        gap_velocity = self.velocity * self.stage.hot_side.pitch / (self.stage.hot_side.pitch - self.stage.hot_side.outer_diameter)
        return self.density * gap_velocity * self.stage.hot_side.outer_diameter / self.dynamic_viscosity
    
    @property
    def mass_flux(self) -> Q_:
        return self.mass_flow_rate / self.drum.flow_area

    # --- Two-Phase ---
    @property
    def Re_lo(self):
        return (self.mass_flux * (1 - self.quality) * self.stage.cold_side.hydraulic_diameter / WaterProps.sat_liq(self).mu).to_base_units().magnitude
    
    @property
    def xtt(self):
        return ((1 - self.quality) / self.quality) ** 0.9 * (WaterProps.sat_vap(self).rho / WaterProps.sat_liq(self).rho) ** 0.5 * (WaterProps.sat_liq(self).mu / WaterProps.sat_vap(self).mu) ** 0.1
    
    @property
    def F_factor(self):
        return (1 / self.xtt + 0.213) ** 0.736
    
    @property
    def S_factor(self):
        return 1 / (1 + 2.53e-6 * self.Re_lo ** 1.17)
    
@dataclass
class Film:
    bulk: "WaterStream"

    # --- Properties ---
    @property
    def temperature(self) -> Q_:
        return (self.bulk.wall_temperature + self.bulk.temperature) / 2
    
    @property
    def pressure(self) -> Q_:
        return self.bulk.pressure

    @property
    def composition(self):
        return self.bulk.composition

    @property
    def enthalpy(self) -> Q_:
        return WaterProps.enthalpy_from_temperature(self)

    @property
    def density(self) -> Q_:
        return WaterProps.density(self)

    @property
    def specific_heat(self) -> Q_:
        return WaterProps.specific_heat_cp(self)

    @property
    def thermal_conductivity(self) -> Q_:
        return WaterProps.thermal_conductivity(self)

    @property
    def dynamic_viscosity(self) -> Q_:
        return WaterProps.dynamic_viscosity(self)

    # --- Flow ---
    @property
    def prandtl_number(self) -> Q_:
        return self.dynamic_viscosity * self.specific_heat / self.thermal_conductivity

    @property
    def velocity(self) -> Q_:
        return self.bulk.velocity

    @property
    def hydraulic_diameter(self) -> Q_:
        return self.bulk.stage.cold_side.hydraulic_diameter

    @property
    def reynolds_number(self) -> Q_:
        return self.density * self.velocity * self.hydraulic_diameter / self.dynamic_viscosity

    @property
    def reynolds_gap(self) -> Q_:
        hs = self.bulk.stage.hot_side
        gap_velocity = self.velocity * hs.pitch / (hs.pitch - hs.outer_diameter)
        return self.density * gap_velocity * hs.outer_diameter / self.dynamic_viscosity
