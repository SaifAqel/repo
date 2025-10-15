from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, Optional, Literal
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
    def HX_area(self) -> Q_:
        p = pi * self.inner_diameter
        return p * self.inner_length
    
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
    def HX_area(self) -> Q_:
        p = pi * self.inner_diameter
        return p * self.inner_length
    
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
    def outer_diameter(self): return self.inner_diameter + 2*self.wall.thickness

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


######################### Streams #########################

@dataclass
class GasStream:
        mass_flow_rate: Q_
        temperature: Q_
        pressure: Q_
        composition: Dict[str, Q_]
        spectroscopic_data: Dict[str, Q_]

        stage: FirePass | SmokePass | Reversal
        gas_props: GasProps

        @property
        def density(self) -> Q_:
            return self.gas_props.density(self)
        
        @property
        def specific_heat(self) -> Q_:
            return self.gas_props.specific_heat(self)
        
        @property
        def dynamic_viscosity(self) -> Q_:
            return self.gas_props.viscosity(self)
        
        @property
        def thermal_conductivity(self) -> Q_:
            return self.gas_props.thermal_conductivity(self)
        
        @property
        def mass_flux(self) -> Q_:
            return self.mass_flow_rate / self.stage.hot_side.flow_area

        @property
        def velocity(self) -> Q_:
            return self.mass_flow_rate / (self.density * self.stage.hot_side.flow_area)

        @property
        def reynolds_number(self) -> Q_:
            return (self.density * self.velocity * self.stage.hot_side.hydraulic_diameter) / self.dynamic_viscosity
        
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
            return self.nusselt_number * self.thermal_conductivity / self.stage.hot_side.inner_diameter
        
        @property
        def absorption_coefficient (self) -> Q_:
            return sum(self.composition[s] * self.spectroscopic_data[s] for s in self.composition)

        @property
        def emissivity(self) -> Q_:
            return Q_(1.0 - exp((-self.absorption_coefficient * self.stage.hot_side.path_length).magnitude), "dimensionless")
            
        def radiation_coefficient(self, T_wall) -> Q_:
            sigma = 5.670374419e-8 * ureg.watt / (ureg.meter**2 * ureg.kelvin**4)
            mean_temperature = (self.temperature + T_wall) / 2
            return 4.0 * sigma * (mean_temperature**3) * self.emissivity

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

    stage: FirePass | SmokePass | Reversal
    water_props: WaterProps

    @property
    def quality(self) -> Q_:
        return self.water_props.quality(self)
    
    @property
    def molecular_weight(self) -> Q_:
        return Q_(18.0, "kg/kmol")
 
    @property
    def saturation_temperature(self) -> Q_:
        return self.water_props.saturation_temperature(self)
    
    @property
    def liquid_saturation_enthalpy(self) -> Q_:
        return self.water_props.saturation_enthalpy_liquid(self)

    @property
    def temperature(self) -> Q_:
        return self.water_props.temperature(self)

    @property
    def density(self) -> Q_:
        return self.water_props.density(self)

    @property
    def specific_heat(self) -> Q_:
        return self.water_props.specific_heat(self)

    @property
    def thermal_conductivity(self) -> Q_:
        return self.water_props.thermal_conductivity(self)

    @property
    def dynamic_viscosity(self) -> Q_:
        return self.water_props.dynamic_viscosity(self)

    @property
    def surface_tension(self) -> Q_:
        return self.water_props.surface_tension(self)

    @property
    def latent_heat_of_vaporization(self) -> Q_:
        return self.water_props.latent_heat_of_vaporization(self)

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