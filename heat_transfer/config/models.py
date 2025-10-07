from __future__ import annotations
from dataclasses import dataclass
from common.units import ureg, Q_
from typing import Dict

@dataclass
class Wall:
    thickness: Q_
    conductivity: Q_

@dataclass
class Surface:
    roughness: Q_
    emissivity: Q_
    fouling_thickness: Q_
    fouling_conductivity: Q_

@dataclass
class Surfaces:
    inner: Surface
    outer: Surface

@dataclass
class DrumGeometry:
    inner_diameter: Q_
    inner_length: Q_
    wall: Wall

@dataclass
class PassGeometry:
    inner_diameter: Q_
    inner_length: Q_
    number_of_tubes: Q_
    layout: str
    pitch: Q_
    wall: Wall

@dataclass
class ReversalGeometry:
    inner_diameter: Q_
    inner_length: Q_
    wall: Wall

@dataclass
class Nozzle:
    diameter: Q_
    length: Q_


@dataclass
class ReversalNozzles:
    inlet: Nozzle
    outlet: Nozzle


@dataclass
class Drum:
    geometry: DrumGeometry
    surfaces: Surfaces


@dataclass
class Pass:
    geometry: PassGeometry
    surfaces: Surfaces


@dataclass
class Reversal:
    geometry: ReversalGeometry
    nozzles: ReversalNozzles
    surfaces: Surfaces


@dataclass
class Stages:
    drum: Drum
    pass1: Pass
    reversal1: Reversal
    pass2: Pass
    reversal2: Reversal
    pass3: Pass


@dataclass
class GasStream:
    mass_flow_rate: Q_
    temperature: Q_
    pressure: Q_
    composition: Dict[str, Q_]
    spectroscopic_data: Dict[str, Q_]
    z: Q_

@dataclass
class Water:
    mass_flow_rate: Q_
    temperature: Q_
    pressure: Q_
    composition: Dict[str, Q_]
    z: Q_


@dataclass
class Environment:
    ambient_temperature: Q_
    radiative_temperature: Q_
    external_emissivity: Q_
    external_h: Q_
    radiation_view_factor_external: Q_

@dataclass
class Config:
    gas_inlet: GasStream
    water_inlet: Water
    environment: Environment
    stages: Stages