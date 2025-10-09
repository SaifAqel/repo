from dataclasses import dataclass
from pint import Quantity as    Q_
from math import pi

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
class TubeBank:
    geom: TubeGeometry
    layout: BankLayout

@dataclass(frozen=True)
class ReversalChamber:
    geometry: ReversalChamberGeometry
    nozzles: Nozzles







@dataclass(frozen=True)
class FirePass:
    geom: TubeGeometry
    wall: Wall
    shell: Shell

@dataclass(frozen=True)
class SmokePass:
    geom: TubeBank
    wall: Wall
    shell: Shell

@dataclass(frozen=True)
class Reversal:
    geom: ReversalChamber
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

