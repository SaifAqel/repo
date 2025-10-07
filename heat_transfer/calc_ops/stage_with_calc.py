from heat_transfer.config.models import (Config, Stages, Drum, Pass, Reversal)
from dataclasses import dataclass
from common.units import ureg, Q_
from math import pi

@dataclass
class PassWithCalc(Pass):

    @property
    def tube_inner_flow_area(self) -> Q_:
        di = self.geometry.inner_diameter
        return self.geometry.number_of_tubes * pi * (di/2)**2

    @property
    def tube_inner_perimeter(self) -> Q_:
        di = self.geometry.inner_diameter
        return self.geometry.number_of_tubes * pi * di
    
    @property
    def tube_outer_perimeter(self) -> Q_:
        do = self.outer_diameter
        return pi * do

    @property
    def tube_inner_heat_transfer_area(self) -> Q_:
        Pi = self.tube_inner_perimeter
        L = self.geometry.inner_length
        return Pi * L

    @property
    def outer_diameter(self) -> Q_:
        di = self.geometry.inner_diameter
        k = self.geometry.wall.thickness
        return di + (2 * k)

    @property
    def cross_section_outer_area(self) -> Q_:
        do = self.outer_diameter
        return pi * (do/2)**2
    
    @property
    def rel_roughness(self) -> Q_:
        return self.surfaces.inner.roughness / self.geometry.inner_diameter
    
    @property
    def path_length(self) -> Q_:
        """Characteristic path length for radiation calculations."""
        return self.geometry.inner_diameter * 0.9

@dataclass
class DrumWithCalc(Drum):

    @property
    def cross_section_inner_area(self) -> Q_:
        di = self.geometry.inner_diameter
        return pi * (di/2)**2
    
    @property
    def outer_diameter(self) -> Q_:
        di = self.geometry.inner_diameter
        k = self.geometry.wall.thickness
        return di + (2 * k)

    @property
    def cross_section_outer_area(self) -> Q_:
        do = self.outer_diameter
        return pi * (do/2)**2
    
@dataclass
class ReversalWithCalc(Reversal):
    @property
    def tube_inner_flow_area(self) -> Q_:
        di = self.geometry.inner_diameter
        return pi * (di/2)**2  # Single large chamber

    @property
    def tube_inner_perimeter(self) -> Q_:
        di = self.geometry.inner_diameter
        return pi * di
    
    @property
    def tube_outer_perimeter(self) -> Q_:
        do = self.outer_diameter
        return pi * do

    @property
    def tube_inner_heat_transfer_area(self) -> Q_:
        Pi = self.tube_inner_perimeter
        L = self.geometry.inner_length
        return Pi * L

    @property
    def outer_diameter(self) -> Q_:
        di = self.geometry.inner_diameter
        k = self.geometry.wall.thickness
        return di + (2 * k)

    @property
    def cross_section_outer_area(self) -> Q_:
        do = self.outer_diameter
        return pi * (do/2)**2
    
    @property
    def rel_roughness(self) -> Q_:
        return self.surfaces.inner.roughness / self.geometry.inner_diameter
    
    @property
    def path_length(self) -> Q_:
        return self.geometry.inner_diameter  # Mean beam ~di for chamber

@dataclass
class StagesWithCalc(Stages):
    drum: DrumWithCalc
    pass1: PassWithCalc
    reversal1: ReversalWithCalc
    pass2: PassWithCalc
    reversal2: ReversalWithCalc
    pass3: PassWithCalc


@dataclass
class ConfigWithCalc(Config):
    stages: StagesWithCalc


def with_calc(cfg: Config) -> ConfigWithCalc:
    s = cfg.stages
    stages_calc = StagesWithCalc(
        drum=DrumWithCalc(geometry=s.drum.geometry, surfaces=s.drum.surfaces),
        pass1=PassWithCalc(geometry=s.pass1.geometry, surfaces=s.pass1.surfaces),
        reversal1=ReversalWithCalc(geometry=s.reversal1.geometry, surfaces=s.reversal1.surfaces, nozzles=s.reversal1.nozzles),
        pass2=PassWithCalc(geometry=s.pass2.geometry, surfaces=s.pass2.surfaces),
        reversal2=ReversalWithCalc(geometry=s.reversal2.geometry, surfaces=s.reversal2.surfaces, nozzles=s.reversal2.nozzles),
        pass3=PassWithCalc(geometry=s.pass3.geometry, surfaces=s.pass3.surfaces),
    )
    return ConfigWithCalc(
        gas_inlet=cfg.gas_inlet,
        water_inlet=cfg.water_inlet,
        environment=cfg.environment,
        stages=stages_calc,
    )