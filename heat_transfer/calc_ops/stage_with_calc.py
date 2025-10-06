from heat_transfer.config.schemas import (Config, Stages, Drum, Pass, Reversal)
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
        return self.geometry.inner_length

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
class StagesWithCalc(Stages):
    drum: DrumWithCalc
    pass1: PassWithCalc
    reversal1: Reversal
    pass2: PassWithCalc
    reversal2: Reversal
    pass3: PassWithCalc


@dataclass
class ConfigWithCalc(Config):
    stages: StagesWithCalc


def with_calc(cfg: Config) -> ConfigWithCalc:
    s = cfg.stages
    stages_calc = StagesWithCalc(
        drum=DrumWithCalc(geometry=s.drum.geometry, surfaces=s.drum.surfaces),
        pass1=PassWithCalc(geometry=s.pass1.geometry, surfaces=s.pass1.surfaces),
        reversal1=s.reversal1,
        pass2=PassWithCalc(geometry=s.pass2.geometry, surfaces=s.pass2.surfaces),
        reversal2=s.reversal2,
        pass3=PassWithCalc(geometry=s.pass3.geometry, surfaces=s.pass3.surfaces),
    )
    return ConfigWithCalc(
        gas_side=cfg.gas_side,
        water_side=cfg.water_side,
        environment=cfg.environment,
        stages=stages_calc,
    )