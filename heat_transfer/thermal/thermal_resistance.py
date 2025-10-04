from math import log, pi
from common.units import ureg, Q_

class ThermalResistance:
    def __init__(self, PassWithCalc, GasStream, PoolBoilingHTC, h_rad):
        self.h_boil = PoolBoilingHTC.rohsenow_h
        self.h_conv = GasStream.Convective_coefficient
        self.h_rad = h_rad
        self.r_i = PassWithCalc.r_i
        self.r_o = PassWithCalc.r_o
        self.k = PassWithCalc.wall.k_wall
        self.thickness_inner = PassWithCalc.thickness_inner
        self.k_foul_inner = PassWithCalc.k_foul_inner
        self.thickness_outer = PassWithCalc.thickness_outer
        self.k_foul_outer = PassWithCalc.k_foul_outer

    def water_resistance(self):
        return (1 / self.h_boil).to("kelvin*meter**2/watt")

    def gas_resistance(self):
        return (1 / (self.h_rad + self.h_conv)).to("kelvin*meter**2/watt")

    def wall_resistance(self):
        value = log((self.r_o / self.r_i).to_base_units().magnitude)
        return Q_(value, "") / (2 * pi * self.k * self.r_i).to("watt/kelvin")

    def fouling_resistance_inner(self):
        return (self.thickness_inner / self.k_foul_inner).to("kelvin*meter**2/watt")

    def fouling_resistance_outer(self):
        return (self.thickness_outer / self.k_foul_outer).to("kelvin*meter**2/watt")

    def total_resistance(self):
        return (
            self.water_resistance()
            + self.wall_resistance()
            + self.gas_resistance()
            + self.fouling_resistance_inner()
            + self.fouling_resistance_outer()
        ).to("kelvin*meter**2/watt")


