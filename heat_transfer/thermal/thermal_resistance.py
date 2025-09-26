from math import log, pi
from common.units import ureg, Q_

class ThermalResistance:
    def __init__(self, h_boil, h_rad, h_conv, r_i, r_o, k, thickness_inner, k_foul_inner, thickness_outer, k_foul_outer):
        self.h_boil = h_boil
        self.h_rad = h_rad
        self.h_conv = h_conv
        self.r_i = r_i
        self.r_o = r_o
        self.k = k
        self.thickness_inner = thickness_inner
        self.k_foul_inner = k_foul_inner
        self.thickness_outer = thickness_outer
        self.k_foul_outer = k_foul_outer

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


