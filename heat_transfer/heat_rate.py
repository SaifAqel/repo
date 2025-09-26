from math import log, pi
from common.units import ureg, Q_

class WaterResistance:
    def __init__(self, h_boil):
        self.h_boil = h_boil

    def resistance_per_area(self):
        return (1 / self.h_boil).to("kelvin*meter**2/watt")


class GasResistance:
    def __init__(self, h_rad, h_conv):
        self.h_rad = h_rad
        self.h_conv = h_conv

    def resistance_per_area(self):
        return (1 / (self.h_rad + self.h_conv)).to("kelvin*meter**2/watt")


class WallResistance:
    def __init__(self, r_i, r_o, k):
        self.r_i = r_i
        self.r_o = r_o
        self.k = k

    def cyl_per_inner_area(self):
        value = log((self.r_o / self.r_i).to_base_units().magnitude)
        return Q_(value, "") / (2 * pi * self.k * self.r_i).to("watt/kelvin")


class FoulingResistance:
    def __init__(self, thickness, k_foul):
        self.thickness = thickness
        self.k_foul = k_foul

    def resistance_per_area(self):
        return (self.thickness / self.k_foul).to("kelvin*meter**2/watt")


class TotalResistance:
    def __init__(self, r_water, r_wall, r_gas, r_foul_inner, r_foul_outer):
        self.r_water = r_water
        self.r_wall = r_wall
        self.r_gas = r_gas
        self.r_foul_inner = r_foul_inner
        self.r_foul_outer = r_foul_outer

    def resistance_per_area(self):
        return (self.r_water + self.r_wall + self.r_gas + self.r_foul_inner + self.r_foul_outer).to("kelvin*meter**2/watt")


class Heat_Flux:
    def __init__(self, delta_T, r_total_per_area):
        self.delta_T = delta_T
        self.r_total_per_area = r_total_per_area

    def resistance_per_area(self):
        return (self.delta_T / self.r_total_per_area).to("watt/meter**2")


class Heat_Rate:
    def compute(self, flux, area):
        return (flux * area).to("watt")
