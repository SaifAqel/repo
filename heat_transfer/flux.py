class WaterResistance:
    @staticmethod
    def per_area(h_boil):
        return 1 / h_boil

class GasResistance:
    @staticmethod
    def per_area(h_rad, h_conv):
        return 1 / (h_rad + h_conv)

class WallResistance:
    @staticmethod
    def cyl_logmean_per_area(r_i, r_o, k):
        return (r_o - r_i) / k

class FoulingResistance:
    @staticmethod
    def per_area(r_foul):
        return r_foul

class TotalResistance:
    @staticmethod
    def per_area(r_water, r_wall, r_gas, r_foul_inner, r_foul_outer):
        return r_water + r_wall + r_gas + r_foul_inner + r_foul_outer

class Flux:
    @staticmethod
    def per_area(delta_T, r_total_per_area):
        return delta_T / r_total_per_area