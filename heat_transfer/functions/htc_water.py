# water_htc.py
from typing import Literal
from common.units import Q_
from heat_transfer.config.models import WaterStream, GasStream, FirePass, SmokePass, Reversal, Economiser

Zone  = Literal["firepass", "smokepass", "reversal", "economiser"]

class WaterHTC:

    ######################### Nusselt Number #########################
    @staticmethod
    def zone(water) -> Zone:
        return water.stage.__class__.__name__.lower()
    
    @staticmethod
    def Nu_zukauskas(water):
        C = 0.27
        m = 0.63
        F_row = 1.0
        F_arr = 1.0
        Nu = C * (water.film.reynolds_gap.to("dimensionless").magnitude**m) * (water.film.prandtl_number.to("dimensionless").magnitude**0.36) * F_row * F_arr
        return Nu
    
    @staticmethod
    def Nu_churchill_bernstein(water):
        Nu = 0.3 + (
            (0.62 * (water.film.reynolds_number.to("dimensionless").magnitude**0.5) * (water.film.prandtl_number.to("dimensionless").magnitude**(1/3)))
            / ((1 + (0.4 / water.film.prandtl_number.to("dimensionless").magnitude)**(2/3))**0.25)
        ) * ((1 + (water.film.reynolds_number.to("dimensionless").magnitude / 282000)**(5/8))**(4/5))
        return Nu
    
    @staticmethod
    def Nu_sieder_tate(water):
        Nu = 0.027 * (water.film.reynolds_number.to("dimensionless").magnitude**0.8) * (water.film.prandtl_number.to("dimensionless").magnitude**(1/3))
        return Nu
    
    @staticmethod
    def calc_Nu(water):
        zone = WaterHTC.zone(water)
        if zone == "smokepass":
            return WaterHTC.Nu_zukauskas(water)
        elif zone == "firepass":
            return WaterHTC.Nu_churchill_bernstein(water)
        elif zone == "reversal":
            Nu_cb = WaterHTC.Nu_churchill_bernstein(water)
            De = water.film.reynolds_number * (water.stage.hot_side.inner_diameter / (2 * water.stage.hot_side.curvature_radius))
            return Nu_cb * (1 + 0.15 * (De ** 0.5))
        elif zone == "economiser":
            return WaterHTC.Nu_sieder_tate(water)
        else:
            raise ValueError(f"Unknown zone type: {zone}")
        
    ######################### HTC #########################
    @staticmethod
    def htc_conv(water):
        return WaterHTC.calc_Nu(water) * water.film.thermal_conductivity / water.stage.cold_side.hydraulic_diameter
    
    @staticmethod
    def htc_nb(water):
        p_crit = Q_(22.064e6, "Pa") 
        h_nb = 55 * ((water.pressure / p_crit).magnitude ** -0.12) * (water.molecular_weight.magnitude ** -0.5) * (water.q_flux.magnitude ** 0.67)
        return h_nb


    def calc_htc(water):
        if water.enthalpy < water.liquid_saturation_enthalpy:
            return WaterHTC.htc_conv(water)
        else:
            return ( water.S_factor * WaterHTC.htc_conv(water) ) + ( water.F_factor * WaterHTC.htc_nb(water) )
