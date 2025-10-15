# water_htc.py
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Literal, List

from common.units import Q_
from heat_transfer.config.models import WaterStream, GasStream, FirePass, SmokePass, Reversal, Economiser

Zone  = Literal["firepass", "smokepass", "reversal", "economiser"]

class WaterHTC:
    def __init__(self, stage: FirePass | SmokePass | Reversal | Economiser, water: WaterStream, gas: GasStream):
        self.stage = stage
        self.water = water
        self.gas = gas

    def zone(self) -> Zone:
        return self.stage.__class__.__name__.lower()

    def Nu_zukauskas(self):
        """Zukauskas correlation for crossflow over tube banks (turbulent)."""
        C = 0.27
        m = 0.63
        F_row = 1.0
        F_arr = 1.0
        Nu = C * (self.water.reynolds_gap**m) * (self.water.prandtl_number**0.36) * F_row * F_arr
        return Nu

    def Nu_churchill_bernstein(self):
        Nu = 0.3 + (
            (0.62 * (self.water.reynolds_number**0.5) * (self.water.prandtl_number**(1/3)))
            / ((1 + (0.4 / self.water.prandtl_number)**(2/3))**0.25)
        ) * ((1 + (self.water.reynolds_number / 282000)**(5/8))**(4/5))
        return Nu

    def Nu_sieder_tate(self):
        Nu = 0.027 * (self.water.reynolds_number**0.8) * (self.water.prandtl_number**(1/3))
        return Nu

    def calc_Nu(self):
        zone = self.zone
        if zone == "smokepass":
            return self.Nu_zukauskas()

        elif zone == "firepass":
            return self.Nu_churchill_bernstein()

        elif zone == "reversal":
            Nu_cb = self.Nu_churchill_bernstein()
            De = self.water.reynolds_number * (self.stage.hot_side.inner_diameter / (2 * self.stage.hot_side.curvature_radius))
            return Nu_cb * (1 + 0.15 * (De ** 0.5))

        elif zone == "economiser":
            return self.Nu_sieder_tate()

        else:
            raise ValueError(f"Unknown zone type: {zone}")

    def htc_conv(self):
        return self.calc_Nu * self.water.thermal_conductivity / self.stage.cold_side.hydraulic_diameter
    
    def htc_nb(self, q_flux):
        """
        Cooper nucleate-boiling correlation for water.
        Returns nucleate boiling heat transfer coefficient [W/m2-K].
        """
        p_crit = Q_(22.064e6, "Pa") 

        h_nb = 55 * ((self.water.pressure / p_crit) ** -0.12) * (self.water.molecular_weight ** -0.5) * (q_flux ** 0.67)
        return h_nb


    def calc_htc(self):
        if self.water.enthalpy < self.water.liquid_saturation_enthalpy:
            return self.htc_conv

        else:
            return self.htc_conv + self.htc_nb
