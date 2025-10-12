# water_htc.py
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Literal, List

from common.units import Q_
from heat_transfer.config.models import WaterStream, GasStream, FirePass, SmokePass, Reversal, Economiser

# -------- types --------
Zone  = Literal["firepass", "smokepass", "reversal", "economiser"]
Phase = Literal["subcooled", "twophase", "superheated"]
Regime = Literal["laminar", "turbulent"]


@dataclass(frozen=True)
class BankBand:
    """Correlation band for Zukauskas-style tube-bank crossflow."""
    re_min: float
    re_max: float
    c: float
    m: float
    n: float


class WaterHTC:
    """
    Compute water-side HTC with:
      - Phase switching: subcooled/superheated -> single-phase; twophase -> Gungor-Winterton.
      - Zone switching: firepass/smokepass/economiser (tube-bank crossflow) or reversal (curved crossflow).
      - Regime switching: laminar vs turbulent attributes where provided, with fallback to generic.
    Attributes required on `self` are accessed dynamically as <base>_<zone>[_<regime>].
    """

    def __init__(self, stage: FirePass | SmokePass | Reversal | Economiser, water: WaterStream, gas: GasStream):
        self.stage = stage
        self.water = water
        self.gas = gas

    # ---------- phase / regime / zone ----------
    def phase_case(self) -> Phase:
        x = self.water.quality.magnitude
        if x < 0:
            return "subcooled"
        if x <= 1:
            return "twophase"
        return "superheated"

    def flow_regime(self) -> Regime:
        return "turbulent" if self.water.reynolds_number >= Q_(1200, "dimensionless") else "laminar"

    def zone(self) -> Zone:
        return self.stage.__class__.__name__.lower()

    # ---------- public API ----------
    def compute_htc(self):
        z = self.zone()
        p = self.phase_case()
        if p in {"subcooled", "superheated"}:
            return self._h_singlephase(z)
        if p == "twophase":
            return self._h_twophase(z)
        raise ValueError(f"Unknown phase: {p}")

    # ---------- regime-aware attribute picker ----------
    def _pick(self, *, zone: Zone, base: str):
        """
        Return regime-specific attribute if present, else generic.
        Looks for <base>_<zone>_<regime> then <base>_<zone>.
        """
        regime = self.flow_regime()
        name_regime = f"{base}_{zone}_{regime}"
        if hasattr(self, name_regime):
            return getattr(self, name_regime)
        name_generic = f"{base}_{zone}"
        if hasattr(self, name_generic):
            return getattr(self, name_generic)
        raise AttributeError(f"Missing attribute for base='{base}' zone='{zone}' "
                             f"(looked for '{name_regime}' and '{name_generic}')")

    # ---------- correlations (static) ----------
    def nusselt_churchill_bernstein(self):
        term1 = 0.3
        term2 = (0.62 * self.water.reynolds_number**0.5 * self.water.prandtl_number**(1/3)) / ((1 + (0.4 / self.water.prandtl_number)**(2/3))**0.25)
        term3 = (1 + (self.water.reynolds_number / 282000)**(5/8))**(4/5)
        return term1 + term2 * term3

    @staticmethod
    def heat_transfer_coefficient_from_nusselt(nusselt, thermal_conductivity, characteristic_diameter):
        return nusselt * thermal_conductivity / characteristic_diameter

    @staticmethod
    def nusselt_zukauskas(reynolds_gap, prandtl, prandtl_surface,
                           coefficient_c, exponent_m, exponent_n,
                           row_correction_factor, arrangement_correction_factor):
        return (
            coefficient_c
            * reynolds_gap**exponent_m
            * prandtl**exponent_n
            * (prandtl / prandtl_surface)**0.25
            * row_correction_factor
            * arrangement_correction_factor
        )

    @staticmethod
    def curvature_correction_factor_dean(reynolds, diameter, curvature_radius, coefficient_a, exponent_b):
        dean_number = reynolds * math.sqrt(diameter / (2.0 * curvature_radius))
        return 1.0 + coefficient_a * dean_number**exponent_b

    @staticmethod
    def nusselt_churchill_bernstein_with_curvature(reynolds, prandtl, diameter, curvature_radius, coefficient_a, exponent_b):
        nu_base = WaterHTC.nusselt_churchill_bernstein(reynolds, prandtl)
        factor = WaterHTC.curvature_correction_factor_dean(reynolds, diameter, curvature_radius, coefficient_a, exponent_b)
        return nu_base * factor

    @staticmethod
    def gungor_winterton_boiling(h_singlephase, h_nucleate, mass_flux, heat_flux, latent_heat, vapor_quality,
                                 liquid_density, vapor_density, liquid_viscosity, vapor_viscosity, surface_tension, diameter):
        # Units must be consistent (Quantities OK).
        reynolds_liquid_only = mass_flux * (1.0 - vapor_quality) * diameter / liquid_viscosity
        boiling_number = heat_flux / (mass_flux * latent_heat)
        x_tt = ((1.0 - vapor_quality) / vapor_quality)**0.9 \
               * (vapor_density / liquid_density)**0.5 \
               * (liquid_viscosity / vapor_viscosity)**0.1
        enhancement_factor = 1.0 + 3000.0 * boiling_number**0.86 + 1.12 * x_tt**(-1.86)
        suppression_factor = 1.0 / (1.0 + 1.15e-6 * reynolds_liquid_only**1.17)
        return suppression_factor * h_nucleate + enhancement_factor * h_singlephase

    @staticmethod
    def _nusselt_bank_banded(reynolds_gap, prandtl, prandtl_surface,
                             bands: List[BankBand], row_factor, arrangement_factor):
        for b in bands:
            if b.re_min <= reynolds_gap < b.re_max:
                return (
                    b.c
                    * reynolds_gap**b.m
                    * prandtl**b.n
                    * (prandtl / prandtl_surface)**0.25
                    * row_factor
                    * arrangement_factor
                )
        raise ValueError(f"Re={reynolds_gap} not covered by any correlation band")

    # ---------- engines ----------
    def _h_singlephase(self, zone: Zone):
        if zone in {"firepass", "smokepass", "economiser"}:
            re  = self._pick(zone=zone, base="reynolds_gap")
            pr  = self._pick(zone=zone, base="prandtl")
            prs = self._pick(zone=zone, base="prandtl_surface")
            bands = self._pick(zone=zone, base="zukauskas_bands")  # List[BankBand]
            rf = self._pick(zone=zone, base="row_correction_factor")
            af = self._pick(zone=zone, base="arrangement_correction_factor")
            nu = self._nusselt_bank_banded(re, pr, prs, bands, rf, af)
            k  = self._pick(zone=zone, base="thermal_conductivity")
            d  = self._pick(zone=zone, base="characteristic_diameter")
            return self.heat_transfer_coefficient_from_nusselt(nu, k, d)

        if zone == "reversal":
            re  = self._pick(zone=zone, base="reynolds") if hasattr(self, "reynolds_reversal") else self.reynolds_reversal
            pr  = self._pick(zone=zone, base="prandtl")  if hasattr(self, "prandtl_reversal")  else self.prandtl_reversal
            d_h = self._pick(zone=zone, base="diameter") if hasattr(self, "diameter_reversal") else self.diameter_reversal
            rc  = self._pick(zone=zone, base="curvature_radius") if hasattr(self, "curvature_radius_reversal") else self.curvature_radius_reversal
            a   = self._pick(zone=zone, base="coefficient_a")
            b   = self._pick(zone=zone, base="exponent_b")
            nu  = self.nusselt_churchill_bernstein_with_curvature(re, pr, d_h, rc, a, b)
            k   = self._pick(zone=zone, base="thermal_conductivity") if hasattr(self, "thermal_conductivity_reversal") else self.thermal_conductivity_reversal
            d_c = self._pick(zone=zone, base="characteristic_diameter") if hasattr(self, "characteristic_diameter_reversal") else self.characteristic_diameter_reversal
            return self.heat_transfer_coefficient_from_nusselt(nu, k, d_c)

        raise ValueError(f"Unknown zone: {zone}")

    def _h_twophase(self, zone: Zone):
        h_sp = self._h_singlephase(zone)
        h_nb = self._pick(zone=zone, base="h_nucleate")
        G    = self._pick(zone=zone, base="mass_flux")
        qpp  = self._pick(zone=zone, base="heat_flux")
        h_lv = self._pick(zone=zone, base="latent_heat")
        x    = self._pick(zone=zone, base="vapor_quality")
        rho_l= self._pick(zone=zone, base="liquid_density")
        rho_v= self._pick(zone=zone, base="vapor_density")
        mu_l = self._pick(zone=zone, base="liquid_viscosity")
        mu_v = self._pick(zone=zone, base="vapor_viscosity")
        sigma= self._pick(zone=zone, base="surface_tension")
        d    = self._pick(zone=zone, base="diameter")
        return self.gungor_winterton_boiling(h_sp, h_nb, G, qpp, h_lv, x, rho_l, rho_v, mu_l, mu_v, sigma, d)
