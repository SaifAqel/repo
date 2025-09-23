# io/case_loader.py
from __future__ import annotations

import tomllib
from typing import Dict, Callable, List, Tuple, Optional

from common.units import ureg, Q_

from heat_transfer1.thermo.provider import CoolPropThermoProvider
from heat_transfer1.geometry.tube import SingleTubeGeom
from heat_transfer1.geometry.tube_bank import TubeBankGeom
from heat_transfer1.geometry.reversal import ReversalChamberGeom
from heat_transfer1.models.axial_marcher import AxialMarcher
from heat_transfer1.models.stage_runner import StageRunner
from heat_transfer1.transport.htc_gas import (
    GnielinskiModel,
    DittusBoelterModel,
    ChurchillBernsteinModel,
)
from heat_transfer1.transport.emissivity import (
    LecknerModel,
    SmithModel,
    WSGGM_SimpleModel,
)
from heat_transfer1.transport.htc_water_boil import (
    ChenBoilingModel,
    ThomModel,
    GungorWintertonModel,
)
from heat_transfer1.flow.friction import ColebrookWhite, Churchill, Haaland
from heat_transfer1.core.types import GasState


# ------------ helpers ------------

def _pick(model_name: str, table: Dict[str, Callable]):
    if model_name not in table:
        raise ValueError(f"unknown model: {model_name}")
    return table[model_name]()

def _get(d: dict, *keys, default=None):
    for k in keys:
        if k in d:
            return d[k]
    return default

def _q(val, unit: str) -> Optional[Q_]:
    if val is None:
        return None
    return Q_(val, unit)

def _q_dimless(val) -> Q_:
    return Q_(val, "dimensionless")

def _q_list(vals, unit: str) -> Optional[List[Q_]]:
    if vals is None:
        return None
    return [Q_(v, unit) for v in vals]

def _map_geom_keys(stype: str, g: dict) -> dict:
    """Translate legacy *_m, *_m2, *_m2K_W keys to new names and wrap with units."""
    if stype == "single_tube":
        return {
            "L": _q(_get(g, "L", "L_m"), "m"),
            "Di": _q(_get(g, "Di", "Di_m"), "m"),
            "Do": _q(_get(g, "Do", "Do_m"), "m"),
            "roughness": _q(_get(g, "roughness", "roughness_m"), "m"),
            "thickness": _q(_get(g, "thickness", "thickness_m"), "m"),
            "R_fg": _q(_get(g, "R_fg", "R_fg_m2K_W"), "m^2*K/W"),
            "R_fw": _q(_get(g, "R_fw", "R_fw_m2K_W"), "m^2*K/W"),
        }
    if stype == "tube_bank":
        return {
            "L": _q(_get(g, "L", "L_m"), "m"),
            "Di": _q(_get(g, "Di", "Di_m"), "m"),
            "Do": _q(_get(g, "Do", "Do_m"), "m"),
            "Ntubes": g.get("Ntubes"),
            "pitch": _q(_get(g, "pitch", "pitch_m"), "m"),
            "arrangement": g.get("arrangement"),
            "roughness": _q(_get(g, "roughness", "roughness_m"), "m"),
            "thickness": _q(_get(g, "thickness", "thickness_m"), "m"),
            "frontal_width": _q(_get(g, "frontal_width", "frontal_width_m"), "m"),
            "frontal_height": _q(_get(g, "frontal_height", "frontal_height_m"), "m"),
        }
    if stype == "reversal_chamber":
        return {
            "L_eq": _q(_get(g, "L_eq", "L_eq_m"), "m"),
            "area_flow": _q(_get(g, "area_flow", "area_flow_m2"), "m^2"),
            "wetted_perimeter": _q(_get(g, "wetted_perimeter", "wetted_perimeter_m"), "m"),
            "ht_perimeter": _q(_get(g, "ht_perimeter", "area_ht_per_length_m"), "m"),
        }
    raise ValueError(f"unknown geometry type: {stype}")

def _build_geom(stype: str, g: dict):
    if stype == "single_tube":
        return SingleTubeGeom(**_map_geom_keys(stype, g))
    if stype == "tube_bank":
        return TubeBankGeom(**_map_geom_keys(stype, g))
    if stype == "reversal_chamber":
        return ReversalChamberGeom(**_map_geom_keys(stype, g))
    raise ValueError(stype)

def _wrap_losses(losses_raw: List[Tuple[str, float | int | str]]) -> List[Tuple[str, Q_]]:
    # Accept [("elbow", 0.9), ...] and wrap to dimensionless Quantity
    return [(name, _q_dimless(k)) for name, k in losses_raw]


# ------------ loader ------------

def load_case(settings_toml_path: str):
    with open(settings_toml_path, "rb") as f:
        cfg = tomllib.load(f)

    backend = (
        cfg.get("property_backend", {}).get("coolprop", {}).get("backend")
        or cfg.get("property_backend", {}).get("backend")
        or "HEOS"
    )
    thermo = CoolPropThermoProvider(backend)

    htc_gas_map = {
        "gnielinski": GnielinskiModel,
        "dittus_boelter": DittusBoelterModel,
        "churchill_bernstein": ChurchillBernsteinModel,
    }
    eps_map = {
        "leckner": LecknerModel,
        "smith": SmithModel,
        "wsggm_simple": WSGGM_SimpleModel,
    }
    boil_map = {
        "chen": ChenBoilingModel,
        "thom": ThomModel,
        "gungor_winterton": GungorWintertonModel,
    }
    fr_map = {
        "colebrook_white": ColebrookWhite,
        "churchill": Churchill,
        "haaland": Haaland,
    }

    gnum = cfg.get("global_numerics", {})
    N_cells_default = int(gnum.get("N_cells_default", 50))

    defaults = cfg.get("defaults", {})
    defaults_foul = defaults.get("fouling", {})
    R_fg_default = _q(_get(defaults_foul, "R_fg", "R_fg_m2K_W", default=0.0), "m^2*K/W")
    R_fw_default = _q(_get(defaults_foul, "R_fw", "R_fw_m2K_W", default=0.0), "m^2*K/W")

    if not cfg.get("stages"):
        raise ValueError("settings.toml must include at least one [[stages]] entry")
    if not cfg.get("drums"):
        raise ValueError("settings.toml must include at least one [[drums]] entry")

    # Inlet gas
    gi = cfg.get("inlet_gas", {})
    T_in = _q(_get(gi, "T_in", "T_in_K"), "K")
    P_in = _q(_get(gi, "P_in", "P_in_Pa"), "Pa")
    m_dot = _q(_get(gi, "m_dot", "m_dot_kg_s"), "kg/s")
    if T_in is None or P_in is None or m_dot is None:
        raise ValueError("inlet_gas requires T_in_K, P_in_Pa, and m_dot_kg_s")

    composition_raw = gi.get("composition", {})
    composition = {sp: _q_dimless(v) for sp, v in composition_raw.items()}

    props = thermo.gas_props(T_in, P_in, composition)

    # Geometry of first stage for velocity
    first_stage = cfg["stages"][0]
    first_geom = _build_geom(first_stage["type"], first_stage["geometry"])
    A0 = first_geom.flow_area()
    rho0 = props["rho"]
    u0 = (m_dot / (rho0 * A0)).to("m/s")

    def gas_inlet() -> GasState:
        return GasState(
            T=T_in,
            P=P_in,
            m_dot=m_dot,
            y=composition,
            rho=rho0,
            cp=props["cp"],
            mu=props["mu"],
            k=props["k"],
            Pr=props["Pr"],
            u=u0,
        )

    # Drum saturation temperature from first drum unless overridden
    drum0 = cfg["drums"][0]
    P_drum = _q(_get(drum0, "pressure", "pressure_Pa"), "Pa")
    if P_drum is None:
        raise ValueError("drums[].pressure_Pa is required")
    Tsat_global = thermo.water_saturation(P_drum)["Tsat"]

    # Build stages
    stages = []
    for s in cfg["stages"]:
        stype = s["type"]
        geom = _build_geom(stype, s["geometry"])
        flow = s.get("flow", "co_current")

        corr = s.get("correlations", {})
        h_gas_name = corr.get("h_gas", "gnielinski")
        eps_name = corr.get("epsilon_gas", "leckner")
        h_boil_name = corr.get("h_water_boiling", "chen")
        fr_name = corr.get("friction", "churchill")

        N = int(s.get("discretization", {}).get("N_cells", N_cells_default))

        f = s.get("fouling", {})
        R_fg = _q(_get(f, "R_fg", "R_fg_m2K_W", default=R_fg_default.magnitude), "m^2*K/W")
        R_fw = _q(_get(f, "R_fw", "R_fw_m2K_W", default=R_fw_default.magnitude), "m^2*K/W")

        losses = _wrap_losses(s.get("losses", []))

        # Water-side inputs
        w = s.get("water", {})
        G_water = _q(_get(w, "G", "G_kg_m2_s"), "kg/m^2/s")
        if G_water is None:
            raise ValueError("stages[].water.G_kg_m2_s is required")
        Dh_water = _q(_get(s, "Dh_water", "Dh_water_m", default=None), "m")
        x_prof   = _q_list(s.get("x_water_profile", None), "dimensionless")
        qpp_prof = _q_list(s.get("qpp_init_profile", None), "W/m^2")

        # Optional per-stage Tsat override
        Tsat = _q(_get(s, "Tsat", "Tsat_K", default=None), "K") or Tsat_global

        marcher = AxialMarcher(
            geom=geom,
            thermo=thermo,
            htc_gas=_pick(h_gas_name, htc_gas_map),
            eps_model=_pick(eps_name, eps_map),
            htc_water=_pick(h_boil_name, boil_map),
            friction=_pick(fr_name, fr_map),
            L=geom.length(),
            N=N,
            flow=flow,
            R_fg=R_fg,
            R_fw=R_fw,
            losses=losses,
        )

        stages.append({
            "runner": StageRunner(geom, marcher),
            "Tsat": Tsat,
            "G_water": G_water,
            "Dh_water": Dh_water,
            "x_water_profile": x_prof,
            "qpp_init_profile": qpp_prof,
        })

    return cfg, {"gas_inlet": gas_inlet, "stages": stages}
