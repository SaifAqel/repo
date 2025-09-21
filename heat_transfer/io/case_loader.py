# io/case_loader.py
from __future__ import annotations

import tomllib
from typing import Dict, Callable

from heat_transfer.thermo.provider import CoolPropThermoProvider
from heat_transfer.geometry.tube import SingleTubeGeom
from heat_transfer.geometry.tube_bank import TubeBankGeom
from heat_transfer.geometry.reversal import ReversalChamberGeom
from heat_transfer.models.axial_marcher import AxialMarcher
from heat_transfer.models.stage_runner import StageRunner
from heat_transfer.transport.htc_gas import (
    GnielinskiModel,
    DittusBoelterModel,
    ChurchillBernsteinModel,
)
from heat_transfer.transport.emissivity import (
    LecknerModel,
    SmithModel,
    WSGGM_SimpleModel,
)
from heat_transfer.transport.htc_water_boil import (
    ChenBoilingModel,
    ThomModel,
    GungorWintertonModel,
)
from heat_transfer.flow.friction import ColebrookWhite, Churchill, Haaland
from heat_transfer.core.types import GasState


def _pick(model_name: str, table: Dict[str, Callable]):
    if model_name not in table:
        raise ValueError(model_name)
    return table[model_name]()


def _get(d: dict, *keys, default=None):
    """Fetch first present key in d from keys."""
    for k in keys:
        if k in d:
            return d[k]
    return default


def load_case(settings_toml_path: str):
    # parse settings.toml
    with open(settings_toml_path, "rb") as f:
        cfg = tomllib.load(f)

    # property backend
    backend = (
        cfg.get("property_backend", {})
           .get("coolprop", {})
           .get("backend")
        or cfg.get("property_backend", {}).get("backend")
        or "HEOS"
    )
    thermo = CoolPropThermoProvider(backend)

    # correlation registries
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
    fr_map = {"colebrook_white": ColebrookWhite, "churchill": Churchill, "haaland": Haaland}

    # defaults
    gnum = cfg.get("global_numerics", {})
    N_cells_default = int(gnum.get("N_cells_default", 50))
    defaults = cfg.get("defaults", {})
    defaults_foul = defaults.get("fouling", {})
    R_fg_default = defaults_foul.get("R_fg", 0.0)
    R_fw_default = defaults_foul.get("R_fw", 0.0)

    # geometry builders
    geom_builders = {
        "single_tube": lambda g: SingleTubeGeom(**g),
        "tube_bank": lambda g: TubeBankGeom(**g),
        "reversal_chamber": lambda g: ReversalChamberGeom(**g),
    }

    # inlet gas (units as per TOML keys)
    gi = cfg.get("inlet_gas", {})
    T_in = _get(gi, "T_in_K", "T_in", default=None)
    P_in = _get(gi, "P_in_Pa", "P_in", default=None)
    m_dot = _get(gi, "m_dot_kg_s", "m_dot", default=None)
    composition = gi.get("composition", {})
    if T_in is None or P_in is None or m_dot is None:
        raise ValueError("inlet_gas requires T_in_K, P_in_Pa, and m_dot_kg_s")

    # thermophysical properties at inlet
    props = thermo.gas_props(T_in, P_in, composition)

    # helper: first stage flow area to compute inlet velocity
    if not cfg.get("stages"):
        raise ValueError("settings.toml must include at least one [[stages]] entry")

    first_stage = cfg["stages"][0]
    first_geom = geom_builders[first_stage["type"]](first_stage["geometry"])
    A0 = first_geom.flow_area()
    rho0 = props["rho"]
    u0 = m_dot / (rho0 * A0)

    def gas_inlet():
        return GasState(
            T_in,
            P_in,
            m_dot,
            composition,
            rho0,
            props["cp"],
            props["mu"],
            props["k"],
            props["Pr"],
            u0,
        )

    # drum saturation temperature (use first drum)
    if not cfg.get("drums"):
        raise ValueError("settings.toml must include at least one [[drums]] entry")
    drum0 = cfg["drums"][0]
    P_drum = _get(drum0, "pressure_Pa", "pressure", default=None)
    if P_drum is None:
        raise ValueError("drums[].pressure_Pa is required")
    Tsat = thermo.water_saturation(P_drum)["Tsat"]

    # build stages
    stages = []
    for s in cfg["stages"]:
        stype = s["type"]
        geom = geom_builders[stype](s["geometry"])
        flow = s.get("flow", "co_current")

        # correlations with sane defaults if omitted
        corr = s.get("correlations", {})
        h_gas_name = corr.get("h_gas", "gnielinski")
        eps_name = corr.get("epsilon_gas", "leckner")
        h_boil_name = corr.get("h_water_boiling", "chen")
        fr_name = corr.get("friction", "churchill")

        # discretization
        N = int(s.get("discretization", {}).get("N_cells", N_cells_default))

        # fouling (accept either unit-tagged or bare)
        f = s.get("fouling", {})
        R_fg = _get(f, "R_fg", "R_fg_m2K_W", default=R_fg_default)
        R_fw = _get(f, "R_fw", "R_fw_m2K_W", default=R_fw_default)

        losses = s.get("losses", [])

        marcher = AxialMarcher(
            geom,
            thermo,
            _pick(h_gas_name, htc_gas_map),
            _pick(eps_name, eps_map),
            _pick(h_boil_name, boil_map),
            _pick(fr_name, fr_map),
            geom.length(),
            N,
            flow,
            R_fg,
            R_fw,
            losses,
        )
        stages.append({"runner": StageRunner(geom, marcher), "Tsat": Tsat, "G_water": 1.0})

    # return unified config and constructed callables
    return cfg, {"gas_inlet": gas_inlet, "stages": stages}
