# pint + TOML loader that returns your boiler config with SI units on every value,
# including dimensionless. Mirrors your Settings example.

from __future__ import annotations
from typing import Dict, Any, List
import tomllib, pathlib
from thermo.core.units import ureg, Q_  # same import style you used

def _q1(x: float) -> Any:
    return Q_(x, "1")

def load_boiler_config(path: str) -> Dict[str, Any]:
    s = tomllib.loads(pathlib.Path(path).read_text())

    # inlet_gas
    inlet_gas = {
        "T_in": Q_(s["inlet_gas"]["T_in_K"], "K"),
        "P_in": Q_(s["inlet_gas"]["P_in_Pa"], "Pa"),
        "m_dot": Q_(s["inlet_gas"]["m_dot_kg_s"], "kg/s"),
        "composition": {k: _q1(v) for k, v in s["inlet_gas"]["composition"].items()},
    }

    # drums
    drums: List[Dict[str, Any]] = []
    for d in s["drums"]:
        drums.append({
            "id": d["id"],
            "pressure": Q_(d["pressure_Pa"], "Pa"),
            "level": _q1(d["level"]),
            "geometry": {
                "diameter": Q_(d["geometry"]["diameter_m"], "m"),
                "length": Q_(d["geometry"]["length_m"], "m"),
                "nozzle_data": {
                    "steam_outlet": {"Dn": Q_(d["geometry"]["nozzle_data"]["steam_outlet"]["Dn_m"], "m")},
                    "water_inlet": {"Dn": Q_(d["geometry"]["nozzle_data"]["water_inlet"]["Dn_m"], "m")},
                },
            },
            "initial_water_inventory": Q_(d["initial_water_inventory_m3"], "m^3"),
        })

    # stages
    stages: List[Dict[str, Any]] = []
    for st in s["stages"]:
        base = {
            "type": st["type"],
            "flow": st["flow"],
            "correlations": st["correlations"],
            "losses": [{"element": L["element"], "K": _q1(L["K"])} for L in st.get("losses", [])],
            "fouling": {
                "R_fg": Q_(st["fouling"]["R_fg_m2K_W"], "m^2*K/W"),
                "R_fw": Q_(st["fouling"]["R_fw_m2K_W"], "m^2*K/W"),
            },
            "discretization": {"N_cells": _q1(st["discretization"]["N_cells"])},
        }

        geo = st["geometry"]
        if st["type"] == "tube_bank":
            geometry = {
                "L": Q_(geo["L_m"], "m"),
                "Di": Q_(geo["Di_m"], "m"),
                "Do": Q_(geo["Do_m"], "m"),
                "Ntubes": _q1(geo["Ntubes"]),
                "pitch": Q_(geo["pitch_m"], "m"),
                "arrangement": geo["arrangement"],
                "roughness": Q_(geo["roughness_m"], "m"),
                "thickness": Q_(geo["thickness_m"], "m"),
                "frontal_width": Q_(geo["frontal_width_m"], "m"),
                "frontal_height": Q_(geo["frontal_height_m"], "m"),
            }
        elif st["type"] == "reversal_chamber":
            geometry = {
                "L_eq": Q_(geo["L_eq_m"], "m"),
                "area_flow": Q_(geo["area_flow_m2"], "m^2"),
                "wetted_perimeter": Q_(geo["wetted_perimeter_m"], "m"),
                "area_ht_per_length": Q_(geo["area_ht_per_length_m"], "m"),  # A_ht / L
            }
        elif st["type"] == "single_tube":
            geometry = {
                "L": Q_(geo["L_m"], "m"),
                "Di": Q_(geo["Di_m"], "m"),
                "Do": Q_(geo["Do_m"], "m"),
                "roughness": Q_(geo["roughness_m"], "m"),
                "thickness": Q_(geo["thickness_m"], "m"),
            }
        else:
            raise ValueError(f"Unknown stage type: {st['type']}")

        stages.append({**base, "geometry": geometry})

    return {
        "inlet_gas": inlet_gas,
        "drums": drums,
        "stages": stages,
    }
