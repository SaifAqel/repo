from dataclasses import dataclass
from typing import Dict
import tomllib, pathlib
from common.units import ureg, Q_


@dataclass(frozen=True)
class Settings:
    species_molar_masses: Dict[str, Q_]
    species_cp_fluids_map: Dict[str, str]
    formation_enthalpies: Dict[str, Q_]
    latent_heat_H2O: Q_
    stoich_O2_per_mol: Dict[str, float]
    air_composition_mol: Dict[str, float]
    air_T: Q_
    air_P: Q_
    fuel_composition_mass: Dict[str, float]
    fuel_T: Q_
    fuel_P: Q_
    fuel_mass_flow: Q_
    T_ref: Q_
    ambient_T: Q_
    excess_air_ratio: float

def load_settings(path: str) -> Settings:
    s = tomllib.loads(pathlib.Path(path).read_text())
    return Settings(
        species_molar_masses={k: Q_(v, "kg/mol") for k, v in s["species"]["molar_masses"].items()},
        species_cp_fluids_map=s["species"]["cp_fluids_map"],
        formation_enthalpies={k: Q_(v, "kJ/mol") for k, v in s["thermo"]["formation_enthalpies"].items()},
        latent_heat_H2O=Q_(s["thermo"]["latent_heat_H2O_kJ_per_kg"], "kJ/kg"),
        stoich_O2_per_mol=s["stoich"]["O2_per_mol"],
        air_composition_mol=s["air"]["composition_mol"],
        air_T=Q_(s["air"]["T_C"], "degC"),
        air_P=Q_(s["air"]["P_Pa"], "Pa"),
        fuel_composition_mass=s["fuel"]["composition_mass"],
        fuel_T=Q_(s["fuel"]["T_C"], "degC"),
        fuel_P=Q_(s["fuel"]["P_Pa"], "Pa"),
        fuel_mass_flow=Q_(s["fuel"]["mass_flow_kg_s"], "kg/s"),
        T_ref=Q_(s["reference"]["T_ref"], "K"),
        ambient_T=Q_(s["reference"]["ambient_T"], "K"),
        excess_air_ratio=s["operation"]["excess_air_ratio"],
    )
