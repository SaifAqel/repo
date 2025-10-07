# schema.py
from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, Any, Mapping, Optional
import fnmatch
import tomllib
from common.units import ureg, Q_


def _to_float(v: Any) -> float:
    # No validation. Best-effort numeric cast.
    if isinstance(v, (int, float)):
        return float(v)
    if isinstance(v, str):
        return float(v.strip())
    # Fallback for types like Decimal
    return float(v)


def _quantity(v: Any, unit: str) -> Q_:
    return ureg.Quantity(_to_float(v), unit)


# -------------------------
# Data Model (mirrors TOML)
# -------------------------
@dataclass
class Wall:
    thickness: Q_
    conductivity: Q_


@dataclass
class Surface:
    roughness: Q_
    emissivity: Q_
    fouling_thickness: Q_
    fouling_conductivity: Q_


@dataclass
class Surfaces:
    inner: Surface
    outer: Surface


@dataclass
class DrumGeometry:
    inner_diameter: Q_
    inner_length: Q_
    wall: Wall


@dataclass
class PassGeometry:
    inner_diameter: Q_
    inner_length: Q_
    number_of_tubes: Q_
    layout: str
    pitch: Q_
    wall: Wall


@dataclass
class ReversalGeometry:
    inner_diameter: Q_
    inner_length: Q_
    wall: Wall


@dataclass
class Nozzle:
    diameter: Q_
    length: Q_


@dataclass
class ReversalNozzles:
    inlet: Nozzle
    outlet: Nozzle


@dataclass
class Drum:
    geometry: DrumGeometry
    surfaces: Surfaces


@dataclass
class Pass:
    geometry: PassGeometry
    surfaces: Surfaces


@dataclass
class Reversal:
    geometry: ReversalGeometry
    nozzles: ReversalNozzles
    surfaces: Surfaces


@dataclass
class Stages:
    drum: Drum
    pass1: Pass
    reversal1: Reversal
    pass2: Pass
    reversal2: Reversal
    pass3: Pass


@dataclass
class GasSide:
    mass_flow_rate: Q_
    inlet_temperature: Q_
    inlet_pressure: Q_
    composition: Dict[str, Q_]
    spectroscopic_data: Dict[str,Q_]


@dataclass
class WaterSide:
    mass_flow_rate: Q_
    inlet_temperature: Q_
    outlet_temperature: Q_
    pressure: Q_
    mu_exp: Q_
    C_sf: Q_
    n: Q_
    composition: Dict[str, Q_]


@dataclass
class Environment:
    ambient_temperature: Q_
    radiative_temperature: Q_
    external_emissivity: Q_
    external_h: Q_
    radiation_view_factor_external: Q_


@dataclass
class Config:
    gas_side: GasSide
    water_side: WaterSide
    environment: Environment
    stages: Stages


# -------------------------
# Loader
# -------------------------
def _dot_join(parent: Optional[str], key: str) -> str:
    return key if not parent else f"{parent}.{key}"


def _resolve_unit(path: str, unit_map: Mapping[str, str]) -> Optional[str]:
    # Exact match first
    if path in unit_map:
        return unit_map[path]
    # Wildcard patterns
    for pattern, unit in unit_map.items():
        if "*" in pattern or "?" in pattern or "[" in pattern:
            if fnmatch.fnmatch(path, pattern):
                return unit
    return None


def _map_leaf_as_quantity(
    path: str, value: Any, unit_map: Mapping[str, str]
) -> Any:
    unit = _resolve_unit(path, unit_map)
    if unit is None:
        return value
    return _quantity(value, unit)


def _load_wall(d: Mapping[str, Any], base: str, um: Mapping[str, str]) -> Wall:
    return Wall(
        thickness=_map_leaf_as_quantity(f"{base}.thickness", d["thickness"], um),
        conductivity=_map_leaf_as_quantity(f"{base}.conductivity", d["conductivity"], um),
    )


def _load_surface(d: Mapping[str, Any], base: str, um: Mapping[str, str]) -> Surface:
    return Surface(
        roughness=_map_leaf_as_quantity(f"{base}.roughness", d["roughness"], um),
        emissivity=_map_leaf_as_quantity(f"{base}.emissivity", d["emissivity"], um),
        fouling_thickness=_map_leaf_as_quantity(f"{base}.fouling_thickness", d["fouling_thickness"], um),
        fouling_conductivity=_map_leaf_as_quantity(f"{base}.fouling_conductivity", d["fouling_conductivity"], um),
    )


def _load_surfaces(d: Mapping[str, Any], base: str, um: Mapping[str, str]) -> Surfaces:
    return Surfaces(
        inner=_load_surface(d["inner"], f"{base}.inner", um),
        outer=_load_surface(d["outer"], f"{base}.outer", um),
    )


def _load_drum_geometry(d: Mapping[str, Any], base: str, um: Mapping[str, str]) -> DrumGeometry:
    return DrumGeometry(
        inner_diameter=_map_leaf_as_quantity(f"{base}.inner_diameter", d["inner_diameter"], um),
        inner_length=_map_leaf_as_quantity(f"{base}.inner_length", d["inner_length"], um),
        wall=_load_wall(d["wall"]["properties"], f"{base}.wall.properties", um),
    )


def _load_pass_geometry(d: Mapping[str, Any], base: str, um: Mapping[str, str]) -> PassGeometry:
    return PassGeometry(
        inner_diameter=_map_leaf_as_quantity(f"{base}.inner_diameter", d["inner_diameter"], um),
        inner_length=_map_leaf_as_quantity(f"{base}.inner_length", d["inner_length"], um),
        number_of_tubes=_map_leaf_as_quantity(f"{base}.number_of_tubes", d["number_of_tubes"], um),
        layout=str(d["layout"]),
        pitch=_map_leaf_as_quantity(f"{base}.pitch", d["pitch"], um),
        wall=_load_wall(d["wall"]["properties"], f"{base}.wall.properties", um),
    )


def _load_reversal_geometry(d: Mapping[str, Any], base: str, um: Mapping[str, str]) -> ReversalGeometry:
    return ReversalGeometry(
        inner_diameter=_map_leaf_as_quantity(f"{base}.inner_diameter", d["inner_diameter"], um),
        inner_length=_map_leaf_as_quantity(f"{base}.inner_length", d["inner_length"], um),
        wall=_load_wall(d["wall"]["properties"], f"{base}.wall.properties", um),
    )


def _load_nozzle(d: Mapping[str, Any], base: str, um: Mapping[str, str]) -> Nozzle:
    return Nozzle(
        diameter=_map_leaf_as_quantity(f"{base}.diameter", d["diameter"], um),
        length=_map_leaf_as_quantity(f"{base}.length", d["length"], um),
    )


def _load_reversal_nozzles(d: Mapping[str, Any], base: str, um: Mapping[str, str]) -> ReversalNozzles:
    return ReversalNozzles(
        inlet=_load_nozzle(d["inlet"], f"{base}.inlet", um),
        outlet=_load_nozzle(d["outlet"], f"{base}.outlet", um),
    )


def _load_drum(d: Mapping[str, Any], base: str, um: Mapping[str, str]) -> Drum:
    return Drum(
        geometry=_load_drum_geometry(d["geometry"], f"{base}.geometry", um),
        surfaces=_load_surfaces(d["surfaces"], f"{base}.surfaces", um),
    )


def _load_pass(d: Mapping[str, Any], base: str, um: Mapping[str, str]) -> Pass:
    return Pass(
        geometry=_load_pass_geometry(d["geometry"], f"{base}.geometry", um),
        surfaces=_load_surfaces(d["surfaces"], f"{base}.surfaces", um),
    )


def _load_reversal(d: Mapping[str, Any], base: str, um: Mapping[str, str]) -> Reversal:
    return Reversal(
        geometry=_load_reversal_geometry(d["geometry"], f"{base}.geometry", um),
        nozzles=_load_reversal_nozzles(d["nozzles"], f"{base}.nozzles", um),
        surfaces=_load_surfaces(d["surfaces"], f"{base}.surfaces", um),
    )


def _load_composition(d: Mapping[str, Any], base: str, um: Mapping[str, str]) -> Dict[str, Q_]:
    out: Dict[str, Q_] = {}
    for k, v in d.items():
        path = f"{base}.{k}"
        out[k] = _map_leaf_as_quantity(path, v, um)
    return out

def _load_spectroscopic(d: Mapping[str, Any], base: str, um: Mapping[str, str]) -> Dict[str, Q_]:
    out: Dict[str, Q_] = {}
    for k, v in d.items():
        path = f"{base}.{k}"
        out[k] = _map_leaf_as_quantity(path, v, um)
    return out


def _load_gas_side(d: Mapping[str, Any], base: str, um: Mapping[str, str]) -> GasSide:
    return GasSide(
        mass_flow_rate=_map_leaf_as_quantity(f"{base}.mass_flow_rate", d["mass_flow_rate"], um),
        inlet_temperature=_map_leaf_as_quantity(f"{base}.inlet_temperature", d["inlet_temperature"], um),
        inlet_pressure=_map_leaf_as_quantity(f"{base}.inlet_pressure", d["inlet_pressure"], um),
        composition=_load_composition(d.get("composition", {}), f"{base}.composition", um),
        spectroscopic_data=_load_spectroscopic(d.get("spectroscopic_data", {}), f"{base}.spectroscopic_data", um),
    )


def _load_water_side(d: Mapping[str, Any], base: str, um: Mapping[str, str]) -> WaterSide:
    return WaterSide(
        mass_flow_rate=_map_leaf_as_quantity(f"{base}.mass_flow_rate", d["mass_flow_rate"], um),
        inlet_temperature=_map_leaf_as_quantity(f"{base}.inlet_temperature", d["inlet_temperature"], um),
        outlet_temperature=_map_leaf_as_quantity(f"{base}.outlet_temperature", d["outlet_temperature"], um),
        pressure=_map_leaf_as_quantity(f"{base}.pressure", d["pressure"], um),
        mu_exp = _map_leaf_as_quantity(f"{base}.mu_exp", d["mu_exp"], um),
        C_sf = _map_leaf_as_quantity(f"{base}.C_sf", d["C_sf"], um),
        n = _map_leaf_as_quantity(f"{base}.n", d["n"], um),
        composition=_load_composition(d.get("composition", {}), f"{base}.composition", um),

    )


def _load_environment(d: Mapping[str, Any], base: str, um: Mapping[str, str]) -> Environment:
    return Environment(
        ambient_temperature=_map_leaf_as_quantity(f"{base}.ambient_temperature", d["ambient_temperature"], um),
        radiative_temperature=_map_leaf_as_quantity(f"{base}.radiative_temperature", d["radiative_temperature"], um),
        external_emissivity=_map_leaf_as_quantity(f"{base}.external_emissivity", d["external_emissivity"], um),
        external_h=_map_leaf_as_quantity(f"{base}.external_h", d["external_h"], um),
        radiation_view_factor_external=_map_leaf_as_quantity(f"{base}.radiation_view_factor_external", d["radiation_view_factor_external"], um),
    )


def _load_stages(d: Mapping[str, Any], base: str, um: Mapping[str, str]) -> Stages:
    return Stages(
        drum=_load_drum(d["drum"], f"{base}.drum", um),
        pass1=_load_pass(d["pass1"], f"{base}.pass1", um),
        reversal1=_load_reversal(d["reversal1"], f"{base}.reversal1", um),
        pass2=_load_pass(d["pass2"], f"{base}.pass2", um),
        reversal2=_load_reversal(d["reversal2"], f"{base}.reversal2", um),
        pass3=_load_pass(d["pass3"], f"{base}.pass3", um),
    )


def _flatten_units(d: Mapping[str, Any], prefix: Optional[str] = None) -> Dict[str, str]:
    flat: Dict[str, str] = {}
    for k, v in d.items():
        path = _dot_join(prefix, k)
        if isinstance(v, dict):
            flat.update(_flatten_units(v, path))
        else:
            flat[path] = str(v)
    return flat


def load_config(config_path: str, units_path: str) -> Config:
    with open(config_path, "rb") as f:
        config_data = tomllib.load(f)
    with open(units_path, "rb") as f:
        units_data = tomllib.load(f)

    unit_map = _flatten_units(units_data)

    return Config(
        gas_side=_load_gas_side(config_data["gas_side"], "gas_side", unit_map),
        water_side=_load_water_side(config_data["water_side"], "water_side", unit_map),
        environment=_load_environment(config_data["environment"], "environment", unit_map),
        stages=_load_stages(config_data["stages"], "stages", unit_map),
    )
