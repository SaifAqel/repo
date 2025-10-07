from __future__ import annotations
from pathlib import Path
from typing import Any, Dict, Optional
import yaml

from common.units import ureg, Q_
# import your dataclasses here (assumes they are in the same module or adjust import)
from heat_transfer.config.models import (
    Wall, Surface, Surfaces, DrumGeometry, PassGeometry, ReversalGeometry,
    Nozzle, ReversalNozzles, Drum, Pass, Reversal, Stages,
    GasStream, Water, Environment, Config
)


class ConfigLoader:
    @staticmethod
    def _qty(node: Dict[str, Any]) -> Q_:
        """
        Convert a YAML node { value: X, unit: "..." } to a pint Quantity Q_.
        If unit is "-" or empty treat as dimensionless.
        """
        if node is None:
            raise ValueError("Expected {{ value, unit }} node, got None")
        value = node.get("value")
        unit = node.get("unit", "")
        if unit is None or unit.strip() == "" or unit.strip() == "-":
            return Q_(value)
        return Q_(value, ureg(unit))

    @classmethod
    def _build_wall(cls, node: Dict[str, Any]) -> Wall:
        props = node.get("properties", node)  # handle either wall: properties: ... or wall: { ... }
        return Wall(
            thickness=cls._qty(props["thickness"]),
            conductivity=cls._qty(props["conductivity"]),
        )

    @classmethod
    def _build_surface(cls, node: Dict[str, Any]) -> Surface:
        return Surface(
            roughness=cls._qty(node["roughness"]),
            emissivity=cls._qty(node["emissivity"]),
            fouling_thickness=cls._qty(node["fouling_thickness"]),
            fouling_conductivity=cls._qty(node["fouling_conductivity"]),
        )

    @classmethod
    def _build_surfaces(cls, node: Dict[str, Any]) -> Surfaces:
        return Surfaces(
            inner=cls._build_surface(node["inner"]),
            outer=cls._build_surface(node["outer"]),
        )

    @classmethod
    def _build_drum_geometry(cls, node: Dict[str, Any]) -> DrumGeometry:
        return DrumGeometry(
            inner_diameter=cls._qty(node["inner_diameter"]),
            inner_length=cls._qty(node["inner_length"]),
            wall=cls._build_wall(node["wall"]),
        )

    @classmethod
    def _build_pass_geometry(cls, node: Dict[str, Any]) -> PassGeometry:
        # number_of_tubes may be missing for drum-like objects; expect it's present for passes
        return PassGeometry(
            inner_diameter=cls._qty(node["inner_diameter"]),
            inner_length=cls._qty(node["inner_length"]),
            number_of_tubes=cls._qty(node["number_of_tubes"]),
            layout=node.get("layout", ""),
            pitch=cls._qty(node["pitch"]),
            wall=cls._build_wall(node["wall"]),
        )

    @classmethod
    def _build_reversal_geometry(cls, node: Dict[str, Any]) -> ReversalGeometry:
        return ReversalGeometry(
            inner_diameter=cls._qty(node["inner_diameter"]),
            inner_length=cls._qty(node["inner_length"]),
            wall=cls._build_wall(node["wall"]),
        )

    @classmethod
    def _build_nozzle(cls, node: Dict[str, Any]) -> Nozzle:
        return Nozzle(
            diameter=cls._qty(node["diameter"]),
            length=cls._qty(node["length"]),
        )

    @classmethod
    def _build_reversal_nozzles(cls, node: Dict[str, Any]) -> ReversalNozzles:
        return ReversalNozzles(
            inlet=cls._build_nozzle(node["inlet"]),
            outlet=cls._build_nozzle(node["outlet"]),
        )

    @classmethod
    def _build_drum(cls, node: Dict[str, Any]) -> Drum:
        return Drum(
            geometry=cls._build_drum_geometry(node["geometry"]),
            surfaces=cls._build_surfaces(node["surfaces"]),
        )

    @classmethod
    def _build_pass(cls, node: Dict[str, Any]) -> Pass:
        return Pass(
            geometry=cls._build_pass_geometry(node["geometry"]),
            surfaces=cls._build_surfaces(node["surfaces"]),
        )

    @classmethod
    def _build_reversal(cls, node: Dict[str, Any]) -> Reversal:
        return Reversal(
            geometry=cls._build_reversal_geometry(node["geometry"]),
            nozzles=cls._build_reversal_nozzles(node["nozzles"]),
            surfaces=cls._build_surfaces(node["surfaces"]),
        )

    @classmethod
    def _build_stages(cls, node: Dict[str, Any]) -> Stages:
        return Stages(
            drum=cls._build_drum(node["drum"]),
            pass1=cls._build_pass(node["pass1"]),
            reversal1=cls._build_reversal(node["reversal1"]),
            pass2=cls._build_pass(node["pass2"]),
            reversal2=cls._build_reversal(node["reversal2"]),
            pass3=cls._build_pass(node["pass3"]),
        )

    @classmethod
    def _build_gas_stream(cls, node: Dict[str, Any]) -> GasStream:
        composition = {k: cls._qty(v) for k, v in node.get("composition", {}).items()}
        spectro = {k: cls._qty(v) for k, v in node.get("spectroscopic_data", {}).items()}
        return GasStream(
            mass_flow_rate=cls._qty(node["mass_flow_rate"]),
            temperature=cls._qty(node["temperature"]),
            pressure=cls._qty(node["pressure"]),
            composition=composition,
            spectroscopic_data=spectro,
            z=cls._qty(node["z"])
        )

    @classmethod
    def _build_water(cls, node: Dict[str, Any]) -> Water:
        composition = {k: cls._qty(v) for k, v in node.get("composition", {}).items()}
        return Water(
            mass_flow_rate=cls._qty(node["mass_flow_rate"]),
            temperature=cls._qty(node["temperature"]),
            pressure=cls._qty(node["pressure"]),
            composition=composition,
            z=cls._qty(node["z"])

        )

    @classmethod
    def _build_environment(cls, node: Dict[str, Any]) -> Environment:
        return Environment(
            ambient_temperature=cls._qty(node["ambient_temperature"]),
            radiative_temperature=cls._qty(node["radiative_temperature"]),
            external_emissivity=cls._qty(node["external_emissivity"]),
            external_h=cls._qty(node["external_h"]),
            radiation_view_factor_external=cls._qty(node["radiation_view_factor_external"]),
        )

    @classmethod
    def load_from_path(cls, path: str | Path) -> Config:
        path = Path(path)
        with path.open("r", encoding="utf-8") as fh:
            doc = yaml.safe_load(fh)
        return cls.load_from_dict(doc)

    @classmethod
    def load_from_string(cls, s: str) -> Config:
        doc = yaml.safe_load(s)
        return cls.load_from_dict(doc)

    @classmethod
    def load_from_dict(cls, doc: Dict[str, Any]) -> Config:
        return Config(
            gas_inlet=cls._build_gas_stream(doc["gas_inlet"]),
            water_inlet=cls._build_water(doc["water_inlet"]),
            environment=cls._build_environment(doc["environment"]),
            stages=cls._build_stages(doc["stages"]),
        )
