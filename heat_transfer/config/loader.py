from __future__ import annotations
from pathlib import Path
from typing import Any, Dict, Optional
import yaml
import cantera as ct
from common.units import ureg, Q_
# import your dataclasses here (assumes they are in the same module or adjust import)
from heat_transfer.config.models import (Wall, Surface, Surfaces, Nozzle, TubeGeometry, ReversalGeometry,
                                         Nozzles, ShellGeometry, FirePass, SmokePass, Reversal, BankGeometry,
                                         Economiser, Stages, GasStream, WaterStream, GasProps, WaterProps,
                                         EconomiserHot, EconomiserCold)

class ConfigLoader:
    @staticmethod
    def _qty(node: Dict[str, Any]) -> Q_:
        if node is None:
            raise ValueError("Expected {{ value, unit }} node, got None")
        value = node.get("value")
        unit = node.get("unit", "")
        if unit is None or unit.strip() == "" or unit.strip() == "-":
            return Q_(value)
        return Q_(value, ureg(unit))
    

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
    def _build_wall(cls, node: Dict[str, Any]) -> Wall:
        props = node.get("properties", node)  # handle either wall: properties: ... or wall: { ... }
        return Wall(
            thickness=cls._qty(props["thickness"]),
            conductivity=cls._qty(props["conductivity"]),
            surfaces=cls._build_surfaces(node["surfaces"])
        )

    @classmethod
    def _build_nozzle(cls, node: Dict[str, Any]) -> Nozzle:
        return Nozzle(
            k=cls._qty(node["k"])
        )
    
    @classmethod
    def _build_nozzles(cls, node: Dict[str, Any]) -> Nozzles:
        return Nozzles(
            inlet=cls._build_nozzle(node["inlet"]),
            outlet=cls._build_nozzle(node["outlet"]),
        )



    @classmethod
    def _build_tube_geometry(cls, node: Dict[str, Any]) -> TubeGeometry:
        # number_of_tubes may be missing for drum-like objects; expect it's present for passes
        return TubeGeometry(
            inner_diameter=cls._qty(node["inner_diameter"]),
            inner_length=cls._qty(node["inner_length"]),
            wall=cls._build_wall(node["wall"])
        )

    @classmethod
    def _build_bank_geometry(cls, node: Dict[str, Any]) -> BankGeometry:
        return BankGeometry(
            inner_diameter=cls._qty(node["inner_diameter"]),
            inner_length=cls._qty(node["inner_length"]),
            tubes_number=cls._qty(node["tubes_number"]),
            layout=node["layout"],
            pitch=cls._qty(node["pitch"]),
            wall=cls._build_wall(node["wall"])

        )

    @classmethod
    def _build_reversal_geometry(cls, node: Dict[str, Any]) -> ReversalGeometry:
        return ReversalGeometry(
            inner_length=cls._qty(node["inner_length"]),
            inner_diameter=cls._qty(node["inner_diameter"]),
            nozzles=cls._build_nozzles(node["nozzles"]),
            wall=cls._build_wall(node["wall"])

        )

    @classmethod
    def _build_economiser_hot(cls, node: Dict[str, Any]) -> EconomiserHot:
        return EconomiserHot(
            inner_length=cls._qty(node["inner_length"]),
            inner_diameter=cls._qty(node["inner_diameter"]),
            wall=cls._build_wall(node["wall"])
        )

    @classmethod
    def _build_economiser_cold(cls, node: Dict[str, Any]) -> EconomiserCold:
        return EconomiserCold(
            inner_length=cls._qty(node["inner_length"]),
            inner_diameter=cls._qty(node["inner_diameter"]),
            wall=cls._build_wall(node["wall"])
        )

    @classmethod
    def _build_shell_geometry(cls, node: Dict[str, Any]) -> ShellGeometry:
        return ShellGeometry(
            flow_area=cls._qty(node["flow_area"]),
            wetted_perimeter=cls._qty(node["wetted_perimeter"]),
            wall=cls._build_wall(node["wall"])
        )

    @classmethod
    def _build_fire_pass(cls, node: Dict[str, Any]) -> FirePass:
        return FirePass(
            hot_side=cls._build_tube_geometry(node["hot_side"]),
            cold_side=cls._build_shell_geometry(node["cold_side"]),
        )
    
    @classmethod
    def _build_smoke_pass(cls, node: Dict[str, Any]) -> SmokePass:
        return SmokePass(
            hot_side=cls._build_bank_geometry(node["hot_side"]),
            cold_side=cls._build_shell_geometry(node["cold_side"])
        )

    @classmethod
    def _build_reversal(cls, node: Dict[str, Any]) -> Reversal:
        return Reversal(
            hot_side=cls._build_reversal_geometry(node["hot_side"]),
            cold_side=cls._build_shell_geometry(node["cold_side"])
        )
    


        

    @classmethod
    def _build_economiser(cls, node: Dict[str, Any]) -> Economiser:
        return Economiser(
            hot_side=cls._build_economiser_hot(node["hot_side"]),
            cold_side=cls._build_economiser_cold(node["cold_side"])
        )

    @classmethod
    def _build_stages(cls, node: Dict[str, Any]) -> Stages:
        return Stages(
            HX_1=cls._build_fire_pass(node["HX_1"]),
            HX_2=cls._build_reversal(node["HX_2"]),
            HX_3=cls._build_smoke_pass(node["HX_3"]),
            HX_4=cls._build_reversal(node["HX_4"]),
            HX_5=cls._build_smoke_pass(node["HX_5"]),
            HX_6=cls._build_economiser(node["HX_6"])
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
            stage= None,
            gas_props=GasProps
        )

    @classmethod
    def _build_water(cls, node: Dict[str, Any]) -> WaterStream:
        composition = {k: cls._qty(v) for k, v in node.get("composition", {}).items()}
        return WaterStream(
            mass_flow_rate=cls._qty(node["mass_flow_rate"]),
            enthalpy=cls._qty(node["enthalpy"]),
            pressure=cls._qty(node["pressure"]),
            composition=composition,
            stage=None,
            water_props=WaterProps
        )

    @classmethod
    def load_stages(cls, path: str | Path) -> Stages:
        path = Path(path)
        with path.open("r", encoding="utf-8") as fh:
            doc = yaml.safe_load(fh)
        return cls._build_stages(doc["stages"])

    
    @classmethod
    def load_gas_stream(cls, path: str | Path) -> GasStream:
        path = Path(path)
        with path.open("r", encoding="utf-8") as fh:
            doc = yaml.safe_load(fh)
        return cls._build_gas_stream(doc["gas_stream"])

    
    @classmethod
    def load_water_stream(cls, path: str | Path) -> WaterStream:
        path = Path(path)
        with path.open("r", encoding="utf-8") as fh:
            doc = yaml.safe_load(fh)
        return cls._build_water(doc["water_stream"])

