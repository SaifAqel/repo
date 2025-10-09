# loader.py
from __future__ import annotations

import yaml
from typing import Dict, Tuple, Any, Optional
import cantera as ct
from common.units import ureg, Q_
from heat_transfer.config.models import (
    Surface, Surfaces, Wall, Nozzle, Nozzles,
    TubeGeometry, ReversalChamberGeometry, BankLayout,
    ShellGeometry, Shell,
    FirePass, SmokePass, Reversal, HotSide, ColdSide, Economiser, Stages,
    WaterStream, GasStream,
)
from heat_transfer.functions.fluid_props import GasProps, WaterProps

# -----------------
# quantity helpers
# -----------------

def _q(node: Optional[Dict]) -> Optional[Q_]:
    if not node:
        return None
    v = node.get("value", None)
    u = node.get("unit", "")
    if v is None:
        return None
    return Q_(v, u)

def _qd(mapping: Optional[Dict[str, Dict]]) -> Dict[str, Q_]:
    return {k: _q(v) for k, v in (mapping or {}).items()}

# -----------------
# streams
# -----------------

def load_streams_from_yaml(yaml_path: str) -> Tuple[GasStream, WaterStream]:
    with open(yaml_path, "r", encoding="utf-8") as f:
        data = yaml.safe_load(f)

    g = data["gas_stream"]
    gas = GasStream(
        mass_flow_rate=_q(g.get("mass_flow_rate")),
        temperature=_q(g.get("temperature")),
        pressure=_q(g.get("pressure")),
        composition=_qd(g.get("composition")),
        spectroscopic_data=_qd(g.get("spectroscopic_data")),
        stage=None,
        gas_props = GasProps(ct.Solution("heat_transfer/config/flue_cantera.yaml"))

    )

    w = data["water_stream"]
    water = WaterStream(
        mass_flow_rate=_q(w.get("mass_flow_rate")),
        temperature=_q(w.get("temperature")),
        pressure=_q(w.get("pressure")),
        composition=_qd(w.get("composition")),
        stage=None,            # you can attach later
        water_props=WaterProps,      
    )

    return gas, water

# -----------------
# stage pieces
# -----------------

def _surface(n: Dict[str, Any]) -> Surface:
    return Surface(
        roughness=_q(n.get("roughness")),
        emissivity=_q(n.get("emissivity")),
        fouling_thickness=_q(n.get("fouling_thickness")),
        fouling_conductivity=_q(n.get("fouling_conductivity")),
    )

def _wall(n: Dict[str, Any]) -> Wall:
    s = n.get("surfaces") or {}
    inner = _surface(s.get("inner", {}))
    outer = _surface(s.get("outer", {}))
    return Wall(
        thickness=_q(n.get("thickness")),
        conductivity=_q(n.get("conductivity")),
        surfaces=Surfaces(inner=inner, outer=outer),
    )

def _shell(n: Dict[str, Any]) -> Shell:
    g = n.get("wall", {})
    geom = n.get("geometry") or {}
    geometry = ShellGeometry(
        flow_area=_q(geom.get("flow_area")),
        wetted_perimeter=_q(geom.get("wetted_perimeter")),
    )
    return Shell(geometry=geometry, wall=_wall(g))

def _tube_geom(n: Dict[str, Any]) -> TubeGeometry:
    return TubeGeometry(
        inner_diameter=_q(n.get("inner_diameter")),
        inner_length=_q(n.get("inner_length")),
    )

def _bank_layout(n: Dict[str, Any]) -> BankLayout:
    return BankLayout(
        tubes_number=_q(n.get("tubes_number")),
        shape=(n.get("shape") or ""),
        pitch=_q(n.get("pitch")),
    )

def _nozzles(n: Dict[str, Any]) -> Nozzles:
    # YAML provides inner/outer objects but no fields yet
    return Nozzles(inner=Nozzle(), outer=Nozzle())

# -----------------
# stages
# -----------------

def load_stages_from_yaml(yaml_path: str) -> Stages:
    with open(yaml_path, "r", encoding="utf-8") as f:
        data = yaml.safe_load(f) or {}
    stages = data.get("stages") or {}

    # HX_1 -> FirePass (single tube pass)
    hx1 = stages.get("HX_1", {}) or {}
    firepass = FirePass(
        geometry=_tube_geom(hx1.get("geometry", {}) or {}),
        wall=_wall(hx1.get("wall", {}) or {}),
        shell=_shell(hx1.get("shell", {}) or {}),
    )

    # HX_2 -> Reversal
    hx2 = stages.get("HX_2", {}) or {}
    reversal2 = Reversal(
        geometry=ReversalChamberGeometry(),                        # no params in YAML
        nozzles=_nozzles(hx2.get("nozzles", {}) or {}),
        wall=_wall(hx2.get("wall", {}) or {}),
        shell=_shell(hx2.get("shell", {}) or {}),
    )

    # HX_3 -> SmokePass (tube bank)
    hx3 = stages.get("HX_3", {}) or {}
    smokepass3 = SmokePass(
        geometry=_tube_geom(hx3.get("geometry", {}) or {}),
        layout=_bank_layout(hx3.get("layout", {}) or {}),
        wall=_wall(hx3.get("wall", {}) or {}),
        shell=_shell(hx3.get("shell", {}) or {}),
    )

    # HX_4 -> Reversal
    hx4 = stages.get("HX_4", {}) or {}
    reversal4 = Reversal(
        geometry=ReversalChamberGeometry(),
        nozzles=_nozzles(hx4.get("nozzles", {}) or {}),
        wall=_wall(hx4.get("wall", {}) or {}),
        shell=_shell(hx4.get("shell", {}) or {}),
    )

    # HX_5 -> SmokePass (tube bank)
    hx5 = stages.get("HX_5", {}) or {}
    smokepass5 = SmokePass(
        geometry=_tube_geom(hx5.get("geometry", {}) or {}),
        layout=_bank_layout(hx5.get("layout", {}) or {}),
        wall=_wall(hx5.get("wall", {}) or {}),
        shell=_shell(hx5.get("shell", {}) or {}),
    )

    # HX_6 -> Economiser
    hx6 = stages.get("HX_6", {}) or {}
    economiser = Economiser(
        hot_side=HotSide(area=_q((hx6.get("hot_side") or {}).get("area"))),
        cold_side=ColdSide(area=_q((hx6.get("cold_side") or {}).get("area"))),
    )

    return Stages(
        HX_1=firepass,
        HX_2=reversal2,
        HX_3=smokepass3,
        HX_4=reversal4,
        HX_5=smokepass5,
        HX_6=economiser,
    )
