# loader.py
from __future__ import annotations
import yaml
from typing import Dict, Tuple, Any
from common.units import ureg, Q_
from heat_transfer.config.models import (
    Surface, Surfaces, Wall, Nozzle, Nozzles,
    TubeGeometry, ReversalChamberGeometry, BankLayout,
    ShellGeometry, Shell, TubeBank, ReversalChamber,
    FirePass, SmokePass, Reversal, HotSide, ColdSide, Economiser, Stages,
    WaterStream, GasStream
)

def _q(node: Dict) -> Q_:
    if not node:  # allow None
        return None
    v = node.get("value")
    u = node.get("unit", "")
    if v is None:
        return None
    return Q_(v, u)

def _qd(mapping: Dict[str, Dict]) -> Dict[str, Q_]:
    return {k: _q(v) for k, v in (mapping or {}).items()}

    


def load_streams_from_yaml(yaml_path: str) -> Tuple[GasStream, WaterStream]:

        with open(yaml_path, "r", encoding="utf-8") as f:
            data = yaml.safe_load(f)

        g = data["gas_stream"]
        gas = GasStream(
            mass_flow_rate=_q(g["mass_flow_rate"]),
            temperature=_q(g["temperature"]),
            pressure=_q(g["pressure"]),
            composition=_qd(g["composition"]),
            spectroscopic_data=_qd(g["spectroscopic_data"]),
            stage=None,
        )

        w = data["water_stream"]
        water = WaterStream(
            mass_flow_rate=_q(w["mass_flow_rate"]),
            temperature=_q(w["temperature"]),
            pressure=_q(w["pressure"]),
            composition=_qd(w["composition"]),
            stage=None,
        )

        return gas, water
    
def _surface(n: Dict[str, Any]) -> Surface:
    return Surface(
        roughness=_q(n.get("roughness")),
        emissivity=_q(n.get("emissivity")),
        fouling_thickness=_q(n.get("fouling_thickness")),
        fouling_conductivity=_q(n.get("fouling_conductivity")),
    )

def _wall(n: Dict[str, Any]) -> Wall:
    inner = _surface((n.get("surfaces") or {}).get("inner", {}))
    outer = _surface((n.get("surfaces") or {}).get("outer", {}))
    return Wall(
        thickness=_q(n.get("thickness")),
        conductivity=_q(n.get("conductivity")),
        surfaces=Surfaces(inner=inner, outer=outer),
    )

def _shell(n: Dict[str, Any]) -> Shell:
    geometry = ShellGeometry(
        flow_area=_q((n.get("geometry") or {}).get("flow_area")),
        wetted_perimeter=_q((n.get("geometry") or {}).get("wetted_perimeter")),
    )
    return Shell(geometry=geometry, wall=_wall(n.get("wall", {})))

def _tube_geom(n: Dict[str, Any]) -> TubeGeometry:
    return TubeGeometry(
        inner_diameter=_q(n.get("inner_diameter")),
        inner_length=_q(n.get("inner_length")),
    )

def _bank_layout(n: Dict[str, Any]) -> BankLayout:
    return BankLayout(
        tubes_number=_q(n.get("tubes_number")),
        shape=n.get("shape", "") or "",
        pitch=_q(n.get("pitch")),
    )

def _nozzles(n: Dict[str, Any]) -> Nozzles:
    return Nozzles(inner=Nozzle(), outer=Nozzle())

# ---------- main loader ----------
def load_stages_from_yaml(yaml_path: str) -> Stages:
    with open(yaml_path, "r", encoding="utf-8") as f:
        data = yaml.safe_load(f) or {}
    stages = (data.get("stages") or {})

    # HX_1 -> FirePass
    hx1 = stages.get("HX_1", {})
    firepass = FirePass(
        geom=_tube_geom(hx1.get("geometry", {})),
        wall=_wall(hx1.get("wall", {})),
        shell=_shell(hx1.get("shell", {})),
    )

    # HX_2 -> Reversal
    hx2 = stages.get("HX_2", {})
    reversal2 = Reversal(
        geometry=ReversalChamber(
            geometry=ReversalChamberGeometry(),
            nozzles=_nozzles(hx2.get("nozzles", {})),
        ),
        wall=_wall(hx2.get("wall", {})),
        shell=_shell(hx2.get("shell", {})),
    )

    # HX_3 -> SmokePass
    hx3 = stages.get("HX_3", {})
    smokepass3 = SmokePass(
        geom=TubeBank(
            geom=_tube_geom(hx3.get("geometry", {})),  # YAML uses "geometry"
            layout=_bank_layout(hx3.get("layout", {})),
        ),
        wall=_wall(hx3.get("wall", {})),
        shell=_shell(hx3.get("shell", {})),
    )

    # HX_4 -> Reversal
    hx4 = stages.get("HX_4", {})
    reversal4 = Reversal(
        geom=ReversalChamber(
            geometry=ReversalChamberGeometry(),
            nozzles=_nozzles(hx4.get("nozzles", {})),
        ),
        wall=_wall(hx4.get("wall", {})),
        shell=_shell(hx4.get("shell", {})),
    )

    # HX_5 -> SmokePass (note: YAML uses "geom" key for tube geometry)
    hx5 = stages.get("HX_5", {})
    smokepass5 = SmokePass(
        geom=TubeBank(
            geom=_tube_geom(hx5.get("geom", {})),  # YAML key is "geom"
            layout=_bank_layout(hx5.get("layout", {})),
        ),
        wall=_wall(hx5.get("wall", {})),
        shell=_shell(hx5.get("shell", {})),
    )

    # HX_6 -> Economiser
    hx6 = stages.get("HX_6", {})
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

# Usage:
# stages = load_stages_from_yaml("stages.yaml")