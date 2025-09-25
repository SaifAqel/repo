#!/usr/bin/env python3
from dataclasses import asdict, is_dataclass
import argparse
import json
import sys

from heat_transfer.config.schemas import load_config


def q_to_str(obj):
    """Pint Quantity -> 'value unit'; passthrough otherwise."""
    try:
        return f"{obj.magnitude:g} {obj.units}"
    except Exception:
        return obj


def to_jsonable(obj):
    """Dataclasses -> dict; dict/list -> recurse; Quantities -> strings."""
    if is_dataclass(obj):
        return {k: to_jsonable(v) for k, v in asdict(obj).items()}
    if isinstance(obj, dict):
        return {k: to_jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [to_jsonable(v) for v in obj]
    return q_to_str(obj)


def main() -> int:
    p = argparse.ArgumentParser(description="Load and pretty-print heat transfer config.")
    p.add_argument("--config", default="heat_transfer/config/settings.toml", help="Path to settings.toml")
    p.add_argument("--units", default="heat_transfer/config/units.toml", help="Path to units.toml")
    p.add_argument("--raw", action="store_true", help="Print Python repr instead of JSON")
    args = p.parse_args()

    cfg = load_config(args.config, args.units)

    if args.raw:
        print(cfg)
    else:
        print(json.dumps(to_jsonable(cfg), indent=2))

    # Targeted reads (kept explicit and readable)
    print(q_to_str(cfg.gas_side.mass_flow_rate.to("kg/s")))
    print(q_to_str(cfg.stages.drum.geometry.wall.thickness.to("m")))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
