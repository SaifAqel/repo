from __future__ import annotations
from common.units import Q_
from heat_transfer.config.loader import load_stages_from_yaml, load_streams_from_yaml
from heat_transfer.functions.stages_chain import MinCounterflowChain

def run(stages_path: str, streams_path: str | None = None):

    cfg_stages = load_stages_from_yaml(stages_path)
    gas_in, water_in = load_streams_from_yaml(streams_path)

    solver = MinCounterflowChain(stages=cfg_stages, gas_in=gas_in, water_in=water_in)
    result = solver.run()

    return result
