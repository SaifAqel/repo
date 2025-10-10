from __future__ import annotations
from common.units import Q_
from heat_transfer.config.loader import ConfigLoader
from heat_transfer.functions.stages_chain import MinCounterflowChain

def run(stages_path: str, streams_path: str | None = None):

    cfg_stages = ConfigLoader.load_stages(stages_path)
    gas_in = ConfigLoader.load_gas_stream(streams_path)
    water_in = ConfigLoader.load_water_stream(streams_path)

    solver = MinCounterflowChain(stages=cfg_stages, gas_in=gas_in, water_in=water_in)
    result = solver.run()

    return result
