from __future__ import annotations
from common.units import Q_
from heat_transfer.config.loader import ConfigLoader
from heat_transfer.functions.stages_chain import MinCounterflowChain

def run(stages_path: str, streams_path: str):

    stages = ConfigLoader.load_stages(stages_path)
    gas = ConfigLoader.load_gas_stream(streams_path)
    water = ConfigLoader.load_water_stream(streams_path)

    solver = MinCounterflowChain(stages=stages, gas=gas, water=water)
    result = solver.run()

    return result
