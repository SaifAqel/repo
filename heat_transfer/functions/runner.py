from __future__ import annotations
from common.units import Q_
from heat_transfer.config.loader import ConfigLoader
from heat_transfer.functions.stages_chain import six_stage_counterflow

def run(stages_path: str, streams_path: str):

    stages = ConfigLoader.load_stages(stages_path)
    gas_in = ConfigLoader.load_gas_stream(streams_path)
    water_in = ConfigLoader.load_water_stream(streams_path)

    solver = six_stage_counterflow(stages=stages)
    result = solver.run(gas=gas_in, water=water_in)

    return result
