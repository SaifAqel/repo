import copy
from typing import List, Tuple
from heat_transfer.config.models import Stages, GasStream, WaterStream
from heat_transfer.functions.stage_solver import StageSolver

class six_stage_counterflow:
    def __init__(self, stages: Stages):
        self.stages = stages

    def run(self, gas: GasStream, water: WaterStream) -> Tuple[List[GasStream], List[WaterStream]]:
        gas_hist: List[GasStream] = []
        water_hist: List[WaterStream] = []

        for stage in self.stages:
            solver = StageSolver(stage=stage, gas=gas, water=water)
            g_list, w_list = solver.solve()  # uses your earlier solve() that returns lists of instances
            gas_hist.extend(copy.deepcopy(g_list))
            water_hist.extend(copy.deepcopy(w_list))
            # gas and water are already updated in-place to stage outlet; they feed the next stage

        return gas_hist, water_hist
 