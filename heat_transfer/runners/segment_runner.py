from heat_transfer.thermal.thermal_resistance import *
from heat_transfer.stage_with_calc import *

# ---- 1) SegmentRunner: one cell update using your heat-rate class ----
class SegmentRunner:
    def __init__(self, PassWithCalc):
        self.PassWithCalc = PassWithCalc
        
    def run_segment(self):
        s_L = self.PassWithCalc.segment_heat_transfer_area

        










