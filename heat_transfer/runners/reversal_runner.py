# ---- ReversalRunner: identical API to PassRunner ----
class ReversalRunner:
    def __init__(self, geom: PassGeom, heat_model, flow="counterflow"):
        self.geom = geom
        self.flow = flow
        self.seg = SegmentRunner(heat_model)

    def run(self, hot_in: Stream, cold_in: Stream, N: int):
        cells = slice_geom(self.geom, N)
        Th = hot_in
        TcN = cold_in
        Q_tot = 0.0
        for cell in cells:
            Th, TcN, Q = self.seg.step_counterflow(Th, TcN, cell)
            Q_tot += Q
        return Th, TcN, Q_tot
