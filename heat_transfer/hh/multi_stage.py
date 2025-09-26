class MultiStageRunner:
    def __init__(self, stages):
        # stages: list of dicts, each: {"runner": PassRunner, "L": float, "dx": float}
        self.stages = stages

    def run(self, T0, P0):
        results = []
        Ti, Pi = T0, P0
        for s in self.stages:
            out = s["runner"].run(s["L"], s["dx"], Ti, Pi)
            results.append(out)
            Ti = out["T"][-1]
            Pi = out["P"][-1]
        return {"stages": results, "T_out": Ti, "P_out": Pi}
