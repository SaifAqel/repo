import numpy as np

class MultiPassTrain:
    
    def __init__(self, segments):
        self.segments = segments

    def run(self, T_in, p_in):
        z_all = []
        T_all = []
        p_all = []

        z_offset = 0.0
        T_curr = T_in
        p_curr = p_in

        for seg in self.segments:
            if seg["type"] == "pass":
                z_eval = seg["z_eval"]
                z_span = (0.0, float(z_eval[-1]))
                sol = seg["solver"].solve(z_span, [T_curr, p_curr], t_eval=z_eval)
                z_seg = z_offset + sol.t
                T_seg = sol.y[0]
                p_seg = sol.y[1]

                z_all.append(z_seg)
                T_all.append(T_seg)
                p_all.append(p_seg)

                z_offset = z_seg[-1]
                T_curr = T_seg[-1]
                p_curr = p_seg[-1]

            elif seg["type"] == "nozzle":
                T_curr, p_curr = seg["nozzle"].apply(T_curr, p_curr)
                # no z added for nozzle

        if z_all:
            z_profile = np.concatenate(z_all)
            T_profile = np.concatenate(T_all)
            p_profile = np.concatenate(p_all)
        else:
            z_profile = np.array([0.0])
            T_profile = np.array([T_curr])
            p_profile = np.array([p_curr])

        return {"z": z_profile, "T": T_profile, "p": p_profile}
