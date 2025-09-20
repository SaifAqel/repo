# reports/summary.py
def stage_summary(sr):
    return {
        "inlet_outlet": {"T_out": sr.T_out, "P_out": sr.P_out},
        "duty": sr.Q,
        "dp": sr.dP,
        "wall_T": {"min": sr.Tw_min, "max": sr.Tw_max},
        "averages": {"h_g": sr.h_g_avg, "epsilon_g": sr.eps_g_avg}
    }
def plant_summary(pr):
    return {
        "totals": {"Q": pr.Q_total, "dP": pr.dP_total, "T_final": pr.T_final, "P_final": pr.P_final},
        "stages": [stage_summary(s) for s in pr.stages]
    }
