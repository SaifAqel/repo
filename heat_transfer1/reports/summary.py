# reports/summary.py
from common.units import ureg, Q_

def stage_summary(sr):
    return {
        "inlet_outlet": {
            "T_out": sr.T_out.to("kelvin"),
            "P_out": sr.P_out.to("pascal"),
        },
        "duty": sr.Q.to("watt"),
        "dp": sr.dP.to("pascal"),
        "wall_T": {
            "min": sr.Tw_min.to("kelvin"),
            "max": sr.Tw_max.to("kelvin"),
        },
        "averages": {
            "h_g": sr.h_g_avg.to("watt/(meter**2*kelvin)"),
            "epsilon_g": sr.eps_g_avg.to("dimensionless"),
        },
    }

def plant_summary(pr):
    return {
        "totals": {
            "Q": pr.Q_total.to("watt"),
            "dP": pr.dP_total.to("pascal"),
            "T_final": pr.T_final.to("kelvin"),
            "P_final": pr.P_final.to("pascal"),
        },
        "stages": [stage_summary(s) for s in pr.stages],
    }
