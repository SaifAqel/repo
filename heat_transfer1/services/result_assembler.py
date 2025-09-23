# services/result_assembler.py
from ..core.types import PlantResult

def assemble_summary(plant: PlantResult):
    return {
        "Q_total": plant.Q_total.to("watt"),
        "dP_total": plant.dP_total.to("pascal"),
        "T_final": plant.T_final.to("kelvin"),
        "P_final": plant.P_final.to("pascal"),
        "stages": [
            {
                "Q": s.Q.to("watt"),
                "dP": s.dP.to("pascal"),
                "T_out": s.T_out.to("kelvin"),
                "P_out": s.P_out.to("pascal"),
                "Tw_min": s.Tw_min.to("kelvin"),
                "Tw_max": s.Tw_max.to("kelvin"),
                "h_g_avg": s.h_g_avg.to("watt/(meter**2*kelvin)"),
                "eps_g_avg": s.eps_g_avg.to("dimensionless"),
            }
            for s in plant.stages
        ],
    }
