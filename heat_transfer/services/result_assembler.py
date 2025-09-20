# services/result_assembler.py
from ..core.types import PlantResult
def assemble_summary(plant: PlantResult):
    return {
        "Q_total": plant.Q_total,
        "dP_total": plant.dP_total,
        "T_final": plant.T_final,
        "P_final": plant.P_final,
        "stages":[{"Q":s.Q,"dP":s.dP,"T_out":s.T_out,"P_out":s.P_out,"Tw_min":s.Tw_min,"Tw_max":s.Tw_max,"h_g_avg":s.h_g_avg,"eps_g_avg":s.eps_g_avg} for s in plant.stages]
    }
