# main.py
import json, os, yaml
from heat_transfer.services.plant_runner import MultiStageHEX
from heat_transfer.services.result_assembler import assemble_summary
from heat_transfer.io.profile_writer import export_profiles
from heat_transfer.io.summary_writer import export_summary
from heat_transfer.reports.plots import plot_profiles
import matplotlib.pyplot as plt
import numpy as np, warnings

np.seterr(all='raise')           # turn RuntimeWarnings into errors
warnings.filterwarnings('error') # catch ComplexWarning as error

class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            return float(obj)
        if isinstance(obj, (np.complexfloating,)):
            return {"real": obj.real, "imag": obj.imag}
        if isinstance(obj, (np.ndarray,)):
            return obj.tolist()
        return super().default(obj)

def build_settings_yaml(N_cells_default, tolerances_rel, tolerances_abs, max_iters, backend, wall_material, R_fg_default, R_fw_default):
    return {
        "global_numerics": {
            "tolerances": {"rel": tolerances_rel, "abs": tolerances_abs},
            "max_iters": max_iters,
            "N_cells_default": N_cells_default
        },
        "property_backend": {"backend": "CoolProp", "coolprop": {"backend": backend}},
        "defaults": {"wall_material": wall_material, "fouling": {"R_fg": R_fg_default, "R_fw": R_fw_default}}
    }

def build_case_yaml(T_in, P_in, m_dot, composition, drum_pressure, drum_geometry, drum_level, drum_inventory, stages):
    return {
        "inlet_gas": {"T_in": T_in, "P_in": P_in, "m_dot": m_dot, "composition": composition},
        "drums": [{"id":"drum1","pressure": drum_pressure,"level": drum_level,"geometry": drum_geometry,"initial_water_inventory": drum_inventory}],
        "stages": stages
    }

def main():
    settings_path = "heat_transfer/config/settings.yaml"
    case_path = "heat_transfer/config/case_stages.yaml"

    # Load YAMLs
    with open(settings_path, "r") as f:
        settings_yaml = yaml.safe_load(f) or {}

    with open(case_path, "r") as f:
        case_yaml = yaml.safe_load(f) or {}

    # Output directory: optional in case YAML, else default
    output_dir = case_yaml.get("output_dir", "outputs")
    os.makedirs(output_dir, exist_ok=True)

    # Run model using the YAML files directly
    runner = MultiStageHEX(settings_path, case_path)
    plant = runner.run()

    # Summarize
    summary = assemble_summary(plant)

    # Export per-stage profiles and plots
    for i, sr in enumerate(plant.stages, start=1):
        profile_path = os.path.join(output_dir, f"stage_{i}_profiles.csv")
        export_profiles(sr, profile_path)
        plot_profiles(sr)

    # Save all open figures
    figs = plt.get_fignums()
    for j, num in enumerate(figs, start=1):
        plt.figure(num)
        plt.savefig(os.path.join(output_dir, f"plot_{j}.png"), dpi=150, bbox_inches="tight")
    plt.close("all")

    # Write summary
    summary_path = os.path.join(output_dir, "summary.json")
    export_summary(summary, summary_path)

    # Console outputs (paths + summary)
    print(json.dumps(summary, indent=2, cls=NumpyEncoder))
    print(json.dumps({
        "profiles": [f"stage_{i}_profiles.csv" for i in range(1, len(plant.stages) + 1)],
        "plots":    [f"plot_{j}.png" for j in range(1, len(figs) + 1)],
        "summary":  "summary.json"
    }, indent=2, cls=NumpyEncoder))



if __name__ == "__main__":
    main()
