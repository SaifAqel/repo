# main.py
import os, json, warnings
import numpy as np
np.seterr(all="raise")            # turn NumPy runtime warnings into errors
warnings.filterwarnings("error")  # escalate Python warnings to errors

# Use non-interactive backend for headless runs
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import tomllib

from heat_transfer1.services.plant_runner import MultiStageHEX
from heat_transfer1.services.result_assembler import assemble_summary
from heat_transfer1.io.profile_writer import export_profiles
from heat_transfer1.io.summary_writer import export_summary
from heat_transfer1.reports.plots import plot_profiles


class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.complexfloating):
            return {"real": float(obj.real), "imag": float(obj.imag)}
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super().default(obj)


def main():
    settings_path = "heat_transfer/config/settings.toml"

    # Load settings to get output_dir. Rest is handled by loaders.
    with open(settings_path, "rb") as f:
        settings_toml = tomllib.load(f) or {}

    output_dir = settings_toml.get("output_dir", "outputs")
    os.makedirs(output_dir, exist_ok=True)

    # Run plant
    runner = MultiStageHEX(settings_path)
    plant = runner.run()

    # Summarize
    summary = assemble_summary(plant)

    # Export per-stage profiles and plots
    for i, stage in enumerate(plant.stages, start=1):
        csv_path = os.path.join(output_dir, f"stage_{i}_profiles.csv")
        export_profiles(stage, csv_path)
        plot_profiles(stage)

    # Save all open figures
    figs = plt.get_fignums()
    for j, num in enumerate(figs, start=1):
        plt.figure(num)
        plt.savefig(os.path.join(output_dir, f"plot_{j}.png"), dpi=150, bbox_inches="tight")
    plt.close("all")

    # Write summary JSON
    summary_path = os.path.join(output_dir, "summary.json")
    export_summary(summary, summary_path)

    # Console outputs
    print(json.dumps(summary, indent=2, cls=NumpyEncoder))
    print(json.dumps(
        {
            "profiles": [f"stage_{i}_profiles.csv" for i in range(1, len(plant.stages) + 1)],
            "plots":    [f"plot_{j}.png" for j in range(1, len(figs) + 1)],
            "summary":  "summary.json",
        },
        indent=2,
        cls=NumpyEncoder,
    ))


if __name__ == "__main__":
    main()
