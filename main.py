# run.py
from heat_transfer.build_and_run import build_and_run

if __name__ == "__main__":
    config_path = "heat_transfer/config/settings.toml"
    units_path = "heat_transfer/config/units.toml"
    mech_yaml_path = "heat_transfer/fluid_props/flue_cantera.yaml"

    z, T, p, Q_dot = build_and_run(config_path, units_path, mech_yaml_path)

    print("z:", z)
    print("T:", T)
    print("p:", p)
    print("Q_dot:", Q_dot)
