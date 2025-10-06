from heat_transfer.runner.run import run

if __name__ == "__main__":
    config_path = "heat_transfer/config/settings.toml"
    units_path = "heat_transfer/config/units.toml"
    mech_yaml_path = "heat_transfer/fluid_props/flue_cantera.yaml"

    z, T, p, Q_dot = run(config_path, units_path, mech_yaml_path)

    print("z:", z)
    print("T:", T)
    print("p:", p)
    print("Q_dot:", Q_dot)
