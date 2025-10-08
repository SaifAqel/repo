from heat_transfer.runner.run import run
from post_processing.post_process import run_post_processing

if __name__ == "__main__":
    config_path = "heat_transfer/config/config.yaml"
    mech_yaml_path = "heat_transfer/fluid_props/flue_cantera.yaml"

    gas_profile, water_profile = run(config_path, mech_yaml_path)

    run_post_processing(gas_profile, water_profile)

    print("Gas points:")
    for gp in gas_profile.points[:]:
        print(f"z: {gp.z}, T: {gp.temperature}, p: {gp.pressure}, Q_dot: {gp.heat_transfer_rate}, velocity: {gp.velocity}, Re: {gp.reynolds_number}, density: {gp.density}")

    print("Water points:")
    for wp in water_profile.points[:]:
        print(f"z: {wp.z}, T: {wp.temperature}, h: {wp.enthalpy}, Q_dot: {wp.heat_transfer_rate}, quality: {wp.quality}, phase: {wp.phase}, Tsat: {wp.saturation_temperature}")