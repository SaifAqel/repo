from heat_transfer.config.loader import ConfigLoader
from heat_transfer.calc_ops.stream_with_calc import GasStream, Water
from heat_transfer.fluid_props.GasProps import GasProps
from heat_transfer.fluid_props.WaterProps import WaterProps
from heat_transfer.runner.ChainStages import ChainStages
from heat_transfer.calc_ops.stage_with_calc import with_calc
import cantera as ct


# Changes to heat_transfer/runner/run.py
# Replace the existing run function with this version.
# Add import for the new class at the top:
# from heat_transfer.runner.chain_stages import ChainStages

def run(config_path: str, mech_yaml_path: str):
    cfg = ConfigLoader.load_from_path(config_path)
    cfg_with_calc = with_calc(cfg)

    gas_props = GasProps(ct.Solution(mech_yaml_path))
    water_props = WaterProps()

    gas = cfg.gas_inlet
    gas_stream = GasStream(
        mass_flow_rate=gas.mass_flow_rate,
        temperature=gas.temperature,
        pressure=gas.pressure,
        composition=gas.composition,
        spectroscopic_data=gas.spectroscopic_data,
        z=gas.z
    )
    
    water = cfg.water_inlet
    water_stream = Water(
        mass_flow_rate=water.mass_flow_rate,
        temperature=water.temperature,
        pressure=water.pressure,
        composition=water.composition,
        z=water.z
    )

    chain = ChainStages(
        cfg_with_calc=cfg_with_calc,
        gas_props=gas_props,
        water_stream=water_stream,
        gas_stream=gas_stream,
        water_props=water_props
    )
    z, T, p, Q_dot = chain.run_chain()

    return z, T, p, Q_dot