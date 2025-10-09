from heat_transfer.dump.loader import ConfigLoader
from heat_transfer.functions.GasProps import GasProps
from heat_transfer.functions.WaterProps import WaterProps
from heat_transfer.dump.ChainStages import ChainStages
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

    gas_inlet = cfg.gas_inlet
    water_inlet = cfg.water_inlet

    chain = ChainStages(
        cfg_with_calc=cfg_with_calc,
        gas_props=gas_props,
        water_stream=water_inlet,
        gas_stream=gas_inlet,
        water_props=water_props
    )
    gas_profile, water_profile = chain.run_chain()

    return gas_profile, water_profile