from pathlib import Path
from pprint import pprint
import numpy as np
from common.units import ureg, Q_
from heat_transfer.config.schemas import load_config
from heat_transfer.calc_ops.stage_with_calc import with_calc as stages_with_calc
from heat_transfer.calc_ops.stream_with_calc import GasStream, Water, GasStreamWithCalc, WaterWithCalc
from heat_transfer.calc_ops.stage_with_calc import ReversalWithCalc
from heat_transfer.solver.nozzle import Nozzle
from heat_transfer.fluid_props.GasProps import GasProps
from heat_transfer.fluid_props.WaterProps import WaterProps
from heat_transfer.solver.ODE import FireTubeGasODE
import cantera as ct
from heat_transfer.runner.ChainStages import ChainStages


# Changes to heat_transfer/runner/run.py
# Replace the existing run function with this version.
# Add import for the new class at the top:
# from heat_transfer.runner.chain_stages import ChainStages

def run(config_path: str, units_path: str, mech_yaml_path: str):
    cfg = load_config(config_path, units_path)
    cfg_with_calc = stages_with_calc(cfg)

    gas = cfg.gas_side
    gas_stream = GasStream(
        mass_flow_rate=gas.mass_flow_rate,
        temperature=gas.inlet_temperature,
        pressure=gas.inlet_pressure,
        composition=gas.composition,
        spectroscopic_data=gas.spectroscopic_data
    )

    gas_cantera = ct.Solution(mech_yaml_path)
    gas_props = GasProps(gas_cantera)

    water = cfg.water_side
    water_stream = Water(
        mass_flow_rate=water.mass_flow_rate,
        inlet_temperature=water.inlet_temperature,
        outlet_temperature=water.outlet_temperature,
        pressure=water.pressure,
        superheat=getattr(water, "superheat", None),  # Assuming optional
        mu_exp=water.mu_exp,
        C_sf=water.C_sf,
        n=water.n,
        composition=water.composition
    )

    water_with_calc = WaterWithCalc(
        Water=water_stream,
        WaterProps=WaterProps
    )

    # pprint(...)  # Optional: keep or remove

    # Now use the new ChainStages class
    chain = ChainStages(
        cfg_with_calc=cfg_with_calc,
        gas_props=gas_props,
        water_with_calc=water_with_calc,
        gas_stream=gas_stream  # Passed and updated in-place
    )
    z, T, p, Q_dot = chain.run_chain()

    return z, T, p, Q_dot