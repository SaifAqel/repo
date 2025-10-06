from pathlib import Path
from pprint import pprint

from heat_transfer.config.schemas import load_config
from heat_transfer.calc_ops.stage_with_calc import with_calc as stages_with_calc
from heat_transfer.calc_ops.stream_with_calc import GasStream, Water, GasStreamWithCalc, WaterWithCalc
from heat_transfer.fluid_props.GasProps import GasProps
from heat_transfer.fluid_props.WaterProps import WaterProps
from heat_transfer.solver.ODE import FireTubeGasODE
import cantera as ct


def run(config_path: str, units_path: str, mech_yaml_path: str):
    cfg = load_config(config_path, units_path)
    cfg_with_calc = stages_with_calc(cfg)

    gas = cfg.gas_side
    gas_stream = GasStream(
        mass_flow_rate=gas.mass_flow_rate,
        temperature=gas.inlet_temperature,
        pressure=gas.inlet_pressure,
        composition=gas.composition,
    )

    gas_cantera = ct.Solution('gri30.yaml')

    gas_stream_with_calc = GasStreamWithCalc(
        gas_props=GasProps(gas_cantera),      # your GasProps instance
        gas_stream=gas_stream,
        geometry=cfg_with_calc.stages.pass1  # or the appropriate geometry object
    )

    water = cfg.water_side
    water_stream = Water(
        mass_flow_rate = getattr(water, "mass_flow_rate", None),
        inlet_temperature = getattr(water, "inlet_temperature", None),
        outlet_temperature = getattr(water, "outlet_temperature", None),
        pressure = getattr(water, "pressure", getattr(water, "outlet_pressure", None)),
        superheat = getattr(water, "superheat", None),
        mu_exp = getattr(water, "mu_exp", None),
        C_sf = getattr(water, "C_sf", None),
        n = getattr(water, "n", None),
        composition = getattr(water, "composition", {}),
    )

    water_with_calc = WaterWithCalc(
        Water = water_stream,
        WaterProps = WaterProps
    )

    pprint({"loaded_config": cfg})
    pprint({"stages_with_calc_pass1": cfg_with_calc.stages.pass1})
    pprint({"example_gas_stream": gas_stream_with_calc})
    pprint({"example_water_stream": water_stream})

    ode = FireTubeGasODE(cfg_with_calc.stages.pass1, gas_stream_with_calc, water_with_calc)
    res = ode.rhs()

    # must return z, T, p, Q_dot for caller
    return res.get("z"), res.get("T"), res.get("p"), res.get("Q_dot")