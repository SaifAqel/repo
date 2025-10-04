# heat_transfer/pipeline.py
import numpy as np
import cantera as ct

from heat_transfer.config.schemas import load_config, Drum, Pass, Reversal
from heat_transfer.config.with_calc import with_calc

from heat_transfer.fluid_props.GasProps import HotFlueGas
from heat_transfer.fluid_props.WaterProps import WaterProps

from heat_transfer.thermal.thermal_resistance import ThermalResistance

from heat_transfer.pressure_drop.pressure_drop import TubePressureDrop
from heat_transfer.thermal.heat import HeatSystem

from heat_transfer.solver.ODE import FireTubeGasODE
from heat_transfer.solver.train import MultiPassTrain
from heat_transfer.solver.nozzle import Nozzle


def build_and_run(config_path, units_path, mech_yaml_path):
    # Load configuration
    cfg = load_config(config_path, units_path)
    cfg_calc = with_calc(cfg)

    # Create initial gas stream from composition
    gas = ct.Solution(mech_yaml_path)
    gas.X = {k: float(v) for k, v in cfg.gas_side.composition.items()}

    # Instantiate properties
    gas_props = HotFlueGas(gas)
    water_props = WaterProps()
    m_dot = cfg.gas_side.mass_flow_rate
    T_in = cfg.gas_side.inlet_temperature
    p_in = cfg.gas_side.inlet_pressure

    # Per-pass builder
    friction = TubePressureDrop()

    # Keep track of results
    results = []

    def build_pass(pass_with_calc, gas_in):
        g = pass_with_calc.geometry  # expects: D_i, A, L, rel_rough, segment_z()
        rr = ThermalResistance.total_resistance

        # Heat system setup using proper streams
        heat = HeatSystem(
            PassWithCalc=g,
            total_thermal_resistance=rr,
            GasStream=gas_in,
            WaterStream=water_props,
            HotFlueGas=gas_props
        )

        # Solve ODE for gas
        ode = FireTubeGasODE(
            pass_geom=g,
            friction=friction,
            m_dot=m_dot,
        )
        z0 = 0.0
        z1 = g.inner_length
        y0 = [T_in, p_in]
        sol = ode.solve((z0, z1), y0, dense_output=True)

        # Store results along the pass
        z = g.segment_z()
        T = sol.sol(z)[0]
        p = sol.sol(z)[1]

        # Update gas for next pass (assuming perfect mixing)
        gas_next = gas_in.clone()
        gas_next.TP = T[-1], p[-1]  # update temperature and pressure

        return {"type": "pass", "z": z, "T": T, "p": p, "ode": ode, "heat": heat, "gas_out": gas_next}

    # Convert stages dataclass to a list of its fields
    stage_list = [
        cfg_calc.stages.drum,
        cfg_calc.stages.pass1,
        cfg_calc.stages.reversal1,
        cfg_calc.stages.pass2,
        cfg_calc.stages.reversal2,
        cfg_calc.stages.pass3
    ]

    gas_current = gas
    results = []

    for stage in stage_list:
        # Only handle Pass stages for now
        if isinstance(stage, Pass):
            res = build_pass(stage, gas_current)
            results.append(res)
            gas_current = res["gas_out"]  # propagate gas
        else:
            # For now, skip non-pass stages
            continue

    return results
