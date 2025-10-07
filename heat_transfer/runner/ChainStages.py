# Add this new file: heat_transfer/runner/chain_stages.py

from typing import List, Tuple
import numpy as np
from common.units import ureg, Q_
from heat_transfer.calc_ops.stage_with_calc import PassWithCalc, ReversalWithCalc, ConfigWithCalc
from heat_transfer.calc_ops.stream_with_calc import GasStream, GasStreamWithCalc, WaterWithCalc
from heat_transfer.fluid_props.GasProps import GasProps
from heat_transfer.solver.ODE import FireTubeGasODE
from heat_transfer.solver.nozzle import Nozzle
from dataclasses import dataclass

# Placeholder classes for economizer and superheater.
# These are defined similarly to PassWithCalc for now, assuming multi-tube geometries.
# They can be refined later with specific attributes (e.g., finned tubes for economizer).
# For now, they inherit from PassWithCalc, but you can override properties as needed.

@dataclass
class EconomizerWithCalc(PassWithCalc):
    """
    Placeholder for economizer stage. Assumes similar to a tube pass for gas flow.
    Refine later: e.g., add fin details, external flow correlations if needed.
    Water side might be single-phase (subcooled), so adjust heat transfer accordingly.
    """
    # Add specific fields if needed, e.g., fin_efficiency: Q_
    pass  # Well-defined later

@dataclass
class SuperheaterWithCalc(PassWithCalc):
    """
    Placeholder for superheater stage. Assumes similar to a tube pass for gas flow.
    Refine later: e.g., specific radiation view factors, superheated steam correlations.
    Water side is vapor phase, so use vapor convection instead of boiling.
    """
    # Add specific fields if needed, e.g., coil_arrangement: str
    pass  # Well-defined later

class ChainStages:
    """
    Standalone class to chain all stages, including placeholders for economizer and superheater.
    Handles ODE solving per stage, nozzle pressure drops for reversals, and collects results.
    """
    def __init__(
        self,
        cfg_with_calc: ConfigWithCalc,
        gas_props: GasProps,
        water_with_calc: WaterWithCalc,
        gas_stream: GasStream  # Initial gas stream; updated in-place
    ):
        self.cfg_with_calc = cfg_with_calc
        self.gas_props = gas_props
        self.water_with_calc = water_with_calc
        self.gas_stream = gas_stream  # Mutable; updated across stages
        
        # Define stage sequence. Insert placeholders where appropriate.
        # Assuming: superheater after pass1 (hot zone), economizer after pass3 (cool zone).
        # Adjust sequence based on actual boiler design.
        self.stages = [
            cfg_with_calc.stages.pass1,
            # Superheater placeholder here (after first pass, hottest gas)
            SuperheaterWithCalc(
                geometry=cfg_with_calc.stages.pass1.geometry,  # Placeholder: copy pass1 geom; redefine later
                surfaces=cfg_with_calc.stages.pass1.surfaces
            ),
            cfg_with_calc.stages.reversal1,
            cfg_with_calc.stages.pass2,
            cfg_with_calc.stages.reversal2,
            cfg_with_calc.stages.pass3,
            # Economizer placeholder here (after last pass, cooler gas)
            EconomizerWithCalc(
                geometry=cfg_with_calc.stages.pass3.geometry,  # Placeholder: copy pass3 geom; redefine later
                surfaces=cfg_with_calc.stages.pass3.surfaces
            )
        ]

    def run_chain(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Run the chained simulation.
        Returns: z_total, T_total, p_total, Q_dot_total (all numpy arrays)
        """
        z_total = np.array([])
        T_total = np.array([])
        p_total = np.array([])
        Q_dot_total = np.array([])
        current_z = 0.0

        for stage in self.stages:
            # Handle reversal inlet nozzle if applicable
            if isinstance(stage, ReversalWithCalc):
                # Assume K=0.5 for contraction; refine later
                nozzle_in = Nozzle(
                    geom=stage.nozzles.inlet,
                    gas_stream=self.gas_stream,
                    gas_props=self.gas_props,
                    K=Q_(0.5, ureg.dimensionless)
                )
                T_new, p_new = nozzle_in.apply(
                    self.gas_stream.temperature,
                    self.gas_stream.pressure,

                )
                self.gas_stream.temperature = Q_(T_new)
                self.gas_stream.pressure = Q_(p_new)

            # Create stage-specific gas calc
            gas_calc = GasStreamWithCalc(
                gas_props=self.gas_props,
                gas_stream=self.gas_stream,
                geometry=stage
            )

            # Create ODE instance (note: for economizer/superheater, may need custom HeatSystem later)
            ode = FireTubeGasODE(
                pass_with_calc=stage,  # Works for placeholders since they inherit
                gas_stream_with_calc=gas_calc,
                water_with_calc=self.water_with_calc
            )

            # Initial conditions
            y0 = [self.gas_stream.temperature.magnitude, self.gas_stream.pressure.magnitude]
            z_span = (0, stage.geometry.inner_length.magnitude)  # Use path_length for consistency

            # Solve ODE
            res = ode.solve(z_span, y0)
            if not res.success:
                raise ValueError(f"ODE failed for stage {type(stage).__name__}: {res.message}")

            # Collect results
            z_stage = res.t + current_z
            T_stage = res.y[0]
            p_stage = res.y[1]

            # Compute Q_dot
            Q_dot_stage = []
            for ti in T_stage:
                ode.gas.gas_stream.temperature = Q_(ti, 'kelvin')
                Q_dot_stage.append(ode.heat_system.q_.magnitude)
            Q_dot_stage = np.array(Q_dot_stage)

            # Append to totals
            z_total = np.append(z_total, z_stage)
            T_total = np.append(T_total, T_stage)
            p_total = np.append(p_total, p_stage)
            Q_dot_total = np.append(Q_dot_total, Q_dot_stage)

            # Update current position and inlet for next stage
            current_z += z_span[1]
            self.gas_stream.temperature = Q_(T_stage[-1], 'kelvin')
            self.gas_stream.pressure = Q_(p_stage[-1], 'pascal')

            # Handle reversal outlet nozzle if applicable
            if isinstance(stage, ReversalWithCalc):
                # Assume K=1.0 for expansion; refine later
                nozzle_out = Nozzle(
                    geom=stage.nozzles.outlet,
                    gas_stream=self.gas_stream,
                    gas_props=self.gas_props,
                    K = Q_(0.5, ureg.dimensionless)
                )
                T_new, p_new = nozzle_out.apply(
                    self.gas_stream.temperature,
                    self.gas_stream.pressure
                )
                self.gas_stream.temperature = Q_(T_new, 'kelvin')
                self.gas_stream.pressure = Q_(p_new, 'pascal')

        return z_total, T_total, p_total, Q_dot_total