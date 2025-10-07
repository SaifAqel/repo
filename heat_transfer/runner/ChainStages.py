# Add this new file: heat_transfer/runner/chain_stages.py

from typing import List, Tuple
import numpy as np
from common.units import ureg, Q_
from heat_transfer.calc_ops.stage_with_calc import PassWithCalc, ReversalWithCalc, ConfigWithCalc
from heat_transfer.calc_ops.stream_with_calc import GasStream, GasStreamWithCalc, Water, WaterWithCalc
from heat_transfer.fluid_props.GasProps import GasProps
from heat_transfer.fluid_props.WaterProps import WaterProps
from heat_transfer.solver.ODE import FireTubeGasODE
from heat_transfer.solver.nozzle import Nozzle
from dataclasses import dataclass
from heat_transfer.solver.water_balance import WaterStateConverter, WaterEnergyRHS

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
        water_stream: Water,
        gas_stream: GasStream,  # Initial gas stream; updated in-place
        water_props: WaterProps
    ):
        self.cfg = cfg_with_calc
        self.gas_props = gas_props
        self.water_stream = water_stream
        self.gas_stream = gas_stream  # Mutable; updated across stages
        self.water_props = water_props
        
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

            water_calc = WaterWithCalc(
                water_stream = self.water_stream,
                water_props = self.water_props
            )

            # prepare water helpers for this stage
            converter = WaterStateConverter(water_calc)
            water_rhs = WaterEnergyRHS(water_calc, stage)

            # Create ODE instance
            ode = FireTubeGasODE(
                pass_with_calc=stage,
                gas_stream_with_calc=gas_calc,
                water_with_calc=water_calc,
                water_rhs=water_rhs
            )

            # Initial conditions: gas T, gas p, system-level water enthalpy
            h0_sys = water_calc.enthalpy
            y0 = [self.gas_stream.temperature.magnitude, self.gas_stream.pressure.magnitude, h0_sys.to('J/kg').magnitude]
            z_span = (0, stage.geometry.inner_length.magnitude)

            # Solve ODE
            res = ode.solve(z_span, y0)
            if not res.success:
                raise ValueError(f"ODE failed for stage {type(stage).__name__}: {res.message}")

            # Collect results
            z_stage = res.t + current_z
            T_stage = res.y[0]
            p_stage = res.y[1]
            h_stage = res.y[2]

            # Compute Q_dot per axial point using paired gas T and water h
            Q_dot_stage = []
            for ti, hi in zip(T_stage, h_stage):
                ode.gas.gas_stream.temperature = Q_(ti, 'kelvin')
                # set local water temperature for HeatSystem only (transient)
                T_w_local = converter.T_from_h(Q_(hi, 'J/kg'), water_calc.water_stream.pressure)
                water_calc.water_stream.temperature = T_w_local
                Q_dot_stage.append(ode.heat_system.q_.magnitude)
            Q_dot_stage = np.array(Q_dot_stage)

            # Append to totals
            z_total = np.append(z_total, z_stage)
            T_total = np.append(T_total, T_stage)
            p_total = np.append(p_total, p_stage)
            Q_dot_total = np.append(Q_dot_total, Q_dot_stage)

            # Update current position and inlet for next stage (gas)
            current_z += z_span[1]
            self.gas_stream.temperature = Q_(T_stage[-1], 'kelvin')
            self.gas_stream.pressure = Q_(p_stage[-1], 'pascal')

            # update system-level water enthalpy for next stage
            # Get final enthalpy from ODE solution
            h_end = Q_(h_stage[-1], 'J/kg')

            # Convert enthalpy to temperature using your converter
            T_end = converter.T_from_h(h_end, water_calc.water_stream.pressure)
            water_calc.water_stream.temperature = T_end



        final_h = water_calc.enthalpy
        final_T = WaterStateConverter(water_calc).T_from_h(final_h, water_calc.water_stream.pressure)
        water_calc.water_stream.temperature = final_T


        return z_total, T_total, p_total, Q_dot_total