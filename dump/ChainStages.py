# Add this new file: heat_transfer/runner/chain_stages.py

from typing import List, Tuple
import numpy as np
from common.units import ureg, Q_
from heat_transfer.config.models import WaterStreamProfile, GasStreamProfile
from heat_transfer.calc_ops.stage_with_calc import PassWithCalc, ReversalWithCalc, ConfigWithCalc
from heat_transfer.calc_ops.stream_with_calc import GasStream, GasStreamWithCalc, Water, WaterWithCalc
from heat_transfer.functions.GasProps import GasProps
from heat_transfer.functions.WaterProps import WaterProps
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
        self.water_inlet = water_stream  # Keep inlet immutable
        self.gas_inlet = gas_stream  # Keep inlet immutable
        self.water_props = water_props
        
        # Define stage sequence. Insert placeholders where appropriate.
        # Assuming: superheater after pass1 (hot zone), economizer after pass3 (cool zone).
        # Adjust sequence based on actual boiler design.
        self.stages = [
            cfg_with_calc.stages.pass1,         
            cfg_with_calc.stages.reversal1,
            cfg_with_calc.stages.pass2,
            cfg_with_calc.stages.reversal2,
            cfg_with_calc.stages.pass3,
            cfg_with_calc.stages.economiszer
        ]

    def run_chain(self) -> Tuple[GasStreamProfile, WaterStreamProfile]:
        """
        Run the chained simulation.
        Returns GasStreamProfile and WaterStreamProfile, each containing a list of instances along z.
        """
        gas_points: List[GasStream] = []
        water_points: List[Water] = []
        current_z = 0.0
        current_gas_T = self.gas_inlet.temperature
        current_gas_p = self.gas_inlet.pressure

        # Initial water calc to get starting h
        initial_water_calc = WaterWithCalc(
            water_stream=self.water_inlet,
            water_props=self.water_props
        )
        current_water_h = initial_water_calc.enthalpy

        for stage in self.stages:
            # Handle reversal inlet nozzle if applicable
            nozzle_applied = False
            if isinstance(stage, ReversalWithCalc):
                # Assume K=0.5 for contraction; refine later
                nozzle_in = Nozzle(
                    geom=stage.nozzles.inlet,
                    gas_stream=self.gas_inlet,  # Use inlet for fixed props
                    gas_props=self.gas_props,
                    K=Q_(0.5, ureg.dimensionless)
                )
                T_new, p_new = nozzle_in.apply(
                    current_gas_T,
                    current_gas_p,
                )
                # Add point after nozzle at advanced z
                z_after_nozzle = current_z + stage.nozzles.inlet.length.magnitude
                gas_point_nozzle = GasStream(
                    mass_flow_rate=self.gas_inlet.mass_flow_rate,
                    temperature=Q_(T_new.magnitude, 'kelvin'),
                    pressure=Q_(p_new.magnitude, 'pascal'),
                    composition=self.gas_inlet.composition,
                    spectroscopic_data=self.gas_inlet.spectroscopic_data,
                    z=Q_(z_after_nozzle, 'm'),
                    heat_transfer_rate=None,  # No heat in nozzle
                    velocity=None,
                    reynolds_number=None,
                    density=None
                )
                gas_points.append(gas_point_nozzle)

                # Water point at same z (no change)
                water_point_nozzle = Water(
                    mass_flow_rate=self.water_inlet.mass_flow_rate,
                    temperature=self.water_inlet.temperature,  # Will update later if needed
                    pressure=self.water_inlet.pressure,
                    composition=self.water_inlet.composition,
                    z=Q_(z_after_nozzle, 'm'),
                    enthalpy=current_water_h,
                    heat_transfer_rate=None,
                    quality=None,
                    phase=None,
                    saturation_temperature=None
                )
                water_points.append(water_point_nozzle)

                # Update currents
                current_gas_T = Q_(T_new.magnitude, 'kelvin')
                current_gas_p = Q_(p_new.magnitude, 'pascal')
                current_z = z_after_nozzle
                nozzle_applied = True

            # Create temporary stream objects for calcs
            temp_gas_stream = GasStream(
                mass_flow_rate=self.gas_inlet.mass_flow_rate,
                temperature=current_gas_T,
                pressure=current_gas_p,
                composition=self.gas_inlet.composition,
                spectroscopic_data=self.gas_inlet.spectroscopic_data,
                z=Q_(current_z, 'm')
            )
            temp_water_stream = Water(
                mass_flow_rate=self.water_inlet.mass_flow_rate,
                temperature=self.water_inlet.temperature,  # Temp, updated in loop
                pressure=self.water_inlet.pressure,
                composition=self.water_inlet.composition,
                z=Q_(current_z, 'm')
            )

            # Create stage-specific gas calc
            gas_calc = GasStreamWithCalc(
                gas_props=self.gas_props,
                gas_stream=temp_gas_stream,
                geometry=stage
            )

            water_calc = WaterWithCalc(
                water_stream=temp_water_stream,
                water_props=self.water_props
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
            y0 = [current_gas_T.magnitude, current_gas_p.magnitude, current_water_h.to('J/kg').magnitude]
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
            T_water_stage = []
            for i in range(len(T_stage)):
                ode.gas.gas_stream.temperature = Q_(T_stage[i], 'kelvin')
                # set local water temperature for HeatSystem only (transient)
                hi = Q_(h_stage[i], 'J/kg')
                T_w_local = converter.T_from_h(hi, water_calc.water_stream.pressure)
                water_calc.water_stream.temperature = T_w_local
                Q_dot = ode.heat_system.q_
                Q_dot_stage.append(Q_dot.magnitude)
                T_water_stage.append(T_w_local.magnitude)

                # Temporarily update calcs for additional properties
                temp_gas_stream.temperature = Q_(T_stage[i], 'kelvin')
                temp_gas_stream.pressure = Q_(p_stage[i], 'pascal')
                gas_calc = GasStreamWithCalc(
                    gas_props=self.gas_props,
                    gas_stream=temp_gas_stream,
                    geometry=stage
                )
                water_calc.water_stream.temperature = T_w_local
                water_calc.water_stream.pressure = self.water_inlet.pressure  # Assume constant

                # Create gas point with additional props
                gas_point = GasStream(
                    mass_flow_rate=self.gas_inlet.mass_flow_rate,
                    temperature=Q_(T_stage[i], 'kelvin'),
                    pressure=Q_(p_stage[i], 'pascal'),
                    composition=self.gas_inlet.composition,
                    spectroscopic_data=self.gas_inlet.spectroscopic_data,
                    z=Q_(z_stage[i], 'm'),
                    heat_transfer_rate=Q_(Q_dot_stage[i], 'W/m'),
                    velocity=gas_calc.velocity,
                    reynolds_number=gas_calc.reynolds_number,
                    density=gas_calc.density
                )
                gas_points.append(gas_point)

                # Create water point with additional props
                water_point = Water(
                    mass_flow_rate=self.water_inlet.mass_flow_rate,
                    temperature=T_w_local,
                    pressure=self.water_inlet.pressure,
                    composition=self.water_inlet.composition,
                    z=Q_(z_stage[i], 'm'),
                    enthalpy=hi,
                    heat_transfer_rate=Q_(Q_dot_stage[i], 'W/m'),
                    quality=water_calc.quality,
                    phase=water_calc.phase,
                    saturation_temperature=water_calc.saturation_temperature
                )
                water_points.append(water_point)

            # Update current position and states for next stage
            current_z += z_span[1]
            current_gas_T = Q_(T_stage[-1], 'kelvin')
            current_gas_p = Q_(p_stage[-1], 'pascal')
            current_water_h = Q_(h_stage[-1], 'J/kg')

            # If reversal, handle outlet nozzle similarly
            if isinstance(stage, ReversalWithCalc):
                nozzle_out = Nozzle(
                    geom=stage.nozzles.outlet,
                    gas_stream=self.gas_inlet,
                    gas_props=self.gas_props,
                    K=Q_(0.5, ureg.dimensionless)  # Assume same K
                )
                T_new, p_new = nozzle_out.apply(
                    current_gas_T,
                    current_gas_p,
                )
                z_after_outlet = current_z + stage.nozzles.outlet.length.magnitude
                gas_point_outlet = GasStream(
                    mass_flow_rate=self.gas_inlet.mass_flow_rate,
                    temperature=Q_(T_new.magnitude, 'kelvin'),
                    pressure=Q_(p_new.magnitude, 'pascal'),
                    composition=self.gas_inlet.composition,
                    spectroscopic_data=self.gas_inlet.spectroscopic_data,
                    z=Q_(z_after_outlet, 'm'),
                    heat_transfer_rate=None,
                    velocity=None,
                    reynolds_number=None,
                    density=None
                )
                gas_points.append(gas_point_outlet)

                water_point_outlet = Water(
                    mass_flow_rate=self.water_inlet.mass_flow_rate,
                    temperature=converter.T_from_h(current_water_h, self.water_inlet.pressure),
                    pressure=self.water_inlet.pressure,
                    composition=self.water_inlet.composition,
                    z=Q_(z_after_outlet, 'm'),
                    enthalpy=current_water_h,
                    heat_transfer_rate=None,
                    quality=water_calc.quality,  # Use last calc
                    phase=water_calc.phase,
                    saturation_temperature=water_calc.saturation_temperature
                )
                water_points.append(water_point_outlet)

                current_gas_T = Q_(T_new.magnitude, 'kelvin')
                current_gas_p = Q_(p_new.magnitude, 'pascal')
                current_z = z_after_outlet

        gas_profile = GasStreamProfile(points=gas_points)
        water_profile = WaterStreamProfile(points=water_points)
        return gas_profile, water_profile