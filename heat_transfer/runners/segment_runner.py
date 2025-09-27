from __future__ import annotations
from dataclasses import dataclass, replace
from typing import Callable, Any, Dict
from common.units import Q_

# --- State container for a single segment ---
@dataclass
class SegmentState:
    T_gas_in: Q_
    p_gas_in: Q_
    T_water_in: Q_
    T_wall_inner: Q_
    T_wall_outer: Q_
    T_gas_bulk: Q_
    p_gas_out: Q_
    T_gas_out: Q_
    q_segment: Q_
    h_conv: Q_
    h_rad: Q_
    R_total: Q_
    Re: Q_
    Pr: Q_
    Nu: Q_
    v: Q_
    f: Q_
    dp: Q_

# --- Iteration driver (no validations, no assumptions) ---
class SegmentEquilibrium:
    """
    Minimal fixed-point iterator for one PassWithCalc segment.

    Required callables (you provide these; this class only loops):
      - make_gas_stream(state) -> GasStreamWithCalc
      - radiation_h(state) -> Q_          # returns h_rad
      - thermal_resistance(state) -> Q_   # returns total area-specific R [K*m^2/W]
      - pressure_drop(state) -> Dict[str, Q_] with keys {"f","dp"}  # Darcy f and Δp
      - update_wall_temperatures(state, q) -> Dict[str, Q_] with keys {"T_wall_inner","T_wall_outer"}

    Geometry and constants are provided via 'env' dict, which must include:
      {"A_segment": Q_, "sigma": Q_, "L_segment": Q_}
    """
    def __init__(self,
                 env: Dict[str, Q_],
                 make_gas_stream: Callable[[SegmentState], Any],
                 radiation_h: Callable[[SegmentState], Q_],
                 thermal_resistance: Callable[[SegmentState], Q_],
                 pressure_drop: Callable[[SegmentState], Dict[str, Q_]],
                 update_wall_temperatures: Callable[[SegmentState, Q_], Dict[str, Q_]]):
        self.env = env
        self.make_gas_stream = make_gas_stream
        self.radiation_h = radiation_h
        self.thermal_resistance = thermal_resistance
        self.pressure_drop = pressure_drop
        self.update_wall_temperatures = update_wall_temperatures

    def run(self,
            init: SegmentState,
            tol_T: Q_,
            tol_p: Q_,
            max_iter: int) -> SegmentState:

        state = replace(init)

        for _ in range(max_iter):
            # 1) Update gas-side transport using user-provided GasStreamWithCalc
            gs = self.make_gas_stream(state)
            v = gs.velocity
            Re = gs.reynolds_number
            Pr = gs.prandt_number
            Nu = gs.nusselt_number
            h_conv = gs.Convective_coefficient

            # 2) Radiation coefficient from user callback
            h_rad = self.radiation_h(state)

            # 3) Total thermal resistance (area-specific)
            R_total = self.thermal_resistance(replace(state, h_conv=h_conv, h_rad=h_rad))

            # 4) Heat rate over this segment
            A = self.env["A_segment"]
            dT = state.T_gas_bulk - state.T_water_in
            q_segment = (A * dT / R_total)

            # 5) Update wall temperatures from supplied relation
            w = self.update_wall_temperatures(state, q_segment)
            T_wi = w["T_wall_inner"]
            T_wo = w["T_wall_outer"]

            # 6) Gas temperature drop by energy balance
            m_dot = gs.m_dot
            cp_g = gs.cp
            T_gas_out = state.T_gas_in - q_segment / (m_dot * cp_g)

            # 7) Pressure drop from user callback
            pd = self.pressure_drop(replace(state, v=v, Re=Re))
            f = pd["f"]
            dp = pd["dp"]
            p_out = state.p_gas_in - dp

            # 8) Build next state
            next_state = SegmentState(
                T_gas_in=state.T_gas_in,
                p_gas_in=state.p_gas_in,
                T_water_in=state.T_water_in,
                T_wall_inner=T_wi,
                T_wall_outer=T_wo,
                T_gas_bulk=(state.T_gas_in + T_gas_out)/2,
                p_gas_out=p_out,
                T_gas_out=T_gas_out,
                q_segment=q_segment,
                h_conv=h_conv,
                h_rad=h_rad,
                R_total=R_total,
                Re=Re,
                Pr=Pr,
                Nu=Nu,
                v=v,
                f=f,
                dp=dp,
            )

            # 9) Convergence check (tolerances provided by caller; no hardcoding)
            dT_max = max(
                abs((next_state.T_gas_out - state.T_gas_out).to(tol_T.units)),
                abs((next_state.T_wall_inner - state.T_wall_inner).to(tol_T.units)),
                abs((next_state.T_wall_outer - state.T_wall_outer).to(tol_T.units)),
            )
            dp_res = abs((next_state.p_gas_out - state.p_gas_out).to(tol_p.units))

            state = next_state
            if dT_max <= tol_T and dp_res <= tol_p:
                break

        return state
