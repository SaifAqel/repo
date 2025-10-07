# heat_transfer/solver/water_balance.py
from dataclasses import dataclass
from common.units import Q_, ureg
from heat_transfer.calc_ops.stream_with_calc import WaterWithCalc
from scipy.optimize import root_scalar

from scipy.optimize import root_scalar


@dataclass
class WaterStateConverter:
    water_with_calc: WaterWithCalc

    def h_from_T(self, T: Q_, P: Q_) -> Q_:
        Tsat = self.water_with_calc.water_props.Tsat(P)
        if T <= Tsat:
            return self.water_with_calc.water_props.h_l_sat(P)
        return self.water_with_calc.water_props.h_v_sat(P)

    def T_from_h(self, h: Q_, P: Q_) -> Q_:
        """
        Convert enthalpy h at pressure P to temperature T.
        Handles subcooled liquid, saturated mixture, and superheated vapor.
        Uses IAPWS97(P, T) single-phase solve. Returns Tsat(P) if two-phase.
        Upper bracket extended to 2273 K.
        """
        Tsat = self.water_with_calc.water_props.Tsat(P)
        h_l_sat = self.water_with_calc.water_props.h_l_sat(P)
        h_v_sat = self.water_with_calc.water_props.h_v_sat(P)

        # saturated mixture
        if h_l_sat <= h <= h_v_sat:
            return Tsat

        margin = 5.0  # K

        # ---------------- Subcooled liquid ----------------
        if h < h_l_sat:
            def f(Tk):
                target = self.water_with_calc.water_props.h_single(P, Q_(Tk, 'kelvin')).to('J/kg')
                return (target.magnitude - h.to('J/kg').magnitude)

            low = max(273.15 + 1e-2, (Tsat.to('kelvin').magnitude - 500))
            high = Tsat.to('kelvin').magnitude - margin

            f_low, f_high = f(low), f(high)
            if f_low * f_high > 0:
                return Q_(low if abs(f_low) < abs(f_high) else high, 'kelvin')

            sol = root_scalar(f, bracket=[low, high])
            return Q_(sol.root, 'kelvin')

        # ---------------- Superheated vapor ----------------
        if h > h_v_sat:
            def f2(Tk):
                target = self.water_with_calc.water_props.h_single(P, Q_(Tk, 'kelvin')).to('J/kg')
                return (target.magnitude - h.to('J/kg').magnitude)

            low = Tsat.to('kelvin').magnitude + margin
            high = min(2273.0 - margin, low + 1500.0)

            f_low, f_high = f2(low), f2(high)
            if f_low * f_high > 0:
                return Q_(low if abs(f_low) < abs(f_high) else high, 'kelvin')

            sol = root_scalar(f2, bracket=[low, high])
            return Q_(sol.root, 'kelvin')

    def quality_from_h(self, h: Q_, P: Q_) -> Q_:
        Tsat = self.water_with_calc.water_props.Tsat(P)
        h_l = self.water_with_calc.water_props.h_l_sat(P)
        h_fg = self.water_with_calc.water_props.h_fg(P)
        return (h - h_l) / h_fg



@dataclass
class WaterEnergyRHS:
    water_with_calc: WaterWithCalc
    geom: object

    def m_dot_per_circuit(self):
        m = self.water_with_calc.water_stream.mass_flow_rate
        circuits = getattr(self.geom.geometry, "circuit_count", None)
        if circuits is None:
            circuits = getattr(self.geom.geometry, "number_of_tubes", 1)
        return (m / circuits).to("kg / s")

    def dhdz(self, h_local: Q_, q_per_length: Q_) -> Q_:
        m_dot = self.m_dot_per_circuit()
        return (q_per_length / m_dot).to("J / (kg * m)")
