# heat_transfer/solver/water_balance.py
from dataclasses import dataclass
from common.units import Q_, ureg
from heat_transfer.calc_ops.stream_with_calc import WaterWithCalc
from scipy.optimize import root_scalar

@dataclass
class WaterStateConverter:
    water_with_calc: WaterWithCalc

    def h_from_T(self, T: Q_, P: Q_) -> Q_:
        Tsat = self.water_with_calc.water_props.Tsat(P)
        if T <= Tsat:
            return self.water_with_calc.water_props.h_l(T)
        return self.water_with_calc.water_props.h_v(T)

    def T_from_h(self, h: Q_, P: Q_) -> Q_:
        """
        Convert enthalpy h at pressure P to temperature T.
        Handles:
        - Subcooled liquid
        - Saturated mixture
        - Superheated vapor
        """
        # Saturation properties
        Tsat = self.water_with_calc.water_props.Tsat(P)
        h_l_sat = self.water_with_calc.water_props.h_l(Tsat)
        h_v_sat = self.water_with_calc.water_props.h_v(Tsat)

        # Saturated mixture
        if h_l_sat <= h <= h_v_sat:
            return Tsat

        margin = 5.0  # K, keep away from phase boundaries

        # ---------------- Subcooled liquid ----------------
        if h < h_l_sat:
            def f(Tk):
                return (self.water_with_calc.water_props.h_l(Q_(Tk, 'kelvin')).to('J/kg') - h.to('J/kg')).magnitude

            low = max(273.15 + 1e-2, Tsat.to('kelvin').magnitude - 500)
            high = Tsat.to('kelvin').magnitude - margin

            f_low, f_high = f(low), f(high)
            # If root bracket invalid, return nearest boundary
            if f_low * f_high > 0:
                return Q_(low if abs(f_low) < abs(f_high) else high, 'kelvin')

            sol = root_scalar(f, bracket=[low, high])
            return Q_(sol.root, 'kelvin')

        # ---------------- Superheated vapor ----------------
        if h > h_v_sat:
            def f2(Tk):
                return (self.water_with_calc.water_props.h_v(Q_(Tk, 'kelvin')).to('J/kg') - h.to('J/kg')).magnitude

            low = Tsat.to('kelvin').magnitude + margin
            high = min(1073.15 - margin, low + 500)

            f_low, f_high = f2(low), f2(high)
            if f_low * f_high > 0:
                return Q_(low if abs(f_low) < abs(f_high) else high, 'kelvin')

            sol = root_scalar(f2, bracket=[low, high])
            return Q_(sol.root, 'kelvin')




        def quality_from_h(self, h: Q_, P: Q_) -> Q_:
            Tsat = self.water_with_calc.water_props.Tsat(P)
            h_l = self.water_with_calc.water_props.h_l(Tsat)
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
