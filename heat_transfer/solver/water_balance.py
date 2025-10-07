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
        # Saturation at given P (10 bar)
        Tsat = self.water_with_calc.water_props.Tsat(P)
        h_l = self.water_with_calc.water_props.h_l(Tsat)
        h_v = self.water_with_calc.water_props.h_v(Tsat)

        # Saturated mixture
        if h_l <= h <= h_v:
            return Tsat

        margin = 5.0  # K, keep away from phase boundaries

        if h < h_l:
            # Subcooled liquid
            def f(Tk):
                return (self.water_with_calc.water_props.h_l(Q_(Tk, 'kelvin')).to('J/kg') - h.to('J/kg')).magnitude

            low = 273.15 + margin
            high = Tsat.to('kelvin').magnitude - margin
            if high <= low:
                high = Tsat.to('kelvin').magnitude - 1e-2
                low = 273.15 + 1e-2
            sol = root_scalar(f, bracket=[low, high])
            return Q_(sol.root, 'kelvin')

        # Superheated vapor
        def f2(Tk):
            return (self.water_with_calc.water_props.h_v(Q_(Tk, 'kelvin')).to('J/kg') - h.to('J/kg')).magnitude

        low = Tsat.to('kelvin').magnitude + margin
        high = min(1073.15 - margin, low + 500)
        if high <= low:
            high = low + 10  # fallback
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
