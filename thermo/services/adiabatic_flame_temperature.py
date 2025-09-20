from thermo.core.units import ureg, Q_
from thermo.core.heat_capacity import MixtureCp

class AdiabaticFlameTemperature:
    def __init__(self, cp: "MixtureCp", solver):
        self._cp = cp
        self._solver = solver

    def _flue_sensible_enthalpy(self, T_K: Q_, P_Pa: Q_, mass_frac: dict, 
                                mass_flow_kg_s: Q_, T_ref_K: Q_) -> Q_:
        # integrate_cp_mass should return Δh [J/kg]
        dh = self._cp.integrate_cp_mass(P_Pa, mass_frac, T_ref_K, T_K)  # no extra units
        return (mass_flow_kg_s * dh).to(ureg.watt)

    def _residual(self, T_K: Q_, P_Pa: Q_, mass_frac: dict, mass_flow_kg_s: Q_, 
                  Q_in: Q_, T_ref_K: Q_) -> Q_:
        h_flue = self._flue_sensible_enthalpy(T_K, P_Pa, mass_frac, mass_flow_kg_s, T_ref_K)
        return (Q_in.to(ureg.watt) - h_flue).to(ureg.watt)

    def solve(self, P_Pa: Q_, flue_mass_fracs: dict, flue_mass_flow: Q_, 
              Q_in: Q_, T_ref_K: Q_) -> Q_:
        f = lambda T: self._residual(T * ureg.kelvin, P_Pa, flue_mass_fracs, 
                                     flue_mass_flow, Q_in, T_ref_K).to(ureg.watt).m
        T_val = self._solver(f, (1700, 3000), (), 1e-6)
        return T_val * ureg.kelvin
