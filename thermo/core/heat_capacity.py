import CoolProp.CoolProp as CP
from scipy.integrate import quad
from typing import Dict
from common.units import ureg, Q_


class MixtureCp:
    def __init__(self, fluid_map: dict):
        self._map = fluid_map

    def cp_mass_mixture(self, T_K: Q_, P_Pa: Q_, mass_fractions: Dict[str, Q_]) -> Q_:
        T_val = T_K.to("kelvin").magnitude
        P_val = P_Pa.to("pascal").magnitude

        cp_mix_kJ_per_kgK = 0.0
        for fluid, w in mass_fractions.items():
            w_val = w.to("").magnitude if hasattr(w, "to") else float(w)
            if fluid == "H2O":
                AS = CP.AbstractState("HEOS", "Water")
                try:
                    AS.update(CP.PT_INPUTS, P_val, T_val)
                except ValueError:
                    AS.specify_phase(CP.iphase_gas)
                    AS.update(CP.PT_INPUTS, P_val, T_val + 1e-3)
                cp_i = AS.cpmass()
            else:
                AS = CP.AbstractState("HEOS", self._map[fluid])
                AS.update(CP.PT_INPUTS, P_val, T_val)
                cp_i = AS.cpmass()
            cp_mix_kJ_per_kgK += w_val * (cp_i / 1000.0)

        return cp_mix_kJ_per_kgK * ureg.kilojoule / (ureg.kilogram * ureg.kelvin)

    def integrate_cp_mass(self, P_Pa: Q_, mass_fractions: Dict[str, Q_], T1: Q_, T2: Q_) -> Q_:
        P_val = P_Pa.to("pascal").magnitude
        T1_val = T1.to("kelvin").magnitude
        T2_val = T2.to("kelvin").magnitude

        result = quad(
            lambda T: self.cp_mass_mixture(Q_(T, ureg.kelvin), Q_(P_val, ureg.pascal), mass_fractions).to(
                "kilojoule/(kilogram*kelvin)"
            ).magnitude,
            T1_val,
            T2_val,
        )[0]

        return result * ureg.kilojoule / ureg.kilogram
