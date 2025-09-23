# thermo/provider.py
from abc import ABC, abstractmethod
from typing import Dict, List, Tuple
import CoolProp.CoolProp as CP
from functools import lru_cache
from common.units import Q_, ureg

_ALIAS = {
    "CH4":"Methane","CO2":"CarbonDioxide","N2":"Nitrogen","O2":"Oxygen",
    "H2O":"Water","H2":"Hydrogen","AR":"Argon","Ar":"Argon",
}

class IThermoProvider(ABC):
    @abstractmethod
    def gas_props(self, T: Q_, P: Q_, y: Dict[str, Q_]) -> Dict[str, Q_]: ...
    @abstractmethod
    def water_saturation(self, P: Q_) -> Dict[str, Q_]: ...

class CoolPropThermoProvider(IThermoProvider):
    def __init__(self, backend: str = "HEOS"):
        self.backend = backend

    @staticmethod
    def _alias_and_vector(y: Dict[str, Q_]) -> Tuple[List[str], List[float]]:
        species: List[str] = []
        z: List[float] = []
        for k, v in y.items():
            name = _ALIAS.get(k, k)
            val = float(v.to("dimensionless").magnitude) if hasattr(v, "to") else float(v)
            species.append(name); z.append(val)
        s = sum(z)
        if s <= 0: raise ValueError("sum(y) must be > 0")
        z = [max(0.0, zi/s) for zi in z]
        return species, z

    def gas_props(self, T: Q_, P: Q_, y: Dict[str, Q_]) -> Dict[str, Q_]:
        sp, z = self._alias_and_vector(y)
        st = CP.AbstractState(self.backend, "&".join(sp))
        st.set_mole_fractions(z)
        try:
            st.specify_phase(CP.iphase_gas)
        except Exception:
            pass
        try:
            st.update(CP.PT_INPUTS, P.to("Pa").magnitude, T.to("K").magnitude)
        except ValueError:
            rho_guess = 1.0
            st.update(CP.DmassT_INPUTS, rho_guess, T.to("K").magnitude)
        return {
            "rho": Q_(st.rhomass(), "kg/m^3"),
            "cp":  Q_(st.cpmass(),  "J/kg/K"),
            "mu":  Q_(st.viscosity(), "Pa*s"),
            "k":   Q_(st.conductivity(), "W/m/K"),
            "Pr":  Q_(st.Prandtl(), "dimensionless"),
        }

    def water_saturation(self, P: Q_) -> Dict[str, Q_]:
        Tsat = Q_(CP.PropsSI("T","P",P.to("Pa").magnitude,"Q",0,"Water"), "K")
        return {"Tsat": Tsat}
