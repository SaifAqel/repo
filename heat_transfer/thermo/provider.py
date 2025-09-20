# thermo/provider.py
from abc import ABC, abstractmethod
from typing import Dict, Tuple, List
import CoolProp.CoolProp as CP

_ALIAS = { "CH4": "Methane", "CO2": "CarbonDioxide", "N2": "Nitrogen", "O2": "Oxygen", } 

class IThermoProvider(ABC): 
    @abstractmethod 
    def gas_props(self, T: float, P: float, y: Dict[str,float]) -> Dict[str,float]: ... 
    @abstractmethod 
    def water_saturation(self, P: float) -> Dict[str,float]: ...

class CoolPropThermoProvider(IThermoProvider):
    def __init__(self, backend: str = "HEOS"):
        self.backend = backend  # "HEOS" preferred for transport. "PR"/"SRK" OK for density/rough transport.

    def _canon(self, s: str) -> str:
        return _ALIAS.get(s, s)

    def water_saturation(self, P: float) -> Dict[str, float]:
        Tsat  = CP.PropsSI("T", "P", P, "Q", 0, "Water")
        rho_l = CP.PropsSI("D", "P", P, "Q", 0, "Water")
        rho_v = CP.PropsSI("D", "P", P, "Q", 1, "Water")
        h_fg  = CP.PropsSI("H", "P", P, "Q", 1, "Water") - CP.PropsSI("H", "P", P, "Q", 0, "Water")
        return {"Tsat": Tsat, "rho_l": rho_l, "rho_v": rho_v, "h_fg": h_fg}


    def _mix_lists(self, y: Dict[str,float]) -> Tuple[List[str], List[float]]:
        # deterministic order + canonical names, normalized x
        items = sorted((self._canon(k), v) for k, v in y.items())
        fluids = [k for k, _ in items]
        xs_raw = [v for _, v in items]
        s = sum(xs_raw)
        xs = [v/s for v in xs_raw] if s > 0 else xs_raw
        return fluids, xs

    def _state(self, backend: str, y: Dict[str,float]):
        fluids, xs = self._mix_lists(y)
        AS = CP.AbstractState(backend, "&".join(fluids))
        AS.set_mole_fractions(xs)
        # Force gas to avoid phase-envelope search/critical finder
        AS.specify_phase(CP.iphase_gas)
        return AS

    def _gas_props_with_backend(self, backend: str, T: float, P: float, y: Dict[str,float]) -> Dict[str,float]:
        AS = self._state(backend, y)
        AS.update(CP.PT_INPUTS, P, T)
        rho = AS.rhomass()                     # kg/m^3
        cp  = AS.cpmass()                      # J/(kg·K)
        mu  = AS.viscosity()                   # Pa·s
        k   = AS.conductivity()                # W/(m·K)
        Pr  = cp * mu / k                      # dimensionless
        return {"rho": rho, "cp": cp, "mu": mu, "k": k, "Pr": Pr}

    def _ideal_gas_density(self, T: float, P: float, y: Dict[str,float]) -> float:
        fluids, xs = self._mix_lists(y)
        Ms = [CP.PropsSI("M", f) for f in fluids]  # kg/mol
        Mmix = sum(x*m for x, m in zip(xs, Ms))
        R = 8.314462618   #J/(mol*K)
        return P * Mmix / (R * T)


    def gas_props(self, T: float, P: float, y: Dict[str,float]) -> Dict[str,float]:
        # Try HEOS first with forced gas phase, then SRK, then PR. Last resort: ideal-gas rho + HEOS transport if possible.
        backends = [self.backend, "SRK", "PR"]
        last_err = None
        for be in backends:
            try:
                return self._gas_props_with_backend(be, T, P, y)
            except Exception as e:
                last_err = e

        # Final fallback: ideal-gas density + best-effort transport from HEOS
        try:
            rho = self._ideal_gas_density(T, P, y)
            props = self._gas_props_with_backend("HEOS", T, P, y)
            props["rho"] = rho
            props["Pr"]  = props["cp"] * props["mu"] / props["k"]
            return props
        except Exception:
            # If even transport fails, return only rho_ideal with minimal placeholders
            rho = self._ideal_gas_density(T, P, y)
            return {"rho": rho, "cp": float("nan"), "mu": float("nan"), "k": float("nan"), "Pr": float("nan")}
