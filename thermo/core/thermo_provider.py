from typing import Protocol

class IThermoProvider(Protocol):
    def cp_mass(self, T_K: float, P_Pa: float, species: str) -> float:
        """Return J/(kg·K). Exact behavior matches legacy CoolProp usage."""
