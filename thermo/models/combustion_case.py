from dataclasses import dataclass
from thermo.core.streams import GasStream

@dataclass(frozen=True)
class CombustionCase:
    air: GasStream
    fuel: GasStream
    excess_air_ratio: float
    T_ref: float
