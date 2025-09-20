# thermo/materials.py
from typing import Callable
def k_wall_linear(a: float, b: float) -> Callable[[float],float]:
    def k(T: float) -> float:
        return a + b*T
    return k
