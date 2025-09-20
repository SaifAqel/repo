# flow/minor_losses.py
from typing import List, Tuple
def K_total(elements: List[Tuple[str,float]]) -> float:
    return sum(k for _,k in elements)
