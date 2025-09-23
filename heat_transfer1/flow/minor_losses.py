# flow/minor_losses.py
from typing import List, Tuple
from common.units import ureg, Q_

def K_total(elements: List[Tuple[str, Q_]]) -> Q_:
    return sum(k.to("dimensionless") for _, k in elements)
