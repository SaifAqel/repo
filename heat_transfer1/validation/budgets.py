# validation/budgets.py
from common.units import ureg, Q_

def energy_balance(stage_result, m_dot, cp, Tin, Tout):
    m_dot = Q_(m_dot, ureg.kg / ureg.second)
    cp    = Q_(cp, ureg.joule / (ureg.kg * ureg.kelvin))
    Tin   = Q_(Tin, ureg.kelvin)
    Tout  = Q_(Tout, ureg.kelvin)

    Qcells = sum(c.qpp for c in stage_result.cells) * Q_(1.0, ureg.dimensionless)
    Qbulk  = m_dot * cp * (Tin - Tout)

    return {"Q_cells": Qcells.to(ureg.watt), "Q_bulk": Qbulk.to(ureg.watt)}
