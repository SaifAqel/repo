# validation/budgets.py
def energy_balance(stage_result, m_dot, cp, Tin, Tout):
    Qcells = sum(c.qpp for c in stage_result.cells)*(stage_result.cells[0].x*0+1)
    return {"Q_cells": Qcells, "Q_bulk": m_dot*cp*(Tin-Tout)}
