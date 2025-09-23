# core/types.py
from dataclasses import dataclass
from typing import Dict, List
from common.units import Q_

@dataclass(frozen=True)
class GasState:
    T: Q_          # K
    P: Q_          # Pa
    m_dot: Q_      # kg/s
    y: Dict[str, Q_]   # dimensionless
    rho: Q_        # kg/m^3
    cp: Q_         # J/kg/K
    mu: Q_         # Pa*s
    k: Q_          # W/m/K
    Pr: Q_         # 1
    u: Q_          # m/s

@dataclass(frozen=True)
class CellResult:
    x: Q_          # m
    Tg: Q_         # K
    Pg: Q_         # Pa
    rhog: Q_       # kg/m^3
    cp: Q_         # J/kg/K
    mu: Q_         # Pa*s
    k: Q_          # W/m/K
    Pr: Q_         # 1
    ug: Q_         # m/s
    Twi: Q_        # K
    Two: Q_        # K
    qpp: Q_        # W/m^2
    h_g: Q_        # W/m^2/K
    eps_g: Q_      # 1

@dataclass(frozen=True)
class StageResult:
    cells: List[CellResult]
    Q: Q_          # W
    dP: Q_         # Pa
    T_out: Q_      # K
    P_out: Q_      # Pa
    Tw_min: Q_     # K
    Tw_max: Q_     # K
    h_g_avg: Q_    # W/m^2/K
    eps_g_avg: Q_  # 1

@dataclass(frozen=True)
class PlantResult:
    stages: List[StageResult]
    Q_total: Q_    # W
    dP_total: Q_   # Pa
    T_final: Q_    # K
    P_final: Q_    # Pa
