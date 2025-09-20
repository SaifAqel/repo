# core/types.py
from dataclasses import dataclass
from typing import Dict, List
@dataclass(frozen=True)
class GasState:
    T: float
    P: float
    m_dot: float
    y: Dict[str,float]
    rho: float
    cp: float
    mu: float
    k: float
    Pr: float
    u: float
@dataclass(frozen=True)
class WallState:
    Ti: float
    To: float
    k_wall: float
@dataclass(frozen=True)
class WaterState:
    P_drum: float
    Tsat: float
    hw: float
    quality: float
@dataclass(frozen=True)
class CellResult:
    x: float
    Tg: float
    Pg: float
    rhog: float
    ug: float
    cp: float
    mu: float
    k: float
    Pr: float
    h_g: float
    eps_g: float
    h_rad: float
    h_w: float
    qpp: float
    Twi: float
    Two: float
@dataclass(frozen=True)
class StageResult:
    cells: List[CellResult]
    Q: float
    dP: float
    T_out: float
    P_out: float
    Tw_min: float
    Tw_max: float
    h_g_avg: float
    eps_g_avg: float
@dataclass(frozen=True)
class PlantResult:
    stages: List[StageResult]
    Q_total: float
    dP_total: float
    T_final: float
    P_final: float
