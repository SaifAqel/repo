# models/heat_balance.py
from dataclasses import dataclass
from typing import Callable
from math import isfinite, log

EPS = 1e-12

@dataclass
class LocalHeatBalance:
    k_wall: Callable[[float], float]

    def solve(
        self, h_g: float, h_rad: float, h_w: float,
        R_fg: float, R_fw: float, Di: float, Do: float,
        Tg: float, Tsat: float, guess_Twi: float, guess_Two: float
    ) -> tuple:       
        Twi, Two = guess_Twi, guess_Two
        for _ in range(50):
            Tm = (Twi + Two) / 2.0
            k = self.k_wall(Tm)
            if not isfinite(k) or k <= 0:
                raise ValueError(f"k_wall(T={Tm}) invalid: {k}")

            # Resistances (per same normalization as your original)
            Ri = 1.0 / max(h_w, EPS) + R_fw
            # Correct cylindrical wall resistance with πL normalized out
            Rw = log(Do / Di) / (2.0 * max(k, EPS))
            Ro = 1.0 / max(h_g + h_rad, EPS) + R_fg

            denom = Ro + Rw + Ri
            if not isfinite(denom) or denom <= 0:
                raise FloatingPointError(f"Denominator invalid: Ro={Ro}, Rw={Rw}, Ri={Ri}, sum={denom}")

            qpp = (Tg - Tsat) / denom
            Twi_new = Tsat + qpp * Ri
            Two_new = Tg - qpp * Ro

            if isfinite(Twi_new) and isfinite(Two_new) and abs(Twi_new - Twi) < 1e-6 and abs(Two_new - Two) < 1e-6:
                return qpp, Twi_new, Two_new

            if not isfinite(Twi_new) or not isfinite(Two_new):
                raise FloatingPointError(f"Iterate invalid: Twi_new={Twi_new}, Two_new={Two_new}, qpp={qpp}")

            Twi, Two = Twi_new, Two_new

        # Not converged
        return qpp, Twi, Two
