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
        # Basic input validation
        if not (isfinite(Di) and isfinite(Do) and Di > 0 and Do > 0 and Do > Di):
            raise ValueError(f"Bad geometry: Di={Di}, Do={Do} (require 0 < Di < Do)")
        if not isfinite(Tg) or not isfinite(Tsat):
            raise ValueError(f"Non-finite temperatures: Tg={Tg}, Tsat={Tsat}")
        if not isfinite(h_g) or not isfinite(h_rad) or not isfinite(h_w):
            raise ValueError(f"Non-finite h: h_g={h_g}, h_rad={h_rad}, h_w={h_w}")
        if h_g + h_rad <= 0 or h_w <= 0:
            raise ValueError(f"Non-positive film coefficient: h_g+h_rad={h_g+h_rad}, h_w={h_w}")
        if not isfinite(R_fg) or not isfinite(R_fw) or R_fg < 0 or R_fw < 0:
            raise ValueError(f"Fouling must be finite and ≥0: R_fg={R_fg}, R_fw={R_fw}")

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
