# models/heat_balance.py
from dataclasses import dataclass
from typing import Callable, Tuple
import numpy as np

from common.units import ureg, Q_

TOL_T = Q_(1e-6, ureg.kelvin)

@dataclass
class LocalHeatBalance:
    k_wall: Callable  # callable that returns thermal conductivity

    def solve(
        self,
        h_g: float,            # W/(m^2*K)
        h_rad: float,          # W/(m^2*K)
        h_w: float,            # W/(m^2*K)
        R_fg: float,           # K*m^2/W  (fin/gas-side extra resistance per area)
        R_fw: float,           # K*m^2/W  (fouling/water-side resistance per area)
        Di: float,             # m
        Do: float,             # m
        Tg: float,             # K
        Tsat: float,           # K
        guess_Twi: float,      # K
        guess_Two: float       # K
    ) -> Tuple:
        # Promote all inputs to pint quantities
        h_g   = Q_(h_g,   ureg.watt / (ureg.meter**2 * ureg.kelvin))
        h_rad = Q_(h_rad, ureg.watt / (ureg.meter**2 * ureg.kelvin))
        h_w   = Q_(h_w,   ureg.watt / (ureg.meter**2 * ureg.kelvin))
        R_fg  = Q_(R_fg,  ureg.kelvin * ureg.meter**2 / ureg.watt)
        R_fw  = Q_(R_fw,  ureg.kelvin * ureg.meter**2 / ureg.watt)
        Di    = Q_(Di,    ureg.meter)
        Do    = Q_(Do,    ureg.meter)
        Tg    = Q_(Tg,    ureg.kelvin)
        Tsat  = Q_(Tsat,  ureg.kelvin)
        Twi   = Q_(guess_Twi, ureg.kelvin)
        Two   = Q_(guess_Two, ureg.kelvin)

        # Tiny floor values with proper units
        eps_h = Q_(1e-24, ureg.watt / (ureg.meter**2 * ureg.kelvin))
        eps_k = Q_(1e-30, ureg.watt / (ureg.meter * ureg.kelvin))

        for _ in range(50):
            Tm = (Twi + Two) / Q_(2.0, ureg.dimensionless)

            # Evaluate k_wall; accept either quantity or float return
            k_val = self.k_wall(Tm)
            if isinstance(k_val, ureg.Quantity):
                k = k_val.to(ureg.watt / (ureg.meter * ureg.kelvin))
            else:
                k = Q_(k_val, ureg.watt / (ureg.meter * ureg.kelvin))

            if not np.isfinite(k.magnitude) or k <= eps_k:
                raise ValueError(f"k_wall(T={Tm}) invalid: {k}")

            # Resistances per unit area
            Ri = (Q_(1.0, ureg.dimensionless) / max(h_w, eps_h)) + R_fw
            # Cylindrical wall resistance with πL normalized out: ln(Do/Di)/(2*k)
            ln_ratio = np.log((Do / Di).to(ureg.dimensionless).magnitude)
            Rw = Q_(ln_ratio, ureg.dimensionless) / (Q_(2.0, ureg.dimensionless) * k)
            Ro = (Q_(1.0, ureg.dimensionless) / max(h_g + h_rad, eps_h)) + R_fg

            denom = Ro + Rw + Ri
            if not np.isfinite(denom.magnitude) or denom <= Q_(0.0, denom.units):
                raise FloatingPointError(
                    f"Denominator invalid: Ro={Ro}, Rw={Rw}, Ri={Ri}, sum={denom}"
                )

            qpp = (Tg - Tsat) / denom  # W/m^2
            Twi_new = Tsat + qpp * Ri
            Two_new = Tg   - qpp * Ro

            if (
                np.isfinite(Twi_new.magnitude)
                and np.isfinite(Two_new.magnitude)
                and abs((Twi_new - Twi).to(ureg.kelvin)) < TOL_T
                and abs((Two_new - Two).to(ureg.kelvin)) < TOL_T
            ):
                return qpp, Twi_new, Two_new

            if not np.isfinite(Twi_new.magnitude) or not np.isfinite(Two_new.magnitude):
                raise FloatingPointError(
                    f"Iterate invalid: Twi_new={Twi_new}, Two_new={Two_new}, qpp={qpp}"
                )

            Twi, Two = Twi_new, Two_new

        # Not converged
        return qpp, Twi, Two
