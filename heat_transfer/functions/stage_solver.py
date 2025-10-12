# heat_transfer/functions/stage_solver.py
from __future__ import annotations
from dataclasses import dataclass
from math import pi, log
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import solve_ivp

from common.units import Q_, ureg
from heat_transfer.config.models import (
    FirePass, SmokePass, Reversal, BankGeometry, GasStream, WaterStream
)
from heat_transfer.functions.htc_water import HTCFunctions


def _areas(stage) -> tuple[Q_, Q_, Q_, Q_]:
    ri = (stage.hot_side.inner_diameter / 2).to("m")
    ro = (stage.hot_side.outer_diameter / 2).to("m")
    L  = stage.hot_side.inner_length.to("m")
    Ai = 2 * pi * ri * L
    Ao = 2 * pi * ro * L
    At = Ai
    if isinstance(stage.hot_side, BankGeometry):
        Ai *= stage.hot_side.tubes_number.to("").magnitude
        Ao *= stage.hot_side.tubes_number.to("").magnitude
        At  = Ai
    return ri, ro, Ai, Ao


# ... imports unchanged ...

def _overall_Q_for_T(Tg, Pg, Pw, hw, stage, gas, water) -> tuple[Q_, Q_, Q_]:
    ri, ro, Ai, Ao = _areas(stage)
    L  = stage.hot_side.inner_length
    kw = stage.hot_side.wall.conductivity.to("W/(m*K)")
    Dh_shell = stage.cold_side.hydraulic_diameter

    # local states from ODE
    TgQ = Q_(Tg, "K")
    PgQ = Q_(Pg, "Pa")
    PwQ = Q_(Pw, "Pa")
    hwQ = Q_(hw, "J/kg")

    # sync stream objects to local state BEFORE any property is queried
    gas.temperature = TgQ
    gas.pressure    = PgQ
    water.enthalpy  = hwQ

    Tw_bulk = water.water_props.T_ph(PwQ, hwQ)

    h_g = gas.convective_coefficient.to("W/(m^2*K)")

    def F(TwiK: float) -> float:
        Twi = Q_(TwiK, "K")
        h_r = gas.radiation_coefficient(Twi).to("W/(m^2*K)")
        qph = 2 * pi * ri * (h_g + h_r) * (TgQ - Twi)                # W/m
        Two = Twi - qph * log(ro/ri) / (2 * pi * kw)                 # K
        qpp = (qph / (2 * pi * ro)).to("W/m^2")
        h_w = HTCFunctions.shell(water, Dh_shell, L, Two, qpp)       # W/m^2-K
        qpc = 2 * pi * ro * h_w * (Two - Tw_bulk)                    # W/m
        if isinstance(stage.hot_side, BankGeometry):
            qpc *= stage.hot_side.tubes_number.to("").magnitude
        return (qph - qpc).to("W/m").magnitude

    a = Tw_bulk.to("K").magnitude + 1e-6
    b = TgQ.to("K").magnitude - 1e-6
    xs = np.linspace(a, b, 17)
    Fs = [F(x) for x in xs]
    i0 = next((i for i in range(len(xs)-1) if Fs[i]*Fs[i+1] <= 0.0), 0)
    Twi = Q_(brentq(F, xs[i0], xs[i0+1]), "K")

    h_r = gas.radiation_coefficient(Twi).to("W/(m^2*K)")
    qph = 2 * pi * ri * (h_g + h_r) * (TgQ - Twi)                    # W/m
    Two = Twi - qph * log(ro/ri) / (2 * pi * kw)                     # K
    if isinstance(stage.hot_side, BankGeometry):
        qph *= stage.hot_side.tubes_number.to("").magnitude
    return qph.to("W/m"), Twi, Two



def _dp_gas(stage, gas) -> Q_:
    L  = stage.hot_side.inner_length.to("m")
    D  = stage.hot_side.hydraulic_diameter.to("m")
    v  = gas.velocity.to("m/s")
    rho = gas.density.to("kg/m^3")
    f = gas.friction_factor
    dp_fric = f * (L / D) * 0.5 * rho * v * v                         # Pa
    k_in  = getattr(getattr(stage.hot_side, "nozzles", None), "inlet",  None)
    k_out = getattr(getattr(stage.hot_side, "nozzles", None), "outlet", None)
    Ksum = 0.0
    if k_in is not None:  Ksum += k_in.k.to("").magnitude
    if k_out is not None: Ksum += k_out.k.to("").magnitude
    dp_k = Q_(Ksum, "") * 0.5 * rho * v * v                           # Pa
    return (dp_fric + dp_k).to("Pa")

@dataclass
class StageSolver:
    stage: FirePass | SmokePass | Reversal
    gas: GasStream
    water: WaterStream

    def rhs(self, x, y):
        Tg, pg, hw = y

        # sync local state for all property calls
        self.gas.temperature = Q_(Tg, "K")
        self.gas.pressure    = Q_(pg, "Pa")
        self.water.enthalpy  = Q_(hw, "J/kg")

        qph, _, _ = _overall_Q_for_T(
            Tg, pg, self.water.pressure.to("Pa").magnitude, hw,
            self.stage, self.gas, self.water
        )

        dTgdx = (- qph / (self.gas.mass_flow_rate * self.gas.specific_heat)).to("K/m")
        dpdx  = (- _dp_gas(self.stage, self.gas) / self.stage.hot_side.inner_length).to("Pa/m")
        dhwdx = (  qph / self.water.mass_flow_rate).to("J/(kg*m)")
        return [dTgdx.magnitude, dpdx.magnitude, dhwdx.magnitude]

    def solve(self):
        L = self.stage.hot_side.inner_length.to("m").magnitude
        y0 = [
            self.gas.temperature.to("K").magnitude,
            self.gas.pressure.to("Pa").magnitude,
            self.water.enthalpy.to("J/kg").magnitude,
        ]
        sol = solve_ivp(self.rhs, (0.0, L), y0=y0, dense_output=False)
        return sol
