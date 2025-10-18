from typing import Callable, Dict, List, Any, Optional
from heat_transfer.functions.heat_rate import HeatRate
from heat_transfer.config.models import FirePass, SmokePass, Reversal, Economiser, GasStream, WaterStream
import math
import copy

class StageSolver:

    def __init__(self, stage: FirePass | SmokePass | Reversal | Economiser, gas: GasStream, water: WaterStream):
        self.stage = stage
        self.gas = gas
        self.water = water
        self.qprime = None

    def iterate_wall_temperature(self, *, guess: Optional[float] = None, rtol: float = 1e-4, atol_T: float = 1e-3, atol_q: float = 1e-3, max_iter: int = 50, omega: float = 0.5,) -> Dict[str, Any]:

        Twi = (
            guess
            or getattr(self.gas, "wall_temperature", None)
            or 0.5 * (self.gas.temperature + self.water.temperature)
        )

        Two = getattr(self.water, "wall_temperature", None)        
        qprime = None

        for k in range(1, max_iter + 1):
            self.gas.wall_temperature = Twi
            self.water.wall_temperature = Two

            qprime_new = HeatRate(self.stage, self.gas, self.water).heat_rate_per_length()
            
            walls = self.gas.update_walls(qprime_new)
            Twi_new = walls["Twi"]
            Two_new = walls.get("Two")

            conv_Twi = abs(Twi_new - Twi) <= max(atol_T, rtol * max(abs(Twi_new), 1.0))
            conv_qprime = (qprime is not None) and (
                abs(qprime_new - qprime) <= max(atol_q, rtol * max(abs(qprime_new), 1.0))
            )
            conv_Two = (Two is not None and Two_new is not None) and (
                abs(Two_new - Two) <= max(atol_T, rtol * max(abs(Two_new), 1.0))
            )            

            if conv_Twi and conv_qprime and (conv_Two or Two_new is None):
                self.gas.wall_temperature = Twi_new
                self.water.wall_temperature = Two_new
                self.qprime = qprime_new
                return {
                    "converged": True,
                    "iterations": k,
                    "Twi": Twi_new,
                    "Two": Two_new,
                    "qprime": qprime_new,
                }
            
            Twi = omega * Twi_new + (1.0 - omega) * Twi
            if Two_new is not None and Two is not None:
                Two = omega * Two_new + (1.0 - omega) * Two
            else:
                Two = Two_new
            qprime = qprime_new

        self.gas.wall_temperature = Twi
        self.water.wall_temperature = Two
        self.qprime = qprime
        return {
            "converged": False,
            "iterations": max_iter,
            "Twi": Twi,
            "Two": Two,
            "qprime": qprime,
        }

    def _rhs(self) -> Dict[str, float]:

        dTgdx = - self.qprime / (self.gas.mass_flow_rate * self.gas.specific_heat)
        dhwdx = + self.qprime / self.water.mass_flow_rate
        dpgdx = - self.gas.friction_factor * self.gas.mass_flow_rate**2 / (2.0 * self.stage.hot_side.hydraulic_diameter * self.stage.hot_side.flow_area**2 * self.gas.density)

        return {"dTgdx": dTgdx, "dhwdx": dhwdx, "dpgdx": dpgdx}

    def solve(self, dx_init: float = 0.1, tol_T: float = 2.0):
        nsteps = int(math.ceil(self.stage.hot_side.inner_length / dx))
        dx = dx_init
        gas_list = []
        water_list = []

        while x < self.stage.hot_side.inner_length:
            res = self.iterate_wall_temperature()
            if not res["converged"]:
                raise RuntimeError("Wall iteration failed")

            gas_list.append(copy.deepcopy(self.gas))
            water_list.append(copy.deepcopy(self.water))

            derivs = self._rhs(qprime=res["qprime"])

            dT_est = abs(derivs["dTgdx"]) * dx
            if dT_est > tol_T:      # too large, cut step
                dx *= 0.5
                continue
            if dT_est < 0.25 * tol_T:  # safe, enlarge step
                dx *= 1.2

            self.gas.temperature += derivs["dTgdx"] * dx
            self.gas.pressure    += derivs["dpgdx"] * dx
            self.water.enthalpy  += derivs["dhwdx"] * dx

            x += dx

        return gas_list, water_list