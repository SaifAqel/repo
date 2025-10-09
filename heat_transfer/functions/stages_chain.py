# heat_transfer/runner/min_counterflow_chain.py

from __future__ import annotations
from dataclasses import dataclass
from common.units import Q_
from heat_transfer.config.models import GasStream, WaterStream
from heat_transfer.functions.stage_solver import HeatStageSolver


@dataclass
class MinCounterflowChain:

    stages: list
    gas_in: GasStream
    water_in: WaterStream

    def _terminal_hw(self, hw0: Q_) -> Q_:
        TgQ = self.gas_in.temperature.to("K")
        pgQ = self.gas_in.pressure.to("Pa")
        hwQ = hw0.to("J/kg")
        hw_L_inlet = self.water_in.enthalpy.to("J/kg")

        for geometry in self.stages:
            g = self.gas_in.with_state(temperature=TgQ, pressure=pgQ)
            w = self.water_in.with_state(enthalpy=hwQ)
            sol = HeatStageSolver(geometry=geometry, gas=g, water=w).solve()
            TgQ = Q_(sol.y[0, -1], "K")
            pgQ = Q_(sol.y[1, -1], "Pa")
            hwQ = Q_(sol.y[2, -1], "J/kg")

        return hwQ  # hw at x = L

    def run(self, tol: Q_ = Q_(50.0, "J/kg"), max_iter: int = 10):
        
        hw_L_inlet = self.water_in.enthalpy.to("J/kg")

        h0a = hw_L_inlet.to("J/kg")
        h0b = Q_(h0a.magnitude * 0.99, "J/kg")

        hwLa = self._terminal_hw(h0a)
        fa = (hwLa - hw_L_inlet).to("J/kg").magnitude

        hwLb = self._terminal_hw(h0b)
        fb = (hwLb - hw_L_inlet).to("J/kg").magnitude

        best_res, best_h0 = abs(fa), h0a
        if abs(fb) < best_res:
            best_res, best_h0 = abs(fb), h0b

        for k in range(max_iter):
            if abs(fb) <= tol.to("J/kg").magnitude:
                return {
                    "converged": True,
                    "iterations": k + 2,
                    "hw0": h0b,
                    "residual": Q_(fb, "J/kg"),
                }
            if fb == fa:
                break

            h0c_val = h0b.to("J/kg").magnitude - fb * (
                h0b.to("J/kg").magnitude - h0a.to("J/kg").magnitude
            ) / (fb - fa)
            h0c = Q_(h0c_val, "J/kg")

            hwLc = self._terminal_hw(h0c)
            fc = (hwLc - hw_L_inlet).to("J/kg").magnitude

            h0a, fa = h0b, fb
            h0b, fb = h0c, fc

            if abs(fb) < best_res:
                best_res, best_h0 = abs(fb), h0b

        return {
            "converged": False,
            "iterations": max_iter,
            "hw0": best_h0,
            "residual": Q_(best_res, "J/kg"),
        }
