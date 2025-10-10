from __future__ import annotations
from dataclasses import dataclass
from common.units import Q_
from heat_transfer.config.models import GasStream, WaterStream, Stages
from heat_transfer.functions.stage_solver import HeatStageSolver
from iapws import IAPWS97

@dataclass
class MinCounterflowChain:
    stages: Stages
    gas: GasStream
    water: WaterStream

    def _terminal_hw(self, hw0: Q_) -> Q_:
        TgQ = self.gas.temperature.to("K")
        pgQ = self.gas.pressure.to("Pa")
        hwQ = hw0.to("J/kg")
        PwQ = self.water.pressure

        for stage in self.stages:
            g = GasStream(
                mass_flow_rate=self.gas.mass_flow_rate,
                temperature=TgQ,
                pressure=pgQ,
                composition=self.gas.composition,
                spectroscopic_data=self.gas.spectroscopic_data,
                stage=stage,
                gas_props=self.gas.gas_props,
            )
            w = WaterStream(
                mass_flow_rate=self.water.mass_flow_rate,
                temperature=self.water.temperature,
                pressure=PwQ,
                composition=self.water.composition,
                stage=stage,
                water_props=self.water.water_props,
            )
            sol = HeatStageSolver(stage=stage, gas=g, water=w).solve()
            TgQ = Q_(sol.y[0, -1], "K")
            pgQ = Q_(sol.y[1, -1], "Pa")
            hwQ = Q_(sol.y[2, -1], "J/kg")

        return hwQ

    def _simulate_profile(self, hw0: Q_):
        TgQ = self.gas.temperature.to("K")
        pgQ = self.gas.pressure.to("Pa")
        hwQ = hw0.to("J/kg")
        PwQ = self.water.pressure

        profile = []
        x_offset = 0.0

        for stage in self.stages:
            g0 = GasStream(
                mass_flow_rate=self.gas.mass_flow_rate,
                temperature=TgQ,
                pressure=pgQ,
                composition=self.gas.composition,
                spectroscopic_data=self.gas.spectroscopic_data,
                stage=stage,
                gas_props=self.gas.gas_props,
            )
            w0 = WaterStream(
                mass_flow_rate=self.water.mass_flow_rate,
                temperature=self.water.temperature,
                pressure=PwQ,
                composition=self.water.composition,
                stage=stage,
                water_props=self.water.water_props,
            )

            sol = HeatStageSolver(stage=stage, gas=g0, water=w0).solve()

            for i in range(sol.t.size):
                x_i = x_offset + sol.t[i]
                Tg_i = Q_(sol.y[0, i], "K")
                pg_i = Q_(sol.y[1, i], "Pa")
                hw_i = Q_(sol.y[2, i], "J/kg")

                w_state = IAPWS97(P=PwQ.to("megapascal").magnitude, h=hw_i.to("kJ/kg").magnitude)
                Tw_i = Q_(w_state.T, "K")

                g_i = GasStream(
                    mass_flow_rate=self.gas.mass_flow_rate,
                    temperature=Tg_i,
                    pressure=pg_i,
                    composition=self.gas.composition,
                    spectroscopic_data=self.gas.spectroscopic_data,
                    stage=stage,
                    gas_props=self.gas.gas_props,
                )
                w_i = WaterStream(
                    mass_flow_rate=self.water.mass_flow_rate,
                    temperature=Tw_i,
                    pressure=PwQ,
                    composition=self.water.composition,
                    stage=stage,
                    water_props=self.water.water_props,
                )
                profile.append((Q_(x_i, "m"), g_i, w_i))

            x_offset += sol.t[-1]
            TgQ = Q_(sol.y[0, -1], "K")
            pgQ = Q_(sol.y[1, -1], "Pa")
            hwQ = Q_(sol.y[2, -1], "J/kg")

        return profile

    def run(self, tol: Q_ = Q_(50.0, "J/kg"), max_iter: int = 10):
        hw_L_inlet = self.water._h.to("J/kg")

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
                    "profile": self._simulate_profile(h0b),
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
            "profile": self._simulate_profile(best_h0),
        }
