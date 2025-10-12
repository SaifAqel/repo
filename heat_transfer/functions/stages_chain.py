# heat_transfer/functions/stages_chain.py
from __future__ import annotations
from dataclasses import dataclass
from typing import List, Tuple
from common.units import Q_
from heat_transfer.config.models import Stages, GasStream, WaterStream
from heat_transfer.functions.stage_solver import StageSolver


@dataclass
class MinCounterflowChain:
    stages: Stages
    gas: GasStream
    water: WaterStream

    def run(self, tol: Q_ = Q_(1.0, "J/kg"), max_iter: int = 12):
        stage_list = [self.stages.HX_1, self.stages.HX_2, self.stages.HX_3,
                      self.stages.HX_4, self.stages.HX_5, self.stages.HX_6]

        hw_in = [self.water.enthalpy for _ in stage_list]

        Tg_in = self.gas.temperature
        pg_in = self.gas.pressure

        best = None

        for it in range(max_iter):
            Tg_i, pg_i = Tg_in, pg_in

            profile: List[Tuple[Q_, GasStream, WaterStream]] = []
            x_offset = Q_(0.0, "m")

            cold_out_per_stage: List[Q_] = [None] * 6

            for i, stg in enumerate(stage_list):
                g = GasStream(
                    mass_flow_rate=self.gas.mass_flow_rate,
                    temperature=Tg_i,
                    pressure=pg_i,
                    composition=self.gas.composition,
                    spectroscopic_data=self.gas.spectroscopic_data,
                    stage=stg,
                    gas_props=self.gas.gas_props,
                )
                w = WaterStream(
                    mass_flow_rate=self.water.mass_flow_rate,
                    enthalpy=hw_in[i],
                    pressure=self.water.pressure,
                    composition=self.water.composition,
                    stage=stg,
                    water_props=self.water.water_props,
                )

                sol = StageSolver(stage=stg, gas=g, water=w).solve()

                Tg_i = Q_(sol.y[0, -1], "K")
                pg_i = Q_(sol.y[1, -1], "Pa")
                hw_o = Q_(sol.y[2, -1], "J/kg")
                cold_out_per_stage[i] = hw_o

                for k in range(sol.t.size):
                    x = x_offset + Q_(sol.t[k], "m")
                    gk = GasStream(
                        mass_flow_rate=self.gas.mass_flow_rate,
                        temperature=Q_(sol.y[0, k], "K"),
                        pressure=Q_(sol.y[1, k], "Pa"),
                        composition=self.gas.composition,
                        spectroscopic_data=self.gas.spectroscopic_data,
                        stage=stg,
                        gas_props=self.gas.gas_props,
                    )
                    wk = WaterStream(
                        mass_flow_rate=self.water.mass_flow_rate,
                        enthalpy=Q_(sol.y[2, k], "J/kg"),
                        pressure=self.water.pressure,
                        composition=self.water.composition,
                        stage=stg,
                        water_props=self.water.water_props,
                    )
                    profile.append((x, gk, wk))

                x_offset += Q_(sol.t[-1], "m")

            hw_in_new = hw_in[:]
            hw_in_new[5] = self.water.enthalpy
            hw_in_new[4] = cold_out_per_stage[5]
            hw_in_new[3] = cold_out_per_stage[4]
            hw_in_new[2] = cold_out_per_stage[3]
            hw_in_new[1] = cold_out_per_stage[2]
            hw_in_new[0] = cold_out_per_stage[1]

            res = max(abs((hw_in_new[j] - hw_in[j]).to("J/kg").magnitude) for j in range(6))

            best = {
                "converged": res <= tol.to("J/kg").magnitude,
                "iterations": it + 1,
                "hw0": hw_in_new[5],
                "residual": Q_(res, "J/kg"),
                "profile": profile,
            }

            hw_in = hw_in_new
            if best["converged"]:
                break

            Tg_in, pg_in = Tg_i, pg_i

        return best
