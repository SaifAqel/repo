from dataclasses import dataclass, field
from typing import Dict, Optional, Any
from common.units import Q_

@dataclass(frozen=True)
class Results:
    # Primary results (keep names for drop-in use)
    power_LHV_kW: Q_
    fuel_sensible_kW: Q_
    air_sensible_kW: Q_
    Q_in_total_kW: Q_
    air_molar_flow_mol_s: Q_
    air_mass_flow_kg_s: Q_
    flue_x: Dict[str, Any]
    flue_n_dot: Q_
    flue_mass_flow_kg_s: Q_
    T_ad_K: Q_

    # Optional/expanded fields captured across the program
    fuel_x: Dict[str, Any] = field(default_factory=dict)       # mole fractions
    air_x: Dict[str, Any] = field(default_factory=dict)        # mole fractions
    flue_w: Dict[str, Any] = field(default_factory=dict)       # mass fractions
    excess_air_ratio: Optional[float] = None                   # λ
    O2_req_per_mol_fuel: Optional[float] = None                # mol O2 / mol fuel
    air_cp_mass: Optional[Q_] = None                           # J/(kg*K)
    fuel_cp_mass: Optional[Q_] = None                          # J/(kg*K)
    notes: Dict[str, Any] = field(default_factory=dict)        # arbitrary extras

    # ---- formatting controls ----
    _DEC_FLOW = 3           # mol/s, kg/s
    _DEC_KW = 1             # kW
    _DEC_CP = 1             # J/(kg*K)
    _DEC_TEMP = 2           # K
    _DEC_PCT = 4            # percent points

    def __str__(self) -> str:
        out = []
        out.append("=== Combustion Results ===")

        # Air flow
        out.append("\n[Air flow]")
        out.append(f"ṅ_air: {self._fmt(self.air_molar_flow_mol_s, 'mol/s', self._DEC_FLOW)}")
        out.append(f"ṁ_air: {self._fmt(self.air_mass_flow_kg_s, 'kg/s', self._DEC_FLOW)}")

        # Energy
        out.append("\n[Energy]")
        out.append(f"Power (LHV):   {self._fmt(self.power_LHV_kW, 'kW', self._DEC_KW)}")
        out.append(f"Sensible fuel:  {self._fmt(self.fuel_sensible_kW, 'kW', self._DEC_KW)}")
        out.append(f"Sensible air:   {self._fmt(self.air_sensible_kW, 'kW', self._DEC_KW)}")
        out.append(f"Q_in total:     {self._fmt(self.Q_in_total_kW, 'kW', self._DEC_KW)}")

        # Flue gas
        out.append("\n[Flue gas]")
        out.append(f"ṅ_flue: {self._fmt(self.flue_n_dot, 'mol/s', self._DEC_FLOW)}")
        out.append(f"ṁ_flue: {self._fmt(self.flue_mass_flow_kg_s, 'kg/s', self._DEC_FLOW)}")

        if self.flue_x:
            out.append("\nFlue composition, mole %")
            out += self._fmt_frac_dict(self.flue_x, percent_dec=self._DEC_PCT)
        if self.flue_w:
            out.append("\nFlue composition, mass %")
            out += self._fmt_frac_dict(self.flue_w, percent_dec=self._DEC_PCT)

        # Inputs and props (optional)
        opt_lines = []
        if self.excess_air_ratio is not None:
            opt_lines.append(f"Excess air ratio: {self._fmt_float(self.excess_air_ratio, 4)}")
        if self.O2_req_per_mol_fuel is not None:
            opt_lines.append(f"O2 required: {self._fmt_float(self.O2_req_per_mol_fuel, 4)} mol O2/mol fuel")
        if self.air_x:
            opt_lines.append("\n[Air composition, mole %]")
            opt_lines += self._fmt_frac_dict(self.air_x, percent_dec=self._DEC_PCT)
        if self.fuel_x:
            opt_lines.append("\n[Fuel composition, mole %]")
            opt_lines += self._fmt_frac_dict(self.fuel_x, percent_dec=self._DEC_PCT)
        if self.air_cp_mass is not None or self.fuel_cp_mass is not None:
            opt_lines.append("\n[Heat capacities]")
            if self.air_cp_mass is not None:
                opt_lines.append(f"cp_air:  {self._fmt(self.air_cp_mass, 'J/(kg*K)', self._DEC_CP)}")
            if self.fuel_cp_mass is not None:
                opt_lines.append(f"cp_fuel: {self._fmt(self.fuel_cp_mass, 'J/(kg*K)', self._DEC_CP)}")
        if opt_lines:
            out.append("\n[Inputs]")
            out += opt_lines

        # Flame temperature
        out.append("\n[Adiabatic flame]")
        out.append(f"T_ad: {self._fmt(self.T_ad_K, 'K', self._DEC_TEMP)}")

        return "\n".join(out)

    __repr__ = __str__

    # ---- helpers ----
    @staticmethod
    def _mag(x: Any) -> float:
        try:
            return x.m
        except Exception:
            try:
                return float(x)
            except Exception:
                return x  # let formatting fail loudly for unsupported types

    @staticmethod
    def _fmt_float(val: float, dec: int) -> str:
        return f"{val:,.{dec}f}"

    def _fmt(self, q: Q_, unit: str, dec: int) -> str:
        try:
            mag = q.to(unit).m
        except Exception:
            # already plain number in target unit
            mag = self._mag(q)
        return f"{mag:,.{dec}f} {unit}"

    def _fmt_frac_dict(self, d: Dict[str, Any], percent_dec: int) -> list[str]:
        rows = []
        # Accept floats or dimensionless  quantities
        for k, v in sorted(d.items()):
            mag = self._mag(v)
            pct = mag * 100.0
            rows.append(f"  {k:>4}: {pct:8.{percent_dec}f} %")
        return rows

    # Machine-readable form if needed elsewhere
    def to_dict(self) -> Dict[str, Any]:
        return {
            "power_LHV_kW": self._fmt(self.power_LHV_kW, "kW", self._DEC_KW),
            "fuel_sensible_kW": self._fmt(self.fuel_sensible_kW, "kW", self._DEC_KW),
            "air_sensible_kW": self._fmt(self.air_sensible_kW, "kW", self._DEC_KW),
            "Q_in_total_kW": self._fmt(self.Q_in_total_kW, "kW", self._DEC_KW),
            "air_molar_flow_mol_s": self._fmt(self.air_molar_flow_mol_s, "mol/s", self._DEC_FLOW),
            "air_mass_flow_kg_s": self._fmt(self.air_mass_flow_kg_s, "kg/s", self._DEC_FLOW),
            "flue_n_dot": self._fmt(self.flue_n_dot, "mol/s", self._DEC_FLOW),
            "flue_mass_flow_kg_s": self._fmt(self.flue_mass_flow_kg_s, "kg/s", self._DEC_FLOW),
            "T_ad_K": self._fmt(self.T_ad_K, "K", self._DEC_TEMP),
            "flue_x_pct": {k: float(self._mag(v) * 100.0) for k, v in self.flue_x.items()},
            "flue_w_pct": {k: float(self._mag(v) * 100.0) for k, v in self.flue_w.items()} if self.flue_w else {},
            "fuel_x_pct": {k: float(self._mag(v) * 100.0) for k, v in self.fuel_x.items()} if self.fuel_x else {},
            "air_x_pct": {k: float(self._mag(v) * 100.0) for k, v in self.air_x.items()} if self.air_x else {},
            "excess_air_ratio": self.excess_air_ratio,
            "O2_req_per_mol_fuel": self.O2_req_per_mol_fuel,
            "cp_air_J_per_kgK": None if self.air_cp_mass is None else self._mag(self.air_cp_mass.to("J/(kg*K)")),
            "cp_fuel_J_per_kgK": None if self.fuel_cp_mass is None else self._mag(self.fuel_cp_mass.to("J/(kg*K)")),
            "notes": self.notes,
        }
