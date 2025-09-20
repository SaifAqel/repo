
from thermo.config.schemas import load_settings
from thermo.core.thermo_provider import IThermoProvider
from thermo.core.coolprop_provider import CoolPropThermoProvider
from thermo.core.heat_capacity import MixtureCp
from thermo.core.composition import Composition, mix_molar_mass
from thermo.core.streams import GasStream
from thermo.core.heats import compute_LHV_HHV
from thermo.core.stoch import stoich_O2_required_per_mol_fuel, air_flow_rates
from thermo.core.flue import from_fuel_and_air
from thermo.core.balances import sensible_heat, total_input_heat
from thermo.core.root_solvers import solve_brentq
from thermo.services.adiabatic_flame_temperature import AdiabaticFlameTemperature
from thermo.services.combustor import Combustor
from thermo.core.units import ureg, Q_
import yaml
import math
from pathlib import Path

sigma = 5.670374419e-8 * ureg.watt/(ureg.meter**2*ureg.kelvin**4)  # Stefan-Boltzmann

def load_passes_config(path: str):
    with open(path, "r") as f:
        return yaml.safe_load(f)

def rad_conv_step(T_gas_K: Q_, m_dot_gas: Q_, cp_mix: Q_, Tw_K: Q_, seg: dict):
    """
    Single segment heat removal with radiation + convection.
    seg dict keys:
      - A_rad_m2, eps_gas, A_conv_m2, h_conv_W_m2K
      - optional: name
    Returns: new_T_gas_K, q_rad_W, q_conv_W, q_total_W
    """
    A_rad = (seg.get("A_rad_m2", 0.0) or 0.0) * ureg.meter**2
    eps = seg.get("eps_gas", 0.0) or 0.0
    A_conv = (seg.get("A_conv_m2", 0.0) or 0.0) * ureg.meter**2
    h = (seg.get("h_conv_W_m2K", 0.0) or 0.0) * ureg.watt/(ureg.meter**2*ureg.kelvin)

    # Heat rates at current Tg
    q_rad = eps * sigma * A_rad * (T_gas_K**4 - Tw_K**4)
    q_conv = h * A_conv * (T_gas_K - Tw_K)
    q_total = (q_rad + q_conv).to(ureg.watt)

    # Temperature drop
    dT = (q_total / (m_dot_gas * cp_mix)).to(ureg.kelvin)
    T_out = T_gas_K - dT
    return T_out, q_rad.to(ureg.watt), q_conv.to(ureg.watt), q_total

def main():
    # Inputs
    s = load_settings("thermo/config/settings.toml")
    cfg = load_passes_config("passes.yaml")

    thermo: IThermoProvider = CoolPropThermoProvider(s.species_cp_fluids_map)
    cp = MixtureCp(s.species_cp_fluids_map)
    aft = AdiabaticFlameTemperature(cp, solve_brentq)

    # Fuel and air streams like run_example.py
    M = s.species_molar_masses
    air = GasStream(
        T=s.air_T.to("K"),
        P=s.air_P,
        composition=Composition(s.air_composition_mol, "mole")
    )

    fuel = GasStream(
        T=s.fuel_T.to("K"),
        P=s.fuel_P,
        composition=Composition(s.fuel_composition_mass, "mass"),
        flow_mass=s.fuel_mass_flow
    )

    svc = Combustor(s, thermo, cp,
                    lambda x,M,m,dH,Lat,Mw: compute_LHV_HHV(x,M,m,dH,Lat,Mw),
                    (stoich_O2_required_per_mol_fuel, air_flow_rates),
                    from_fuel_and_air,
                    (sensible_heat, total_input_heat),
                    solve_brentq,
                    aft)

    from thermo.models.combustion_case import CombustionCase
    res = svc.run(CombustionCase(air=air, fuel=fuel, excess_air_ratio=s.excess_air_ratio, T_ref=s.T_ref))

    # Initial gas state after combustor
    T_g = res.T_ad_K  # adiabatic flame temperature
    m_dot = res.flue_mass_flow_kg_s
    w_fracs = Composition(res.flue_x, "mole").to_mass(M).fractions
    cp_mix = cp.cp_mass_mixture(T_g, air.P, w_fracs)

    Tw = Q_(cfg["boiling_water_T_K"], ureg.kelvin)

    segments = []
    # Build ordered list of segments: pass1, rev1, pass2, rev2, pass3
    for key in ["pass1", "rev1", "pass2", "rev2", "pass3"]:
        seg = cfg.get(key, {})
        seg.setdefault("name", key)
        segments.append(seg)

    timeline = []
    for seg in segments:
        T_out, q_rad, q_conv, q_total = rad_conv_step(T_g, m_dot, cp_mix, Tw, seg)
        timeline.append({
            "segment": seg.get("name"),
            "T_in_K": T_g.to("kelvin").m,
            "T_out_K": T_out.to("kelvin").m,
            "q_rad_W": q_rad.to("watt").m,
            "q_conv_W": q_conv.to("watt").m,
            "q_total_W": q_total.to("watt").m
        })
        T_g = T_out  # update
        # update cp at new temperature
        cp_mix = cp.cp_mass_mixture(T_g, air.P, w_fracs)

    print("Adiabatic flame temperature [K]:", res.T_ad_K.to("kelvin").m)
    print("Flue gas mass flow [kg/s]:", res.flue_mass_flow_kg_s.to("kg/s").m)
    print("Water shell temperature [K]:", Tw.to("kelvin").m)
    print("---- Passes summary ----")
    for row in timeline:
        print(f"{row['segment']}: T_in={row['T_in_K']:.1f}K -> T_out={row['T_out_K']:.1f}K  "
              f"q_rad={row['q_rad_W']:.2e} W  q_conv={row['q_conv_W']:.2e} W  q_tot={row['q_total_W']:.2e} W")
    print("Outlet flue gas temperature [K]:", timeline[-1]["T_out_K"] if timeline else res.T_ad_K.to("kelvin").m)

if __name__ == "__main__":
    main()
