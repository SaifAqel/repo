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
from thermo.models.combustion_case import CombustionCase

s = load_settings("thermo/config/settings.toml")
thermo = CoolPropThermoProvider(s.species_cp_fluids_map)
cp = MixtureCp(thermo)
aft = AdiabaticFlameTemperature(cp, solve_brentq)

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

res = svc.run(CombustionCase(air=air, fuel=fuel, excess_air_ratio=s.excess_air_ratio, T_ref=s.T_ref))
print(res)
