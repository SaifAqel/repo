from thermo.core.composition import Composition, mix_molar_mass
from thermo.core.streams import GasStream

class Combustor:
    def __init__(self, settings, thermo, cp, hv, st, flue, balances, solver, aft):
        self.s=settings; self.th=thermo; self.cp=cp; self.hv=hv; self.st=st; self.flue=flue
        self.bal=balances; self.solver=solver; self.aft=aft

    def run(self, case):
        M = self.s.species_molar_masses
        air = case.air
        fuel = case.fuel

        fuel_x = fuel.composition.to_mole(M).fractions
        air_x  = air.composition.to_mole(M).fractions

        M_fuel = mix_molar_mass(fuel_x, M)
        M_air  = mix_molar_mass(air_x, M)

        fuel_n = fuel.molar_flow(M)
        power_LHV_kW = self.hv(fuel_x, M_fuel, fuel.mass_flow(M), self.s.formation_enthalpies,
                               self.s.latent_heat_H2O, M["H2O"])

        O2_req = self.st[0](fuel_x, self.s.stoich_O2_per_mol)
        air_n, air_m = self.st[1](air_x, fuel_n, O2_req, case.excess_air_ratio, M_air)

        flue_x, flue_n = self.flue(fuel_n, air_n, fuel_x, air_x, O2_req)

        # cp at inlets
        air_cp = self.cp.cp_mass_mixture(air.T, air.P, air.as_mass_fraction(M))
        fuel_cp = self.cp.cp_mass_mixture(fuel.T, fuel.P, fuel.as_mass_fraction(M))

        fuel_sens = self.bal[0](fuel.mass_flow(M), fuel_cp, fuel.T, case.T_ref)
        air_sens  = self.bal[0](air_m,     air_cp,  air.T,  case.T_ref)
        Q_in      = self.bal[1](fuel_sens, air_sens, power_LHV_kW)

        flue_w = Composition(flue_x, "mole").to_mass(M).fractions
        flue_m = flue_n * mix_molar_mass(flue_x, M)

        T_ad = self.aft.solve(air.P, flue_w, flue_m, Q_in, case.T_ref)

        from thermo.models.results import Results
        return Results(power_LHV_kW, fuel_sens, air_sens, Q_in, air_n, air_m, flue_x, flue_n, flue_m, T_ad)
