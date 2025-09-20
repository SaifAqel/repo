def stoich_O2_required_per_mol_fuel(fuel_x, table):
    return sum(fuel_x[k]*table.get(k,0.0) for k in fuel_x)

def air_flow_rates(air_x, fuel_n_dot, O2_req, excess, M_air):
    O2_x = air_x["O2"]
    O2_stoich = fuel_n_dot * O2_req
    O2_actual = O2_stoich * excess
    air_n = O2_actual / O2_x
    air_m = air_n * M_air
    return air_n, air_m
