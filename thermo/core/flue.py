def from_fuel_and_air(fuel_n, air_n, fuel_x, air_x, O2_req):
    gf=lambda k:fuel_x.get(k,0.0); ga=lambda k:air_x.get(k,0.0)
    n_CO2 = air_n*ga("CO2")+fuel_n*gf("CO2")+fuel_n*(gf("CH4")+2*gf("C2H6")+3*gf("C3H8")+4*gf("C4H10"))
    n_H2O = fuel_n*gf("H2O")+fuel_n*(2*gf("CH4")+3*gf("C2H6")+4*gf("C3H8")+5*gf("C4H10"))+fuel_n*gf("H2S")+air_n*ga("H2O")
    n_SO2 = fuel_n*gf("H2S")
    n_O2  = air_n*ga("O2") - fuel_n*O2_req
    n_N2  = air_n*ga("N2") + fuel_n*gf("N2")
    n_Ar  = air_n*ga("Ar")
    flows={"CO2":n_CO2,"H2O":n_H2O,"SO2":n_SO2,"O2":n_O2,"N2":n_N2,"Ar":n_Ar}
    n_tot=sum(flows.values())
    x={k:(v/n_tot if n_tot!=0 else 0.0) for k,v in flows.items()}
    return x, n_tot
