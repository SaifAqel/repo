from typing import Dict
from thermo.core.species import parse_CH
from .units import ureg, Q_

def compute_LHV_HHV(
    mole_fractions: Dict[str, float],   # dimensionless
    M_mix,                              # kg/mol  (Quantity)
    mass_flow,                          # kg/s    (Quantity)
    dHf: Dict[str, any],                # kJ/mol  (Quantities)
    latent_H2O,                         # kJ/kg   (Quantity)
    M_H2O                               # kg/mol  (Quantity)
):
    # product enthalpies (kJ/mol)
    H2O_liq = dHf["H2O"]                                # kJ/mol
    H2O_vap = dHf["H2O"] + latent_H2O * M_H2O          # (kJ/kg)*(kg/mol) = kJ/mol

    # initialize as quantities, not floats
    react = 0 * dHf["CO2"]                              # kJ/mol
    HHV_p = 0 * dHf["CO2"]
    LHV_p = 0 * dHf["CO2"]

    for comp, x in mole_fractions.items():
        # use quantity zero as default to avoid dimensionless 0.0
        dh = dHf.get(comp, 0 * dHf["CO2"])              # kJ/mol
        react += x * dh

        C, H = parse_CH(comp)
        if C is not None:
            HHV_p += x * (C * dHf["CO2"] + (H/2) * H2O_liq)
            LHV_p += x * (C * dHf["CO2"] + (H/2) * H2O_vap)
        elif comp == "S":
            HHV_p += x * dHf["SO2"]
            LHV_p += x * dHf["SO2"]
        else:
            HHV_p += x * dh
            LHV_p += x * dh

    LHV_mol = react - LHV_p                             # kJ/mol
    LHV_kg  = LHV_mol / M_mix                           # (kJ/mol)/(kg/mol)=kJ/kg
    power   = LHV_kg * mass_flow                        # (kJ/kg)*(kg/s)=kJ/s = kW

    return power.to(ureg.kilowatt)
