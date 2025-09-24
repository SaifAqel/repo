from dataclasses import dataclass

@dataclass
class WaterProps:
    T_sat: float   # Saturation temperature [K]
    rho_l: float   # Liquid density [kg/m^3]
    rho_v: float   # Vapor density [kg/m^3]
    mu_l: float    # Liquid dynamic viscosity [Pa·s]
    cp_l: float    # Liquid specific heat [J/kg·K]
    h_fg: float    # Latent heat of vaporization [J/kg]
    sigma: float   # Surface tension [N/m]
    Pr_l: float    # Liquid Prandtl number [-]
