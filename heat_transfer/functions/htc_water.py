# heat_transfer/functions/htc_water.py
import math
from common.units import Q_, ureg

class HTCFunctions:
    @staticmethod
    def single_phase(water, D_h: Q_, L: Q_, T_wall: Q_) -> Q_:
        # bulk props
        rho = water.density.to("kg/m^3").magnitude
        mu  = water.dynamic_viscosity.to("Pa*s").magnitude
        k   = water.thermal_conductivity.to("W/(m*K)").magnitude
        cp  = water.specific_heat.to("J/(kg*K)").magnitude
        v   = water.velocity.to("m/s").magnitude
        Dh  = D_h.to("m").magnitude
        Lm  = L.to("m").magnitude

        # wall viscosity (needed for Sieder–Tate)
        Tsat = water.saturation_temperature
        if T_wall < Tsat:
            mu_w = water.water_props.mu_l(water.pressure, T_wall).to("Pa*s").magnitude
        else:
            mu_w = water.water_props.mu_v(water.pressure, T_wall).to("Pa*s").magnitude

        Re = rho*v*Dh/mu
        Pr = mu*cp/k

        if Re < 2300:
            Nu = 1.86*((Re*Pr*Dh/Lm)**(1/3)) * (mu/mu_w)**0.14
        else:
            # Gnielinski with Petukhov friction factor
            f = (0.79*math.log10(Re) - 1.64)**(-2)
            Nu = (f/8)*(Re-1000)*Pr / (1 + 12.7*((f/8)**0.5)*(Pr**(2/3) - 1))

        h = Nu*k/Dh
        return Q_(h, "W/(m^2*K)")
