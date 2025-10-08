import math

class HTCFunctions:
    def sieder_tate(Re, Pr, D_h, L, mu, mu_w, k):
            Nu = 1.86 * (Re * Pr * D_h / L)**(1/3) * (mu / mu_w)**0.14
            h = Nu * k / D_h
            return h

    def gnielinski(Re, Pr, D_h, L, k, f=None):
            if f is None:
                f = (0.79 * math.log(Re) - 1.64)**(-2)
            Nu = (f/8)*(Re - 1000)*Pr / (1 + 12.7*(f/8)**0.5*(Pr**(2/3) - 1))
            h = Nu * k / D_h
            return h

    def rohsenow(q_flux, delta_T):
            h = q_flux / delta_T
            return h

    def chen(h_conv, h_nb, F, S):
            h = S * h_conv + F * h_nb
            return h

    @staticmethod
    def htc_shell(water_with_calc):
        if water_with_calc.phase == "vapor" or water_with_calc.phase == "liquid":
            if water_with_calc.Re < 2300:
                return HTCFunctions.sieder_tate(Re, Pr, D_h, L, mu, mu_w, k)
            else:
                return HTCFunctions.gnielinski(Re, Pr, D_h, k)
        elif water_with_calc.phase == "saturated":
            if water_with_calc.Re < 5000:
                return HTCFunctions.rohsenow(q_flux, T_wall - T_sat)
            else:
                return HTCFunctions.chen(h_conv, h_nb, S, F)