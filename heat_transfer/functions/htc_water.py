# heat_transfer/functions/htc_water.py
import math
from common.units import Q_, ureg

class Correlations:
    @staticmethod
    def sieder_tate(Re: float, Pr: float, Dh_m: float, L_m: float, mu: float, mu_w: float, k: float) -> float:
        Nu = 1.86 * ((Re * Pr * Dh_m / L_m) ** (1.0 / 3.0)) * (mu / mu_w) ** 0.14
        return Nu * k / Dh_m

    @staticmethod
    def petukhov_friction_factor(Re: float) -> float:
        return (0.79 * math.log(Re) - 1.64) ** (-2.0)

    @staticmethod
    def gnielinski(Re: float, Pr: float, Dh_m: float, k: float, f: float) -> float:
        Nu = (f / 8.0) * (Re - 1000.0) * Pr / (1.0 + 12.7 * (f / 8.0) ** 0.5 * (Pr ** (2.0 / 3.0) - 1.0))
        return Nu * k / Dh_m

    @staticmethod
    def gungor_winterton(h_lo: float, G: float, qpp: float, x: float,
                          rho_l: float, rho_v: float, mu_l: float, mu_v: float, h_fg: float) -> float:
        Bo  = qpp / (G * h_fg)
        Xtt = ((1.0 - x) / x) ** 0.9 * (rho_v / rho_l) ** 0.5 * (mu_l / mu_v) ** 0.1 if x > 0.0 else float('inf')
        return h_lo * (1.0 + 3000.0 * Bo ** 0.86 + 1.12 * Xtt ** -0.75 if x > 0.0 else 1.0 + 3000.0 * Bo ** 0.86)


class HTCFunctions:
    @staticmethod
    def single_phase(water, D_h: Q_, L: Q_, T_wall: Q_) -> Q_:
        rho = water.density.to("kg/m^3").magnitude
        mu  = water.dynamic_viscosity.to("Pa*s").magnitude
        k   = water.thermal_conductivity.to("W/(m*K)").magnitude
        cp  = water.specific_heat.to("J/(kg*K)").magnitude
        v   = water.velocity.to("m/s").magnitude
        Dh  = D_h.to("m").magnitude
        Lm  = L.to("m").magnitude

        Tsat = water.saturation_temperature
        mu_w = (water.water_props.mu_l(water.pressure, T_wall) if T_wall < Tsat
                else water.water_props.mu_v(water.pressure, T_wall)).to("Pa*s").magnitude

        Re = rho * v * Dh / mu
        Pr = mu * cp / k

        if Re < 2300.0:
            h = Correlations.sieder_tate(Re, Pr, Dh, Lm, mu, mu_w, k)
        else:
            f = Correlations.petukhov_friction_factor(Re)
            h = Correlations.gnielinski(Re, Pr, Dh, k, f)

        return Q_(h, "W/(m^2*K)")

    @staticmethod
    def boiling(water, D_h: Q_, L: Q_, T_wall: Q_, qpp: Q_, x: float = 0.0) -> Q_:
        rho_l = water.density.to("kg/m^3").magnitude
        mu_l  = water.dynamic_viscosity.to("Pa*s").magnitude
        k_l   = water.thermal_conductivity.to("W/(m*K)").magnitude
        cp_l  = water.specific_heat.to("J/(kg*K)").magnitude
        v     = water.velocity.to("m/s").magnitude
        Dh    = D_h.to("m").magnitude
        Lm    = L.to("m").magnitude

        P = water.pressure
        h_fg = water.latent_heat_of_vaporization.to("J/kg").magnitude
        rho_v = water.water_props.rho_v_sat(P).to("kg/m^3").magnitude
        mu_v  = water.water_props.mu_v_sat(P).to("Pa*s").magnitude

        Re = rho_l * v * Dh / mu_l
        Pr = mu_l * cp_l / k_l

        Tsat = water.saturation_temperature
        mu_w = (water.water_props.mu_l(P, T_wall) if T_wall < Tsat
                else water.water_props.mu_v(P, T_wall)).to("Pa*s").magnitude

        if Re < 2300.0:
            h_lo = Correlations.sieder_tate(Re, Pr, Dh, Lm, mu_l, mu_w, k_l)
        else:
            f = Correlations.petukhov_friction_factor(Re)
            h_lo = Correlations.gnielinski(Re, Pr, Dh, k_l, f)

        A = water.stage.cold_side.flow_area.to("m^2").magnitude
        G = water.mass_flow_rate.to("kg/s").magnitude / A
        qpp_val = qpp.to("W/m^2").magnitude

        h_tp = Correlations.gungor_winterton(h_lo, G, qpp_val, x, rho_l, rho_v, mu_l, mu_v, h_fg)
        return Q_(h_tp, "W/(m^2*K)")

    @staticmethod
    def shell(water, D_h: Q_, L: Q_, T_wall: Q_, qpp: Q_, x: float = 0.0) -> Q_:
        h_lo = HTCFunctions.single_phase(water, D_h, L, T_wall).to("W/(m^2*K)").magnitude
        Tb   = water.temperature.to("K").magnitude
        qval = qpp.to("W/m^2").magnitude
        Tw_sp = Tb + qval / h_lo
        Tsat  = water.saturation_temperature.to("K").magnitude

        if Tw_sp <= Tsat:
            return Q_(h_lo, "W/(m^2*K)")
        return HTCFunctions.boiling(water, D_h, L, T_wall, qpp, x)
