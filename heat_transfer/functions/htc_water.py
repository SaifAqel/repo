# heat_transfer/functions/htc_water.py
import math
from common.units import Q_, ureg
from heat_transfer.config.models import WaterStream

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
        Bo = qpp / max(G * h_fg, 1e-16)  # guard
        if x <= 0.0:
            return h_lo * (1.0 + 3000.0 * Bo**0.86)
        # clamp to avoid singularities near 0 or 1
        x_clamped = min(max(x, 1e-6), 1.0 - 1e-6)
        Xtt = ((1.0 - x_clamped) / x_clamped)**0.9 * (rho_v / rho_l)**0.5 * (mu_l / mu_v)**0.1
        return h_lo * (1.0 + 3000.0 * Bo**0.86 + 1.12 * Xtt**-0.75)


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
        
        # Guard: avoid property calls at nonphysical wall temperatures during bracketing.
        # If T_wall is outside IF97 range, use bulk viscosity for μ_w.
        T_wK = T_wall.to("K").magnitude
        if 273.16 <= T_wK <= 1073.15:
            Tsat = water.saturation_temperature
            mu_w_Q = (water.water_props.mu_l(water.pressure, T_wall)
                      if T_wall <= Tsat else water.water_props.mu_v(water.pressure, T_wall))
            mu_w = mu_w_Q.to("Pa*s").magnitude
        else:
            mu_w = mu  # fallback only affects the (μ/μw)^0.14 correction in Sieder–Tate

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

        # Only use wall μ if T_wall is within IF97 range. Otherwise use bulk μ.
        T_wK = T_wall.to("K").magnitude
        if 273.16 <= T_wK <= 1073.15:
            Tsat = water.saturation_temperature
            mu_w_Q = (water.water_props.mu_l(P, T_wall)
                      if T_wall <= Tsat else water.water_props.mu_v(P, T_wall))
            mu_w = mu_w_Q.to("Pa*s").magnitude
        else:
            mu_w = mu_l

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
    def shell(water:WaterStream, D_h: Q_, L: Q_, T_wall: Q_, qpp: Q_, x: float | None = None) -> Q_:
        # Single-phase baseline
        h_lo = HTCFunctions.single_phase(water, D_h, L, T_wall).to("W/(m^2*K)").magnitude

        # Equilibrium quality from bulk enthalpy
        if x is None:
            P = water.pressure
            h_b = (water.enthalpy).to("J/kg")
            h_l = water.water_props.h_l_sat(P)
            h_v = water.water_props.h_v_sat(P)
            x_eq = ((h_b - h_l) / (h_v - h_l)).to("").magnitude
            x_eq = 0.0 if x_eq < 0.0 else (1.0 if x_eq > 1.0 else x_eq)
        else:
            x_eq = max(0.0, min(1.0, float(x)))

        if x_eq <= 0.0:
            return Q_(h_lo, "W/(m^2*K)")

        # Two-phase HTC with quality
        h_tp = HTCFunctions.boiling(water, D_h, L, T_wall, qpp, x=x_eq).to("W/(m^2*K)").magnitude

        # Smooth blend near onset (optional)
        w = min(1.0, x_eq / 0.02)
        h_eff = (1.0 - w) * h_lo + w * h_tp
        return Q_(h_eff, "W/(m^2*K)")

