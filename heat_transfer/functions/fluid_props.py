from common.units import Q_, Converter
import cantera as ct
from iapws import IAPWS97

class GasProps:
    ######################### Function #########################
    @staticmethod
    def _set_state(gas, film_temperature: Q_ | None):
        sol = ct.Solution("heat_transfer/config/flue_cantera.yaml")
        T = (film_temperature or gas.temperature).to("K").magnitude
        P = gas.pressure.to("Pa").magnitude
        X = {k: v.magnitude for k, v in gas.composition.items()}
        sol.TPX = T, P, X
        return sol
    ######################### Properties #########################
    @staticmethod
    def thermal_conductivity(gas, film_temperature: Q_ | None = None):
        sol = GasProps._set_state(gas, film_temperature)
        return Q_(sol.thermal_conductivity, "W/(m*K)")
    
    @staticmethod
    def viscosity(gas, film_temperature: Q_ | None = None):
        sol = GasProps._set_state(gas, film_temperature)
        return Q_(sol.viscosity, "Pa*s")
    
    @staticmethod
    def density(gas, film_temperature: Q_ | None = None):
        sol = GasProps._set_state(gas, film_temperature)
        return Q_(sol.density, "kg/m^3")
    
    @staticmethod
    def enthalpy(gas, film_temperature: Q_ | None = None):
        sol = GasProps._set_state(gas, film_temperature)
        return Q_(sol.enthalpy_mass, "J/kg")
    
    @staticmethod
    def specific_heat(gas, film_temperature: Q_ | None = None):
        sol = GasProps._set_state(gas, film_temperature)
        return Q_(sol.cp_mass, "J/(kg*K)")


class WaterProps:
    ######################### Functions ######################### 
    @staticmethod
    def sat_liq(water) -> IAPWS97:
        P = Converter._MPa(water.pressure).magnitude
        return IAPWS97(P=P, x=0.0)

    @staticmethod
    def sat_vap(water) -> IAPWS97:
        P = Converter._MPa(water.pressure).magnitude
        return IAPWS97(P=P, x=1.0)

    @staticmethod
    def _state(water) -> IAPWS97:
        P = Converter._MPa(water.pressure).magnitude
        T = getattr(water, "temperature", None)
        h = getattr(water, "enthalpy", None)
        x = getattr(water, "quality", None)

        if x is not None:
            return IAPWS97(P=P, x=Converter._dim(x).magnitude)
        if T is not None:
            return IAPWS97(P=P, T=Converter._K(T).magnitude)
        if h is not None:
            return IAPWS97(P=P, h=Converter._kJkg(h).magnitude)
        raise ValueError("Provide one of: temperature, enthalpy, or quality.")

    ######################### Saturation Properties #########################
    @staticmethod
    def saturation_temperature(water) -> Q_:
        return Q_(WaterProps.sat_liq(water).T, "K")

    @staticmethod
    def saturation_enthalpy_liquid(water) -> Q_:
        return Q_(WaterProps.sat_liq(water).h, "kJ/kg")

    @staticmethod
    def latent_heat(water) -> Q_:
        s_l = WaterProps.sat_liq(water)
        s_v = WaterProps.sat_vap(water)
        return Q_(s_v.h - s_l.h, "kJ/kg")

    @staticmethod
    def surface_tension(water) -> Q_:
        s = WaterProps.sat_liq(water)
        return Q_(s.sigma, "N/m")

    ######################### Properties at specified state #########################
    @staticmethod
    def temperature(water) -> Q_:
        return Q_(WaterProps._state(water).T, "K")

    @staticmethod
    def density(water) -> Q_:
        return Q_(WaterProps._state(water).rho, "kg/m^3")

    @staticmethod
    def dynamic_viscosity(water) -> Q_:
        return Q_(WaterProps._state(water).mu, "Pa*s")

    @staticmethod
    def thermal_conductivity(water) -> Q_:
        return Q_(WaterProps._state(water).k, "W/m/K")

    @staticmethod
    def specific_heat_cp(water) -> Q_:
        return Q_(WaterProps._state(water).cp, "kJ/kg/K")
    
    @staticmethod
    def enthalpy(water) -> Q_:
        return Q_(WaterProps._state(water).h, "kJ/kg")


    ######################### Utilities #########################
    @staticmethod
    def quality_from_h(water) -> Q_:
        h = Converter._kJkg(water.enthalpy).magnitude
        h_f = WaterProps.sat_liq(water).h
        h_g = WaterProps.sat_vap(water).h
        x = (h - h_f) / (h_g - h_f)
        if 0 < x < 1:
            return Q_(x, "dimensionless")
        return None