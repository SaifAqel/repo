import cantera as ct
from common.units import Q_
from iapws import IAPWS97
from typing import TYPE_CHECKING

class GasProps:
   
    _base_sol = None

    @staticmethod
    def _set_state(gas):
        if GasProps._base_sol is None:
            GasProps._base_sol = ct.Solution("heat_transfer/config/flue_cantera.yaml")
        sol = GasProps._base_sol
        T = gas.temperature.to("K").magnitude 
        P = gas.pressure.to("Pa").magnitude
        X = {k: v.magnitude for k, v in gas.composition.items()}
        sol.TPX = T, P, X
        return sol

    @staticmethod
    def thermal_conductivity(gas):
        sol = GasProps._set_state(gas)
        return Q_(sol.thermal_conductivity, "W/(m*K)")
    
    @staticmethod
    def viscosity(gas):
        sol = GasProps._set_state(gas)
        return Q_(sol.viscosity, "Pa*s")
    
    @staticmethod
    def density(gas):
        sol = GasProps._set_state(gas)
        return Q_(sol.density, "kg/m^3")
    
    @staticmethod
    def enthalpy(gas):
        sol = GasProps._set_state(gas)
        return Q_(sol.enthalpy_mass, "J/kg")
    
    @staticmethod
    def specific_heat(gas):
        sol = GasProps._set_state(gas)
        return Q_(sol.cp_mass, "J/(kg*K)")

class WaterProps:

    @staticmethod
    def quality(water) -> Q_:
        P = water.pressure.to("megapascal").magnitude
        h = water.enthalpy.to("kJ/kg").magnitude
        h_l = IAPWS97(P=P, x=0).h
        h_v = IAPWS97(P=P, x=1).h
        x = (h - h_l) / (h_v - h_l)
        return Q_(x, "dimensionless")

    @staticmethod
    def _state(water) -> IAPWS97:
        P = water.pressure.to("megapascal").magnitude
        h = water.enthalpy.to("kJ/kg").magnitude
        x = water.quality
        if 0 <= x <= 1:
            return IAPWS97(P=P, x=x)
        else:
            return IAPWS97(P=P, h=h)

    @staticmethod
    def temperature(water) -> Q_:
        return Q_(WaterProps._state(water).T, "kelvin")

    @staticmethod
    def density(water) -> Q_:
        P = water.pressure.to("megapascal").magnitude
        x = water.quality
        if 0.0 < x < 1.0:
            rho_l = IAPWS97(P=P, x=0).rho
            rho_v = IAPWS97(P=P, x=1).rho
            rho_mix = 1.0 / ((x / rho_v) + ((1.0 - x) / rho_l))
            return Q_(rho_mix, "kg/m^3")
        return Q_(WaterProps._state(water).rho, "kg/m^3")
 
    @staticmethod
    def dynamic_viscosity(water) -> Q_:
        x = water.quality
        if 0.0 < x < 1.0:
            raise ValueError("Viscosity undefined for two-phase mixture without a mixing rule.")
        return Q_(WaterProps._state(water).mu, "kg/m/s")

    @staticmethod
    def thermal_conductivity(water) -> Q_:
        x = water.quality
        if 0.0 < x < 1.0:
            raise ValueError("Thermal conductivity undefined for two-phase mixture without a mixing rule.")
        return Q_(WaterProps._state(water).k, "W/(m*K)")

    @staticmethod
    def specific_heat(water) -> Q_:
        x = water.quality
        if 0.0 < x < 1.0:
            raise ValueError("cp undefined for two-phase mixture without a mixing rule.")
        return Q_(WaterProps._state(water).cp * 1e3, "J/(kg*K)")

    @staticmethod
    def surface_tension(water) -> Q_:
        P = water.pressure.to("megapascal").magnitude
        return Q_(IAPWS97(P=P, x=0).sigma, "N/m")

    @staticmethod
    def specific_enthalpy(water) -> Q_:
        return Q_(WaterProps._state(water).h * 1e3, "J/kg")

    @staticmethod
    def saturation_temperature(water) -> Q_:
        P = water.pressure.to("megapascal").magnitude
        return Q_(IAPWS97(P=P, x=0).T, "kelvin")

    @staticmethod
    def latent_heat_of_vaporization(water) -> Q_:
        P = water.pressure.to("megapascal").magnitude
        wl = IAPWS97(P=P, x=0.0)
        wg = IAPWS97(P=P, x=1.0)
        return Q_((wg.h - wl.h) * 1e3, "J/kg")





