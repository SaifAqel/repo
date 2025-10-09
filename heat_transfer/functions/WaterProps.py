from iapws import IAPWS97
from common.units import ureg, Q_


class WaterProps:

    @staticmethod
    def Tsat(P: Q_) -> Q_:
        w = IAPWS97(P=P.to("megapascal").magnitude, x=0)
        return Q_(w.T, "kelvin")


    # ----- density -----
    @staticmethod
    def rho_l(P: Q_, T: Q_) -> Q_:
    # single-phase liquid at given P and T
        w = IAPWS97(P=P.to("megapascal").magnitude, T=T.to("kelvin").magnitude)
        return Q_(w.rho, "kg/m^3")


    @staticmethod
    def rho_v(P: Q_, T: Q_) -> Q_:
        w = IAPWS97(P=P.to("megapascal").magnitude, T=T.to("kelvin").magnitude)
        return Q_(w.rho, "kg/m^3")


    # Provide saturation fallback helpers
    @staticmethod
    def rho_l_sat(P: Q_) -> Q_:
        w = IAPWS97(P=P.to("megapascal").magnitude, x=0)
        return Q_(w.rho, "kg/m^3")


    @staticmethod
    def rho_v_sat(P: Q_) -> Q_:
        w = IAPWS97(P=P.to("megapascal").magnitude, x=1)
        return Q_(w.rho, "kg/m^3")


    # ----- dynamic viscosity -----
    @staticmethod
    def mu_l(P: Q_, T: Q_) -> Q_:
        w = IAPWS97(P=P.to("megapascal").magnitude, T=T.to("kelvin").magnitude)
        return Q_(w.mu, "Pa*s")


    @staticmethod
    def mu_v(P: Q_, T: Q_) -> Q_:
        w = IAPWS97(P=P.to("megapascal").magnitude, T=T.to("kelvin").magnitude)
        return Q_(w.mu, "Pa*s")


    @staticmethod
    def mu_l_sat(P: Q_) -> Q_:
        w = IAPWS97(P=P.to("megapascal").magnitude, x=0)
        return Q_(w.mu, "Pa*s")


    @staticmethod
    def mu_v_sat(P: Q_) -> Q_:
        w = IAPWS97(P=P.to("megapascal").magnitude, x=1)
        return Q_(w.mu, "Pa*s")


    # ----- thermal conductivity -----
    @staticmethod
    def k_l(P: Q_, T: Q_) -> Q_:
        w = IAPWS97(P=P.to("megapascal").magnitude, T=T.to("kelvin").magnitude)
        return Q_(w.k, "W/(m*K)")


    @staticmethod
    def k_v(P: Q_, T: Q_) -> Q_:
        w = IAPWS97(P=P.to("megapascal").magnitude, T=T.to("kelvin").magnitude)
        return Q_(w.k, "W/(m*K)")


    @staticmethod
    def k_l_sat(P: Q_) -> Q_:
        w = IAPWS97(P=P.to("megapascal").magnitude, x=0)
        return Q_(w.k, "W/(m*K)")


    @staticmethod
    def k_v_sat(P: Q_) -> Q_:
        w = IAPWS97(P=P.to("megapascal").magnitude, x=1)
        return Q_(w.k, "W/(m*K)")


    # ----- specific heat -----
    @staticmethod
    def cp_l(P: Q_, T: Q_) -> Q_:
        w = IAPWS97(P=P.to("megapascal").magnitude, T=T.to("kelvin").magnitude)
        return Q_(w.cp * 1e3, "J/(kg*K)")


    @staticmethod
    def cp_v(P: Q_, T: Q_) -> Q_:
        w = IAPWS97(P=P.to("megapascal").magnitude, T=T.to("kelvin").magnitude)
        return Q_(w.cp * 1e3, "J/(kg*K)")


    @staticmethod
    def cp_l_sat(P: Q_) -> Q_:
        w = IAPWS97(P=P.to("megapascal").magnitude, x=0)
        return Q_(w.cp * 1e3, "J/(kg*K)")


    @staticmethod
    def cp_v_sat(P: Q_) -> Q_:
        w = IAPWS97(P=P.to("megapascal").magnitude, x=1)
        return Q_(w.cp * 1e3, "J/(kg*K)")


    # ----- surface tension -----
    @staticmethod
    def sigma(P: Q_, T: Q_) -> Q_:
    # surface tension typically defined for liquid side. Use single-phase call.
        w = IAPWS97(P=P.to("megapascal").magnitude, T=T.to("kelvin").magnitude)
        return Q_(w.sigma, "N/m")


    @staticmethod
    def sigma_sat(P: Q_) -> Q_:
        w = IAPWS97(P=P.to("megapascal").magnitude, x=0)
        return Q_(w.sigma, "N/m")


    # ----- enthalpies -----
    @staticmethod
    def h_l_sat(P: Q_) -> Q_:
        w = IAPWS97(P=P.to("megapascal").magnitude, x=0)
        return Q_(w.h * 1e3, "J/kg")


    @staticmethod
    def h_v_sat(P: Q_) -> Q_:
        w = IAPWS97(P=P.to("megapascal").magnitude, x=1)
        return Q_(w.h * 1e3, "J/kg")


    @staticmethod
    def h_fg(P: Q_) -> Q_:
        return WaterProps.h_v_sat(P) - WaterProps.h_l_sat(P)


    @staticmethod
    def h_single(P: Q_, T: Q_) -> Q_:
        # single-phase specific enthalpy at P and T
        w = IAPWS97(P=P.to("megapascal").magnitude, T=T.to("kelvin").magnitude)
        return Q_(w.h * 1e3, "J/kg")