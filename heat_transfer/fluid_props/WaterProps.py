from iapws import IAPWS97
from common.units import ureg, Q_


class WaterProps:

    @staticmethod
    def Tsat(P: Q_) -> Q_:
        w = IAPWS97(P=P.to("megapascal").magnitude, x=0)
        return Q_(w.T, "kelvin")

    @staticmethod
    def rho_l(T: Q_) -> Q_:
        w = IAPWS97(T=T.to("kelvin").magnitude, x=0)
        return Q_(w.rho, "kg/m^3")

    @staticmethod
    def rho_v(T: Q_) -> Q_:
        w = IAPWS97(T=T.to("kelvin").magnitude, x=1)
        return Q_(w.rho, "kg/m^3")

    @staticmethod
    def mu_l(T: Q_) -> Q_:
        w = IAPWS97(T=T.to("kelvin").magnitude, x=0)
        return Q_(w.mu, "Pa*s")

    @staticmethod
    def k_l(T: Q_) -> Q_:
        w = IAPWS97(T=T.to("kelvin").magnitude, x=0)
        return Q_(w.k, "W/(m*K)")

    @staticmethod
    def cp_l(T: Q_) -> Q_:
        w = IAPWS97(T=T.to("kelvin").magnitude, x=0)
        return Q_(w.cp * 1e3, "J/(kg*K)")

    @staticmethod
    def sigma(T: Q_) -> Q_:
        w = IAPWS97(T=T.to("kelvin").magnitude, x=0)
        return Q_(w.sigma, "N/m")

    @staticmethod
    def h_fg(P: Q_) -> Q_:
        w_l = IAPWS97(P=P.to("megapascal").magnitude, x=0)
        w_v = IAPWS97(P=P.to("megapascal").magnitude, x=1)
        return Q_((w_v.h - w_l.h) * 1e3, "J/kg")

    @staticmethod
    def h_l(T: Q_) -> Q_:
        w = IAPWS97(T=T.to("kelvin").magnitude, x=0)
        return Q_(w.h * 1e3, "J/kg")

    @staticmethod
    def h_v(T: Q_) -> Q_:
        w = IAPWS97(T=T.to("kelvin").magnitude, x=1)
        return Q_(w.h * 1e3, "J/kg")
