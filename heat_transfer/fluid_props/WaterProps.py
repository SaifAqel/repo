from iapws import IAPWS97
from common.units import ureg, Q_


class WaterProps:
    def Tsat(self, P: Q_) -> Q_:
        w = IAPWS97(P=P.to("megapascal").magnitude, x=0)
        return Q_(w.T, "kelvin")

    def rho_l(self, T: Q_) -> Q_:
        w = IAPWS97(T=T.to("kelvin").magnitude, x=0)
        return Q_(w.rho, "kg/m^3")

    def rho_v(self, T: Q_) -> Q_:
        w = IAPWS97(T=T.to("kelvin").magnitude, x=1)
        return Q_(w.rho, "kg/m^3")

    def mu_l(self, T: Q_) -> Q_:
        w = IAPWS97(T=T.to("kelvin").magnitude, x=0)
        return Q_(w.my, "Pa*s")

    def k_l(self, T: Q_) -> Q_:
        w = IAPWS97(T=T.to("kelvin").magnitude, x=0)
        return Q_(w.k, "W/(m*K)")

    def cp_l(self, T: Q_) -> Q_:
        w = IAPWS97(T=T.to("kelvin").magnitude, x=0)
        return Q_(w.cp, "kJ/(kg*K)")

    def sigma(self, T: Q_) -> Q_:
        w = IAPWS97(T=T.to("kelvin").magnitude, x=0)
        return Q_(w.sigma, "N/m")

    def h_fg(self, P: Q_) -> Q_:
        w_l = IAPWS97(P=P.to("megapascal").magnitude, x=0)
        w_v = IAPWS97(P=P.to("megapascal").magnitude, x=1)
        return Q_(w_v.h - w_l.h, "kJ/kg")

    def h_l(self, T: Q_) -> Q_:
        w = IAPWS97(T=T.to("kelvin").magnitude, x=0)
        return Q_(w.h, "kJ/kg")

    def h_v(self, T: Q_) -> Q_:
        w = IAPWS97(T=T.to("kelvin").magnitude, x=1)
        return Q_(w.h, "kJ/kg")
