from iapws import IAPWS97

_MPA = 1e6
_KJ = 1e3

class WaterProps:
    def Tsat(self, P: float) -> float:
        return IAPWS97(P=P/_MPA, x=0).T

    def rho_l(self, T: float) -> float:
        return IAPWS97(T=T, x=0).rho

    def rho_v(self, T: float) -> float:
        return IAPWS97(T=T, x=1).rho

    def mu_l(self, T: float) -> float:
        return IAPWS97(T=T, x=0).my

    def k_l(self, T: float) -> float:
        return IAPWS97(T=T, x=0).k

    def cp_l(self, T: float) -> float:
        return IAPWS97(T=T, x=0).cp * _KJ

    def sigma(self, T: float) -> float:
        return IAPWS97(T=T, x=0).sigma

    def h_fg(self, P: float) -> float:
        w_l = IAPWS97(P=P/_MPA, x=0)
        w_v = IAPWS97(P=P/_MPA, x=1)
        return (w_v.h - w_l.h) * _KJ

    def h_l(self, T: float) -> float:
        return IAPWS97(T=T, x=0).h * _KJ

    def h_v(self, T: float) -> float:
        return IAPWS97(T=T, x=1).h * _KJ
