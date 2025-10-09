import cantera as ct
from common.units import ureg, Q_
from iapws import IAPWS97

class GasProps:
    def __init__(self, gas: ct.Solution):
        self.gas = gas

    def _set_state(self, T, P, X):
        T_val = Q_(T).to("K").magnitude if not isinstance(T, (int, float)) else T
        P_val = Q_(P).to("Pa").magnitude if not isinstance(P, (int, float)) else P
        self.gas.TPX = T_val, P_val, X

    def thermal_conductivity(self, T, P, X):
        self._set_state(T, P, X)
        return Q_(self.gas.thermal_conductivity, "W/(m*K)")

    def viscosity(self, T, P, X):
        self._set_state(T, P, X)
        return Q_(self.gas.viscosity, "Pa*s")

    def density(self, T, P, X):
        self._set_state(T, P, X)
        return Q_(self.gas.density, "kg/m^3")

    def enthalpy(self, T, P, X):
        self._set_state(T, P, X)
        return Q_(self.gas.enthalpy_mass, "J/kg")

    def cp(self, T, P, X):
        self._set_state(T, P, X)
        return Q_(self.gas.cp_mass, "J/(kg*K)")

class WaterProps:
    @staticmethod
    def Tsat(P: Q_) -> Q_:
        w = IAPWS97(P=P.to("megapascal").magnitude, x=0)
        return Q_(w.T, "kelvin")

    @staticmethod
    def rho_l(P: Q_, T: Q_) -> Q_:
        w = IAPWS97(P=P.to("megapascal").magnitude, T=T.to("kelvin").magnitude)
        return Q_(w.rho, "kg/m^3")

    @staticmethod
    def rho_v(P: Q_, T: Q_) -> Q_:
        w = IAPWS97(P=P.to("megapascal").magnitude, T=T.to("kelvin").magnitude)
        return Q_(w.rho, "kg/m^3")

    @staticmethod
    def rho_l_sat(P: Q_) -> Q_:
        w = IAPWS97(P=P.to("megapascal").magnitude, x=0)
        return Q_(w.rho, "kg/m^3")

    @staticmethod
    def rho_v_sat(P: Q_) -> Q_:
        w = IAPWS97(P=P.to("megapascal").magnitude, x=1)
        return Q_(w.rho, "kg/m^3")

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

    @staticmethod
    def sigma(P: Q_, T: Q_) -> Q_:
        w = IAPWS97(P=P.to("megapascal").magnitude, T=T.to("kelvin").magnitude)
        return Q_(w.sigma, "N/m")

    @staticmethod
    def sigma_sat(P: Q_) -> Q_:
        w = IAPWS97(P=P.to("megapascal").magnitude, x=0)
        return Q_(w.sigma, "N/m")

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