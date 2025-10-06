import cantera as ct
from common.units import ureg, Q_

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