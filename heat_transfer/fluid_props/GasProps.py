import cantera as ct
from common.units import ureg, Q_

class HotFlueGas:
    def __init__(self, gas: ct.Solution):
        self.gas = gas

    def thermal_conductivity(self, T, P, X):
        self.gas.TPX = Q_(T).to("K").m, Q_(P).to("Pa").m, X
        return Q_(self.gas.thermal_conductivity, "W/(m*K)")

    def viscosity(self, T, P, X):
        self.gas.TPX = Q_(T).to("K").m, Q_(P).to("Pa").m, X
        return Q_(self.gas.viscosity, "Pa*s")

    def density(self, T, P, X):
        self.gas.TPX = Q_(T).to("K").m, Q_(P).to("Pa").m, X
        return Q_(self.gas.density, "kg/m^3")

    def enthalpy(self, T, P, X):
        self.gas.TPX = Q_(T).to("K").m, Q_(P).to("Pa").m, X
        return Q_(self.gas.enthalpy_mass, "J/kg")

    def cp(self, T, P, X):
        self.gas.TPX = Q_(T).to("K").m, Q_(P).to("Pa").m, X
        return Q_(self.gas.cp_mass, "J/(kg*K)")