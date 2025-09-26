import cantera as ct

class HotFlueGas:
    def __init__(self, gas: ct.Solution):
        self.gas = gas

    def thermal_conductivity(self, T, P, X):
        self.gas.TPX = T, P, X
        return self.gas.thermal_conductivity

    def viscosity(self, T, P, X):
        self.gas.TPX = T, P, X
        return self.gas.viscosity

    def density(self, T, P, X):
        self.gas.TPX = T, P, X
        return self.gas.density

    def enthalpy(self, T, P, X):
        self.gas.TPX = T, P, X