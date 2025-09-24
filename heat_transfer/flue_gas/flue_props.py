import cantera as ct

class HotFlueGasConductivity:
    def __init__(self, gas: "ct.Solution"):
        self.gas = gas
    def k(self, T, P, X):
        self.gas.TPX = T, P, X
        return self.gas.thermal_conductivity

class HotFlueGasViscosity:
    def __init__(self, gas: "ct.Solution"):
        self.gas = gas
    def mu(self, T, P, X):
        self.gas.TPX = T, P, X
        return self.gas.viscosity

class HotFlueGasDensity:
    def __init__(self, gas: "ct.Solution"):
        self.gas = gas
    def rho(self, T, P, X):
        self.gas.TPX = T, P, X
        return self.gas.density
