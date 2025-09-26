
    
class GasReynoldsNumber:
    def __init__(self, density, velocity, diameter, viscosity):
        self.density = density
        self.velocity = velocity
        self.diameter = diameter
        self.viscosity = viscosity

    def calculate(self):
        return (self.density * self.velocity * self.diameter) / self.viscosity
    
class GasNusselt:
    def __init__(self, Re, Pr, n=0.4):
        self.Re = Re    # Reynolds number
        self.Pr = Pr    # Prandtl number
        self.n = n      # 0.4 for heating, 0.3 for cooling

    def calculate(self):
        return 0.023 * (self.Re ** 0.8) * (self.Pr ** self.n)
    
class PrandtlCalculator:
    def __init__(self, mu: float, cp: float, k: float):
        self.mu = mu
        self.cp = cp
        self.k = k

    def compute(self) -> float:
        return self.mu * self.cp / self.k