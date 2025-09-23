class GasReynoldsNumber:
    def __init__(self, density, velocity, diameter, viscosity):
        self.density = density
        self.velocity = velocity
        self.diameter = diameter
        self.viscosity = viscosity

    def calculate(self):
        return (self.density * self.velocity * self.diameter) / self.viscosity
