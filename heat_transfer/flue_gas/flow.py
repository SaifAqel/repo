class GasVelocityCalculator:
    def __init__(self, mass_flow, density, area):
        self.mass_flow = mass_flow      # kg/s
        self.density = density          # kg/m^3
        self.area = area                # m^2

    def velocity(self):
        return self.mass_flow / (self.density * self.area)