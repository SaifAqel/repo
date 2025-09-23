class WallEmissivity:
    def __init__(self, heat_flux, temperature_wall, temperature_surroundings, stefan_boltzmann):
        self.heat_flux = heat_flux
        self.temperature_wall = temperature_wall
        self.temperature_surroundings = temperature_surroundings
        self.sigma = stefan_boltzmann

    def calculate(self):
        return self.heat_flux / (
            self.sigma * (self.temperature_wall**4 - self.temperature_surroundings**4)
        )
