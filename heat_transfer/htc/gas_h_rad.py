import math

class GasRadiationCoefficient:
    def __init__(self):
        pass

    def compute(self, *,
               temperature,
               pressure,
               path_length,
               composition,
               spectroscopic_data,
               geometry,
               **extra):
        kappa = sum(composition[s]*spectroscopic_data[s] for s in composition)
        emissivity = 1.0 - math.exp(-kappa*path_length)
        return {"absorption_coefficient": kappa, "emissivity": emissivity}

    def h(self, *, mean_temperature, emissivity, sigma):
        return 4.0 * sigma * (mean_temperature**3) * emissivity
