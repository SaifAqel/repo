class Heat_Flux:
    def __init__(self, delta_T, r_total_per_area):
        self.delta_T = delta_T
        self.r_total_per_area = r_total_per_area

    def resistance_per_area(self):
        return (self.delta_T / self.r_total_per_area).to("watt/meter**2")


class Heat_Rate:
    def compute(self, flux, area):
        return (flux * area).to("watt")
    
class HeatLossCalculator:
    def heat_loss(self, heat_rate, velocity, distance):
        return heat_rate * (distance / velocity)

class GasTemperature:
    def new_temperature(self, T1, heat_loss, mass, h, inv_h):
        return inv_h(h(T1) + heat_loss / mass)







    class GasVelocityCalculator:
    def __init__(self, mass_flow, density, area):
        self.mass_flow = mass_flow      # kg/s
        self.density = density          # kg/m^3
        self.area = area                # m^2

    def velocity(self):
        return self.mass_flow / (self.density * self.area)
    

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

            
class PrandtlCalculator:
    def __init__(self, mu: float, cp: float, k: float):
        self.mu = mu
        self.cp = cp
        self.k = k

    def compute(self) -> float:
        return self.mu * self.cp / self.k
        return 0.023 * (self.Re ** 0.8) * (self.Pr ** self.n)