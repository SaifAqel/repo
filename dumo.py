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