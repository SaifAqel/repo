class HeatLossCalculator:
    def heat_loss(self, heat_rate, velocity, distance):
        return heat_rate * (distance / velocity)

class GasTemperature:
    def new_temperature(self, T1, heat_loss, mass, h, inv_h):
        return inv_h(h(T1) + heat_loss / mass)
