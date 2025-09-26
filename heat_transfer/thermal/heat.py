class HeatSystem:
    def __init__(self, PassWithCalc, total_thermal_resistance, GasStream, WaterStream, HotFlueGas):
        self.delta_T = GasStream.inlet_temperature - WaterStream.inlet_temperature
        self.tot_resistance = total_thermal_resistance
        self.PassWithCalc = PassWithCalc
        self.GasStream = GasStream
        self.HotFlueGas = HotFlueGas


    def heat_flux(self):
        return self.delta_T / self.r_total_per_area
    
    def heat_rate(self):
        return self.heat_flux * self.PassWithCalc.segment_heat_transfer_area
    
    def heat_loss(self):
        return self.heat_rate * (self.PassWithCalc.segment_length / self.GasStream.velocity)
    
    def new_gas_temperature(self):
        return self.GasStream.inlet_temperature - (self.heat_loss / (self.GasStream.m_dot * self.HotFlueGas.cp))