class GasConvectiveHeatTransfer:
    def __init__(self, k, d, Re, Pr, Nu):
        self.k = k        # Thermal conductivity [W/m-K]
        self.d = d        # Characteristic length/diameter [m]
        self.Nu = Nu      # Nusselt number

    def h(self):
        return self.Nu * self.k / self.d
