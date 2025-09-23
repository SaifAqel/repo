class GasNusselt:
    def __init__(self, Re, Pr, n=0.4):
        self.Re = Re    # Reynolds number
        self.Pr = Pr    # Prandtl number
        self.n = n      # 0.4 for heating, 0.3 for cooling

    def calculate(self):
        return 0.023 * (self.Re ** 0.8) * (self.Pr ** self.n)
