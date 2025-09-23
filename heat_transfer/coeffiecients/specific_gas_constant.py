class SpecificGasConstant:
    def __init__(self):
        self.Ru = 8.314462618  # J/(mol·K)

    def calculate(self, M):
        return self.Ru / M
