class PrandtlCalculator:
    def __init__(self, mu: float, cp: float, k: float):
        self.mu = mu
        self.cp = cp
        self.k = k

    def compute(self) -> float:
        return self.mu * self.cp / self.k
