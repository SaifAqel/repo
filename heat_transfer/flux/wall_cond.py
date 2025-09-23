class WallHT:
    def __init__(self, k: float, thickness: float):
        """
        k: thermal conductivity of wall (W/m-K)
        thickness: wall thickness (m)
        """
        self.k = k
        self.L = thickness

    def qpp_cond(self, T_hot: float, T_cold: float) -> float:
        """Conductive heat flux (W/m^2) across wall."""
        return self.k * (T_hot - T_cold) / self.L