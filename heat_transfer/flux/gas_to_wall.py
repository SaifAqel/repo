class GasToWallHT:
    def __init__(self, sigma: float):
        self.sigma = sigma

    def qpp_rad(self, T_g: float, T_w: float, eps_g: float, eps_w: float) -> float:
        return self.sigma * (T_g**4 - T_w**4) / (1.0/eps_g + 1.0/eps_w - 1.0)

    def qpp_conv(h: float, T_g: float, T_w: float) -> float:
        return h * (T_g - T_w)
    
    def qpp_total(qpp_rad: float, qpp_conv: float) -> float:
        return qpp_conv + qpp_rad