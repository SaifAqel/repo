class PoolBoilingHTC:
    @staticmethod
    def rohsenow_h(mu_l, h_fg, g, rho_l, rho_v, sigma, cp_l, dT, C_sf, Pr_l, n):
        K = mu_l * h_fg * (g * (rho_l - rho_v) / sigma) ** 0.5 * (cp_l / (C_sf * h_fg * (Pr_l ** n))) ** 3
        return K * (dT ** 2)