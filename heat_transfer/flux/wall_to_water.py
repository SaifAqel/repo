class PoolBoilingHT:
    """
    Pool boiling on a heated wall.
    Nucleate boiling: Rohsenow correlation.
    Critical heat flux: Zuber correlation.
    """
    def __init__(self, fluid: FluidProps, C_sf: float, n: float):
        """
        fluid: saturated properties at system pressure.
        C_sf: surface-fluid constant (e.g., water–polished copper ~ 0.013).
        n: Rohsenow exponent (water ~ 1.7).
        """
        self.f = fluid
        self.C_sf = C_sf
        self.n = n

    def qpp_rohsenow(self, T_wall: float, T_sat: Optional[float] = None) -> float:
        """
        Heat flux q'' [W/m^2] from Rohsenow.
        Valid for nucleate boiling below CHF.
        """
        T_sat = self.f.T_sat if T_sat is None else T_sat
        dT = max(T_wall - T_sat, 0.0)
        if dT <= 0.0:
            return 0.0
        f = self.f
        g_term = math.sqrt(9.80665 * (f.rho_l - f.rho_v) / f.sigma)
        bracket = (f.cp_l * dT) / (self.C_sf * f.h_fg * (f.Pr_l ** self.n))
        return f.mu_l * f.h_fg * g_term * (bracket ** 3)

    def chf_zuber(self) -> float:
        f = self.f
        return 0.131 * f.h_fg * math.sqrt(f.rho_v) * (f.sigma * 9.80665 * (f.rho_l - f.rho_v))**0.25


    def h_boiling(self, T_wall: float, T_sat: Optional[float] = None) -> float:
        """
        Boiling heat-transfer coefficient h [W/m^2·K] = q'' / (T_wall - T_sat).
        Returns 0 if T_wall <= T_sat.
        """
        T_sat = self.f.T_sat if T_sat is None else T_sat
        dT = T_wall - T_sat
        if dT <= 0.0:
            return 0.0
        return self.qpp_rohsenow(T_wall, T_sat) / dT