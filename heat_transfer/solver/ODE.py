from scipy.integrate import solve_ivp

class FireTubeGasODE:

    def __init__(self, pass_geom, props, heat, friction, m_dot):
        self.g = pass_geom
        self.props = props
        self.heat = heat
        self.friction = friction
        self.m_dot = m_dot

    def rhs(self, z, y):
        T, p = y
        rho = self.props.rho(T, p)
        cp  = self.props.cp(T, p)
        mu  = self.props.mu(T, p)

        G   = self.m_dot / self.g.A                   # kg/(m^2·s)
        Re  = (G * self.g.D_i) / mu                   # ρuD/μ with G=ρu
        fD  = self.friction.f_D(Re, self.g.rel_rough)
        q_  = self.heat.q_prime(z, T, p)              # W/m

        dTdz = - q_ / (self.m_dot * cp)
        dpdz = - (fD / (2.0 * self.g.D_i)) * (G*G) / rho
        return [dTdz, dpdz]

    def solve(self, z_span, y0, method="BDF", **kwargs):
        # z_span = (0.0, self.g.L), y0 = [T_in, p_in]
        return solve_ivp(self.rhs, z_span, y0, method=method, **kwargs)
