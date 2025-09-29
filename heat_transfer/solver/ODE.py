from scipy.integrate import solve_ivp

class FireTubeGasODE:
    """
    Minimal ODE for gas-side T(z), p(z) in a fire-tube pass.

    Expects:
      pass_geom: PassWithCalc instance exposing at least:
          .D_i   (m)  inner diameter
          .A     (m^2) flow area
          .P     (m)  inner perimeter
          .L     (m)  length
          .rel_rough  (-) relative roughness (epsilon/D) or equivalent
      props: object with
          .rho(T, p) -> kg/m^3
          .cp(T, p)  -> J/(kg·K)
          .mu(T, p)  -> Pa·s
      heat: object with
          .q_prime(z, T, p) -> W/m   (local heat per unit length, already includes all resistances/radiation/boiling)
      friction: object with
          .f_D(Re, rel_rough) -> (-) Darcy friction factor
      m_dot: float, kg/s (constant)
    """

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
