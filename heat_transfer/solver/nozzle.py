# heat_transfer/solver/nozzle.py
class Nozzle:
    def __init__(self, K, diameter, m_dot, gas_props):
        self.K = K
        self.D = diameter
        self.A = 0.25 * 3.141592653589793 * (self.D ** 2)
        self.m_dot = m_dot
        self.props = gas_props

    def apply(self, T_in, p_in):
        rho = self.props.rho(T_in, p_in)
        v = self.m_dot / (rho * self.A)
        dp = self.K * 0.5 * rho * v * v
        return T_in, p_in - dp
