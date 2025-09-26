# ---- 3) SystemRunner: chain passes and handle outer iteration if needed ----
class SystemRunner:
    def __init__(self, passes):  # list of PassRunner objects in flow order
        self.passes = passes

    def once(self, gas_in: Stream, water_in: Stream, N_each: int):
        gas = gas_in
        water = water_in
        Q_sum = 0.0
        for p in self.passes:
            gas, water, Q = p.run(gas, water, N_each)
            Q_sum += Q
        return gas, water, Q_sum

    def iterate(self, gas_in, water_guess, N_each, alpha=0.4, tol=1e-3, itmax=50):
        water_in = water_guess
        for _ in range(itmax):
            gas_out, water_out, _ = self.once(gas_in, water_in, N_each)
            new_T = alpha*water_out.T + (1-alpha)*water_in.T
            if abs(new_T - water_in.T) < tol:
                return gas_out, Stream(water_in.m_dot, new_T, water_in.p), water_out
            water_in = Stream(water_in.m_dot, new_T, water_in.p)
        return gas_out, water_in, water_out  # last iterate