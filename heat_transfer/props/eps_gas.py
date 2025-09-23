import numpy as np

class WSGG:
    def __init__(self, P_Pa, a, b, active_species):
        self.P_Pa = P_Pa
        self.a = np.array(a, dtype=float)
        self.b = np.array(b, dtype=float)
        self.species = tuple(active_species)

    def u_atm_m(self, y, L_m):
        y_participating = sum(y.get(sp, 0.0) for sp in self.species)
        return (self.P_Pa / ATM_PA) * y_participating * L_m

    def epsilon(self, u):
        return 1.0 - np.sum(self.a * np.exp(-self.b * u))

    def emissivity(self, y, L_m):
        u = self.u_atm_m(y, L_m)
        return self.epsilon(u)
    