import math

class DarcyFriction:
    @staticmethod
    def colebrook(Re, rel_roughness, f0, tol, max_iter):
        f = f0
        for _ in range(max_iter):
            inner = rel_roughness/3.7 + 2.51/(Re*math.sqrt(f))
            g = 1.0/math.sqrt(f) + 2.0*math.log10(inner)
            dg = (-0.5)*f**(-1.5) + (2.0/math.log(10.0))*(1.0/inner)*(2.51/Re)*(-0.5)*f**(-1.5)
            f = f - g/dg
            if abs(g) < tol:
                break
        return f

class TubePressureDrop:
    def __call__(self, L, f, D, rho, v):
        g = f * (rho * v * v / 2) / D
        return g * L
