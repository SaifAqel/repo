# solvers/rootfind.py
from common.units import ureg, Q_
_EPS = Q_(2.220446049250313e-16, ureg.dimensionless)

def bisection(f, a, b, tol, max_iters):
    a = Q_(a)
    b = Q_(b)
    tol = Q_(tol)

    fa = f(a)
    fb = f(b)
    if fa * fb > Q_(0.0, ureg.dimensionless):
        raise ValueError("sign")

    for _ in range(max_iters):
        c = Q_(0.5, ureg.dimensionless) * (a + b)
        fc = f(c)
        if abs(fc) < tol or Q_(0.5, ureg.dimensionless) * (b - a) < tol:
            return c
        if fa * fc < Q_(0.0, ureg.dimensionless):
            b = c
            fb = fc
        else:
            a = c
            fa = fc
    return Q_(0.5, ureg.dimensionless) * (a + b)


def brent(f, a, b, tol, max_iters):
    a = Q_(a)
    b = Q_(b)
    tol = Q_(tol)

    fa = f(a)
    fb = f(b)
    if fa * fb > Q_(0.0, ureg.dimensionless):
        raise ValueError("sign")

    c = a
    fc = fa
    d = e = b - a

    for _ in range(max_iters):
        if fb == Q_(0.0, ureg.dimensionless):
            return b

        if fa * fb > Q_(0.0, ureg.dimensionless):
            a = c
            fa = fc
            d = e = b - a

        if abs(fa) < abs(fb):
            c, b, a = b, a, c
            fc, fb, fa = fb, fa, fc

        m = Q_(0.5, ureg.dimensionless) * (a - b)
        tol1 = _EPS * abs(b) + Q_(0.5, ureg.dimensionless) * tol

        if abs(m) <= tol1:
            return b

        if abs(e) >= tol1 and abs(fc) > abs(fb):
            s = fb / fc
            if a == c:
                p = Q_(2.0, ureg.dimensionless) * m * s
                q = Q_(1.0, ureg.dimensionless) - s
            else:
                q = fc / fa
                r = fb / fa
                p = s * (Q_(2.0, ureg.dimensionless) * m * q * (q - r) - (b - c) * (r - Q_(1.0, ureg.dimensionless)))
                q = (q - Q_(1.0, ureg.dimensionless)) * (r - Q_(1.0, ureg.dimensionless)) * (s - Q_(1.0, ureg.dimensionless))

            if p > Q_(0.0, ureg.dimensionless):
                q = -q
            p = abs(p)
            min1 = Q_(3.0, ureg.dimensionless) * m * q - abs(tol1 * q)
            min2 = abs(e * q)

            if Q_(2.0, ureg.dimensionless) * p < min(min1, min2):
                e = d
                d = p / q
            else:
                d = m
                e = m
        else:
            d = m
            e = m

        c = b
        fc = fb

        if abs(d) > tol1:
            b = b + d
        else:
            b = b + (tol1 if m >= Q_(0.0, ureg.dimensionless) else -tol1)

        fb = f(b)

    return b
