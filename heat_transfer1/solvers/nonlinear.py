# solvers/nonlinear.py
from common.units import ureg, Q_

def newton_2x2(F, J, x0, tol, max_iters):
    # Ensure inputs are quantities
    x = (Q_(x0[0]), Q_(x0[1]))
    tol = Q_(tol)

    for _ in range(max_iters):
        f1, f2 = F(x[0], x[1])

        # Convergence check
        if max(abs(f1), abs(f2)) < tol:
            return x

        a, b, c, d = J(x[0], x[1])

        det = a * d - b * c
        dx = (d * f1 - b * f2) / det
        dy = (-c * f1 + a * f2) / det

        x = (x[0] - dx, x[1] - dy)

    return x
