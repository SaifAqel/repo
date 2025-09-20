from scipy.optimize import root_scalar

def solve_brentq(func, bracket, args, xtol):
    sol = root_scalar(func, bracket=bracket, method="brentq", xtol=xtol, args=args)
    return sol.root
