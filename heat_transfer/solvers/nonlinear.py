# solvers/nonlinear.py
def newton_2x2(F, J, x0, tol, max_iters):
    x = x0
    for _ in range(max_iters):
        f1,f2 = F(x[0],x[1])
        if max(abs(f1),abs(f2))<tol: return x
        a,b,c,d = J(x[0],x[1])
        det = a*d-b*c
        dx = ( d*f1 - b*f2)/det
        dy = (-c*f1 + a*f2)/det
        x = (x[0]-dx, x[1]-dy)
    return x
