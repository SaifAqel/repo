# solvers/rootfind.py
def bisection(f, a, b, tol, max_iters):
    fa = f(a); fb = f(b)
    if fa*fb>0: raise ValueError("sign")
    for _ in range(max_iters):
        c = 0.5*(a+b)
        fc = f(c)
        if abs(fc)<tol or 0.5*(b-a)<tol: return c
        if fa*fc<0: b= c; fb=fc
        else: a=c; fa=fc
    return 0.5*(a+b)
def brent(f,a,b,tol,max_iters):
    from math import copysign
    fa=f(a); fb=f(b)
    if fa*fb>0: raise ValueError("sign")
    c=a; fc=fa; d=e=b-a
    for _ in range(max_iters):
        if fb==0: return b
        if fa*fb>0: a=c; fa=fc; d=e=b-a
        if abs(fa)<abs(fb): c=b; b=a; a=c; fc=fb; fb=fa; fa=fc
        m=0.5*(a-b); tol1=2.220446049250313e-16*abs(b)+0.5*tol
        if abs(m)<=tol1: return b
        if abs(e)>=tol1 and abs(fc)>abs(fb):
            s=fb/fc; 
            if a==c: p=2*m*s; q=1-s
        else:
                q=fc/fa; r=fb/fa; p=s*(2*m*q*(q-r)-(b-c)*(r-1)); q=(q-1)*(r-1)*(s-1)
        if p>0: q=-q
        p=abs(p); min1=3*m*q-abs(tol1*q); min2=abs(e*q)
        if 2*p<min(min1,min2): e=d; d=p/q
        else: d=m; e=m
    else: d=m; e=m
    c=b; fc=fb
    if abs(d)>tol1: b+=d
    else: b+=copysign(tol1,m)
    fb=f(b)
    return b
