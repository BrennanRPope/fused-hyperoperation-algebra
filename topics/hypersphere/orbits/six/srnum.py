import sympy as sp
r, pr, Lm = sp.symbols('r p_r L_m', positive=True)
for Vexpr in [-1/r, -1/r**2*sp.Rational(1,3), -1/r + sp.Rational(3,10)/r**2, sp.log(r)/5]:
    H = sp.sqrt(1 + pr**2 + Lm**2/r**2) + Vexpr
    Hr = sp.diff(H, r).subs(pr, 0); Hrr = sp.diff(H, r, 2).subs(pr, 0); Hpp = sp.diff(H, pr, 2).subs(pr, 0); HL = sp.diff(H, Lm).subs(pr, 0)
    for rv in (0.7, 2.0):
        Ls = sp.nsolve(Hr.subs(r, rv), Lm, 1.0)
        T = float((Hrr*Hpp/HL**2).subs({r: rv, Lm: Ls}))
        v2 = float(((r*HL)**2).subs({r: rv, Lm: Ls}))
        V1 = float(sp.diff(Vexpr, r).subs(r, rv)); V2 = float(sp.diff(Vexpr, r, 2).subs(r, rv))
        clos = 3 + rv*V2/V1
        print(Vexpr, rv, 'T', round(T, 6), 'clos-v2', round(clos - v2, 6), 'clos*(1-v2)', round(clos*(1-v2), 6), 'v2', round(v2,4))
