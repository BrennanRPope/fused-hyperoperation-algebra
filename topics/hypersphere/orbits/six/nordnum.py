import sympy as sp
r, pr, Lm = sp.symbols('r p_r L_m', positive=True)
for mexpr in [1-1/r*sp.Rational(1,2), sp.exp(-1/r), 1 - 1/(3*r**2), 1 + sp.log(r)/7]:
    H = sp.sqrt(mexpr**2 + pr**2 + Lm**2/r**2)
    Hr = sp.diff(H, r).subs(pr, 0); Hrr = sp.diff(H, r, 2).subs(pr, 0); Hpp = sp.diff(H, pr, 2).subs(pr, 0); HL = sp.diff(H, Lm).subs(pr, 0)
    for rv in (1.3, 3.0):
        Ls = sp.nsolve(Hr.subs(r, rv), Lm, 0.5)
        T = float((Hrr*Hpp/HL**2).subs({r: rv, Lm: Ls}))
        E = float(H.subs({pr:0, r: rv, Lm: Ls}))
        v2 = float(((r*HL/1)**2).subs({r: rv, Lm: Ls}))  # coordinate speed = local (conformally flat: c same)
        Phi = sp.log(mexpr)
        N = float((3 + r*sp.diff(Phi, r, 2)/sp.diff(Phi, r)).subs(r, rv))
        print(mexpr, rv, 'T', round(T,6), 'N', round(N,6), 'v2', round(v2,6), 'N/(1-v2)', round(N/(1-v2),6), '(N+v2)/(1-v2)?', round((N - v2)/(1-v2),6), 'N+v2', round(N+v2,6))
