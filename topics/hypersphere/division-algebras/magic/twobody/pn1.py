"""First post-Newtonian two-body orbit turn for bodies that trade gravitons, scalars and vectors.
Harmonic-gauge 1PN Lagrangian, general strengths (G = M = c = 1):
  gravity  (1/r)[1 + 3/2(vA^2+vB^2) - 7/2 vA.vB - 1/2 (n.vA)(n.vB)]           (EIH)
  scalars  (s/r)[1 - 1/2(vA^2+vB^2) + 1/2 vA.vB - 1/2 (n.vA)(n.vB)]           (DEF)
  vectors -(v/r)[1 - 1/2 (vA.vB + (n.vA)(n.vB))]                              (Darwin)
  O(G^2)  -(mA mB/(2 r^2)) (mA C_A + mB C_B)   (each term has one body sourcing twice)
plus the free terms m v^2/2 + m v^4/8.  T = turn/(pi K^2/J^2), K = 1 + s - v, J per reduced mass."""
import sympy as sp
r, rd, ph, eps, nu, s, v, CA, CB, XA = sp.symbols('r rdot phidot epsilon nu s v C_A C_B X_A', real=True)
E, J, u, K = sp.symbols('E J u K', positive=True)
XB = 1 - XA
def L1_rel(nu_, XA_):
    v2 = rd**2 + r**2*ph**2; nv2 = rd**2
    grav = (sp.Rational(3, 2) + nu_/2)*v2 + nu_/2*nv2
    scal = s*(-(1 - nu_)/2*v2 + nu_/2*nv2)
    vect = -v*nu_/2*(v2 + nv2)
    kin = (1 - 3*nu_)/8*v2**2
    stat = -(XA_*CA + (1 - XA_)*CB)/(2*r**2)
    return kin + (grav + scal + vect)/r + stat
def turn_T(nu_, XA_):
    Kv = 1 + s - v
    L1 = L1_rel(nu_, XA_)
    Jexpr = r**2*ph + eps*sp.diff(L1, ph)
    Eexpr = (rd**2 + r**2*ph**2)/2 - Kv/r + eps*(rd*sp.diff(L1, rd) + ph*sp.diff(L1, ph) - L1)
    # Newtonian relations to substitute inside O(eps) terms
    phN = J/r**2; rd2N = 2*E + 2*Kv/r - J**2/r**2
    dLph = sp.expand(sp.diff(L1, ph).subs(ph, phN)); dLph = sp.expand(dLph.subs(rd**2, rd2N)); assert rd not in dLph.free_symbols
    phdot = J/r**2 - eps*dLph/r**2
    # rdot^2 from energy: E = rd^2/2 + r^2 ph^2/2 - K/r + eps*(...)
    corr = (rd*sp.diff(L1, rd) + ph*sp.diff(L1, ph) - L1)
    corrN = sp.expand(corr.subs(ph, phN)).subs(rd**2, rd2N).subs(rd**4, rd2N**2)
    corrN = sp.expand(corrN)
    assert rd not in corrN.free_symbols, corrN
    rd2 = 2*E + 2*Kv/r - r**2*phdot**2 - 2*eps*corrN
    F = sp.series(rd2/(r**4*phdot**2), eps, 0, 2).removeO()
    F = sp.expand(F.subs(r, 1/u))
    F0 = F.coeff(eps, 0); g = sp.expand(F.coeff(eps, 1))
    # F0 = c0 + c1 u - u^2 ; u = p + q cos(theta)
    c0 = sp.expand(F0.coeff(u, 0)); c1 = sp.expand(F0.coeff(u, 1)); assert sp.simplify(F0.coeff(u, 2) + 1) == 0
    p, q2, th = sp.symbols('p q2 theta')
    tot = 0
    for k in range(0, 7):
        gk = g.coeff(u, k)
        if gk == 0: continue
        Ik = sp.integrate(sp.expand((p + sp.sqrt(q2)*sp.cos(th))**k), (th, 0, sp.pi))
        tot += gk*sp.diff(Ik, q2)
    turn = sp.simplify(2*tot.subs({p: c1/2, q2: c0 + c1**2/4}))
    T = sp.simplify(turn/(sp.pi*Kv**2/J**2))
    return sp.simplify(T)
if __name__ == "__main__":
    Tgen = turn_T(nu, XA)
    print("general T =", sp.factor(Tgen))
    print("GR (s=v=0, C=1):", sp.simplify(Tgen.subs({s:0,v:0,CA:1,CB:1})))
    # DEF scalar-tensor: s = aA aB, C_A = (1+s)^2 + bB aA^2  (body A sources twice; normalisation fixed by check_def_probe.py)
    aA,aB,bA,bB=sp.symbols('alpha_A alpha_B beta_A beta_B')
    ss=aA*aB
    TDEF=Tgen.subs({v:0,s:ss,CA:(1+ss)**2+bB*aA**2,CB:(1+ss)**2+bA*aB**2})
    known=6*((1-ss/3)/(1+ss)-(XA*bB*aA**2+(1-XA)*bA*aB**2)/(6*(1+ss)**2))
    print("DEF check (0 = match):", sp.simplify(TDEF-known))
    print("depends on nu?", sp.diff(Tgen,nu))
