"""Independent static (v = 0) two-body Lagrangian to fourth order in the sources (2PN static) for a gravitating
sigma model: g = -e^{2 psi} dt^2 + e^{-2 psi/(d-2)} gamma, flat scalars phi, photons A_k = chi_k dt, in d space
dimensions (dimensional regularisation, d -> 3 at the end).  Static action (K = (d-2) Omega_{d-1} -> 4 pi, G = 1):
   L = (1/4K) int [R(gamma) - 2 G_IJ(Phi) dPhi^I.dPhi^J] - sum_A V_A(Phi(x_A)),
   G = diag(rho, 1.., -e^{-2 psi - 2 c_k.phi}),  rho = (d-1)/(2(d-2)),  V_A = m_A(phi) e^psi - q_A.chi.
The 3-metric gamma couples only to the field stress, so its exchange at this order is
   L_h = (1/K^2) sum n_ab n_cd I_h(ab, cd),  n_ab = m_a m_b (1/rho + alpha_a.alpha_b) - q_a.q_b,
i.e. in bidegree (2,2) a combination of n12^2 and n11 n22 only.  This file computes everything EXCEPT L_h:
Fokker action L4 = S4[Phi1] + 1/2 int int J2 eta^-1 J2 P, with P = r^{2-d}.  Space integrals of r1^A r2^B by
Riesz's formula at d = 3 + eps, numerically, and a Laurent fit in eps (poles must cancel; constant term kept)."""
import itertools, mpmath as mp, sympy as sp
from fractions import Fraction
mp.mp.dps = 110
d = sp.Symbol('d'); Ks = sp.Symbol('K')
R = sp.Rational

# ---- monomial fields: dict {(a0, a1, b0, b1): coeff} meaning coeff * r1^(a0 + a1 d) * r2^(b0 + b1 d)
def mono(a0=0, a1=0, b0=0, b1=0, c=1): return {(a0, a1, b0, b1): sp.sympify(c)}
def add(*fs):
    out = {}
    for f in fs:
        for k, v in f.items(): out[k] = out.get(k, 0) + v
    return out
def scale(f, c): return {k: v*c for k, v in f.items()}
def mul(f, g):
    out = {}
    for k1, v1 in f.items():
        for k2, v2 in g.items():
            k = tuple(x + y for x, y in zip(k1, k2)); out[k] = out.get(k, 0) + v1*v2
    return out
ONE = mono()
def P(c): return mono(2, -1, 0, 0) if c == 1 else mono(0, 0, 2, -1)
def ndot(c, e):
    if c == e: return ONE
    return add(mono(1, 0, -1, 0, R(1, 2)), mono(-1, 0, 1, 0, R(1, 2)), mono(-1, 0, -1, 0, -R(1, 2)))
def rpow(c, x0, x1): return mono(x0, x1, 0, 0) if c == 1 else mono(0, 0, x0, x1)
def gg(c, e):  # grad P_c . grad P_e
    return scale(mul(ndot(c, e), mul(rpow(c, 1, -1), rpow(e, 1, -1))), (d - 2)**2)
def Ginv(c, e):  # inverse Laplacian of gg(c, e) (exact in d dims, separation 1)
    if c == e: return scale(mul(P(c), P(c)), R(1, 2))
    return scale(add(mul(P(1), P(2)), scale(P(1), -1), scale(P(2), -1)), R(1, 2))
def Ginv_at(c, e, A):
    if c == e: return 0 if A == c else R(1, 2)
    return -R(1, 2)
def Pat(c, A): return 0 if c == A else 1

# ---- numeric Riesz integral int d^dx r1^al r2^be (separation 1)
def riesz_num(key, dv):
    a0, a1, b0, b1 = key
    if (b0, b1) == (0, 0) or (a0, a1) == (0, 0): return mp.mpf(0)   # scale-free: zero in dim reg
    g = mp.gamma; rg = mp.rgamma
    def one(lam):
        al = a0 + a1*dv + lam; be = b0 + b1*dv + lam
        return (mp.pi**(dv/2)*g((al + dv)/2)*g((be + dv)/2)*g(-(al + be + dv)/2)
                * rg(-al/2)*rg(-be/2)*rg((al + be + 2*dv)/2))
    # second (analytic) regulator for terms whose exponent tracks d exactly; symmetric average drops its simple pole
    return (one(LAM) + one(-LAM))/2
LAM = mp.mpf(10)**-40
def Knum(dv): return (dv - 2)*2*mp.pi**(dv/2)/mp.gamma(dv/2)

class Expr:
    """sum over terms: coeff(params, d, K) * [integral of monomial key] ; key None = no integral."""
    def __init__(self): self.t = {}
    def addc(self, c):
        self.t[None] = self.t.get(None, 0) + c; return self
    def addint(self, f, c=1):
        for k, v in f.items(): self.t[k] = self.t.get(k, 0) + c*v
        return self
    def __iadd__(self, o):
        for k, v in o.t.items(): self.t[k] = self.t.get(k, 0) + v
        return self

E0 = mp.mpf(10)**-20
JS = (-4, -3, -2, -1, 1, 2, 3, 4)
EPS = [mp.mpf(j)*E0 for j in JS]
POW = list(range(-3, 5))
def at3(ex, params, tol=1e-25):
    """d -> 3 value of an Expr, as a sympy polynomial in params with rational coefficients."""
    # expand each coefficient as polynomial in params
    polys = {}
    for k, c in ex.t.items():
        c = sp.expand(c)
        if c == 0: continue
        pc = sp.Poly(c, *params)
        for mon, cf in pc.terms():
            polys.setdefault(mon, []).append((k, cf))
    out = 0
    rcache = {}
    for mon, lst in polys.items():
        fns = [(k, sp.lambdify((d, Ks), cf, 'mpmath')) for k, cf in lst]
        vals = []
        for e in EPS:
            dv = 3 + e; Kv = Knum(dv); s = mp.mpf(0)
            for k, fn in fns:
                if k is None: s += fn(dv, Kv)
                else:
                    if (k, e) not in rcache: rcache[(k, e)] = riesz_num(k, dv)
                    s += fn(dv, Kv)*rcache[(k, e)]
            vals.append(s)
        M = mp.matrix([[mp.mpf(j)**p for p in POW] for j in JS]); cs = mp.lu_solve(M, mp.matrix(vals))
        cs = [cs[i]/E0**p for i, p in enumerate(POW)]
        poles = [cs[i] for i, p in enumerate(POW) if p < 0]
        c0 = cs[POW.index(0)]
        scale_ = max(1, abs(c0))
        if any(abs(pl)*E0**(-POW[i]) > 1e-10*scale_ for i, pl in enumerate(poles)):
            raise ValueError("pole left in %s: %s" % (mon, [mp.nstr(x, 8) for x in poles]))
        fr = Fraction(mp.nstr(c0, 40)).limit_denominator(10**6)
        if abs(c0 - mp.mpf(fr.numerator)/fr.denominator) > 1e-15:
            raise ValueError("not rational: %s %s" % (mon, mp.nstr(c0, 30)))
        out += sp.Rational(fr.numerator, fr.denominator)*sp.prod([p**e for p, e in zip(params, mon)])
    return sp.expand(out)

# ---- the sigma model
def model(nscal, couplings):
    """fields: index 0 = psi, 1..nscal = phi, then one chi per coupling vector.  Returns eta(d), G1[I][J][K], G2[I][J][K][L]."""
    nf = 1 + nscal + len(couplings); rho = (d - 1)/(2*(d - 2))
    psi = 0; ph = list(range(1, 1 + nscal)); ch = list(range(1 + nscal, nf))
    X = sp.symbols('X0:%d' % nf)
    Gm = sp.zeros(nf, nf); Gm[0, 0] = rho
    for s in ph: Gm[s, s] = 1
    for k, cv in enumerate(couplings):
        Gm[ch[k], ch[k]] = -sp.exp(-2*X[0] - 2*sum(cv[s]*X[ph[s]] for s in range(nscal)))
    z = {x: 0 for x in X}
    eta = Gm.subs(z)
    G1 = [[[sp.diff(Gm[i, j], X[k]).subs(z) for k in range(nf)] for j in range(nf)] for i in range(nf)]
    G2 = [[[[sp.diff(Gm[i, j], X[k], X[l]).subs(z) for l in range(nf)] for k in range(nf)] for j in range(nf)] for i in range(nf)]
    return nf, eta, G1, G2, ch

def lagrangian(nscal, couplings, bodyV):
    """bodyV[A] = (V1, V2, V3, V4): derivative tensors of V_A at Phi = 0 (lists / nested lists, sympy entries).
    Returns (L2, L3, L4_noh) as Exprs (L2 at separation 1)."""
    nf, eta, G1, G2, ch = model(nscal, couplings)
    ei = eta.inv()
    F = range(nf); B = (1, 2)
    s = {A: [-sum(ei[I, J]*bodyV[A][0][J] for J in F) for I in F] for A in B}       # Phi1 = sum_A s_A P_A
    Phi1_at = {A: [sum(s[b][I]*Pat(b, A) for b in B) for I in F] for A in B}
    # L2
    L2 = Expr().addc(sp.Rational(1, 2)*sum(-bodyV[A][0][I]*Phi1_at[A][I] for A in B for I in F))
    # L3: bulk cubic + worldline
    L3 = Expr()
    for A in B:
        V2 = bodyV[A][1]; f = Phi1_at[A]
        L3.addc(-R(1, 2)*sum(V2[I][J]*f[I]*f[J] for I in F for J in F))
    for a, b, c in itertools.product(B, B, B):
        cf = sum(G1[I][J][Kk]*s[a][I]*s[b][J]*s[c][Kk] for I in F for J in F for Kk in F)
        if cf != 0: L3.addint(mul(gg(a, b), P(c)), -cf/(2*Ks))
    # L4: quartic action on Phi1
    L4 = Expr()
    for A in B:
        V3 = bodyV[A][2]; f = Phi1_at[A]
        L4.addc(-R(1, 6)*sum(V3[I][J][Kk]*f[I]*f[J]*f[Kk] for I in F for J in F for Kk in F))
    for a, b, c, e in itertools.product(B, B, B, B):
        cf = sum(G2[I][J][Kk][L]*s[c][I]*s[e][J]*s[a][Kk]*s[b][L] for I in F for J in F for Kk in F for L in F)
        if cf != 0: L4.addint(mul(mul(P(a), P(b)), gg(c, e)), -cf/(4*Ks))
    # J2: bulk {(a,b): vector} for gg(a,b), delta {A: vector}
    Jb = {}; Jd = {A: [0]*nf for A in B}
    for a, b in itertools.product(B, B):
        vec = []
        for Lf in F:
            t1 = sum(G1[I][J][Lf]*s[a][I]*s[b][J] for I in F for J in F)
            t2 = sum(G1[Lf][J][Kk]*s[a][Kk]*s[b][J] for J in F for Kk in F)
            vec.append(-t1/(2*Ks) + t2/Ks)
            if a != b: Jd[b][Lf] += -t2
        Jb[(a, b)] = vec
    for A in B:
        V2 = bodyV[A][1]; f = Phi1_at[A]
        for Lf in F: Jd[A][Lf] += -sum(V2[Lf][I]*f[I] for I in F)
    def qf(x, y): return sum(x[I]*ei[I, J]*y[J] for I in F for J in F)
    # 1/2 int int J eta^-1 J P
    for A in B:
        for C in B:
            if A != C: L4.addc(R(1, 2)*qf(Jd[A], Jd[C]))
    for A in B:
        for (c, e), vec in Jb.items():   # cross terms counted twice
            L4.addc(2*R(1, 2)*qf(Jd[A], vec)*(-Ks)*Ginv_at(c, e, A))
    for (a, b), v1 in Jb.items():
        for (c, e), v2 in Jb.items():
            cf = qf(v1, v2)
            if cf != 0: L4.addint(mul(gg(a, b), Ginv(c, e)), R(1, 2)*cf*(-Ks))
    return L2, L3, L4

def emd_body(m, q, al, be, bp):
    """single scalar, single photon: V = m(phi) e^psi - q chi; m'/m = al, (ln m)'' = be, (ln m)''' = bp.
    fields (psi, phi, chi)."""
    mu = [m, m*al, m*(al**2 + be), m*(al**3 + 3*al*be + bp), None]
    def Vd(n):  # n-th derivative tensor
        def entry(idx):
            if idx.count(2) > 0:
                return -q if (len(idx) == 1) else 0
            k = idx.count(1); return mu[k]
        if n == 1: return [entry((i,)) for i in range(3)]
        if n == 2: return [[entry((i, j)) for j in range(3)] for i in range(3)]
        if n == 3: return [[[entry((i, j, k)) for k in range(3)] for j in range(3)] for i in range(3)]
    return (Vd(1), Vd(2), Vd(3), None)
