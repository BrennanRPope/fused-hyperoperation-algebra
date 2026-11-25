"""Local invariants of e^x - 1 at its double fixed point, tied to the rung height (DRAFT section).

Every germ f(x) = x + a1 x^2 + a2 x^3 + ... (a1 != 0) fades like -1/(a1 n). What tells such germs apart:
  * the iterative residue  resit(f) = 1 - a2/a1^2          (formal invariant; Ecalle; Milnor sec. 12)
  * the Ecalle-Voronin invariants: the two horn maps h+(W), h-(W) built from the attracting and repelling
    Fatou coordinates.  h(W) = W + A0 + sum A_n e^{2 pi i n W}.  f embeds in an analytic flow
    (all real heights analytic at 0) exactly when both horn maps are translations (Ecalle, Voronin).

We work in Z = -1/(a1 x), where f becomes F(Z) = Z + 1 + rho/Z + ...  with rho = resit(f).
Canonical Fatou coordinate:  Phi0(Z) = Z - rho log Z + sum c_j Z^-j  (c_j solved from Phi0(F) = Phi0 + 1);
attracting  Phi_att(Z) = Phi0(F^n Z) - n,   repelling  Phi_rep(Z) = Phi0(F^-n Z) + n, the log branch kept
continuous on each petal (repelling petal: arg in (0, 2 pi)).
Writes horn_invariants.json and prints PASS/FAIL lines.
"""
import os, sys, json, math
import mpmath as mp
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
mp.mp.dps = 34
M = 12          # terms of the asymptotic Fatou series
res = {}
def rep(name, ok, extra=""):
    print(("PASS " if ok else "FAIL ") + name + (("   " + extra) if extra else ""))
    res[name] = bool(ok)

from fractions import Fraction as Fr
x = sp.symbols("x")
N = M + 4
def s_mul(a, b):
    out = [Fr(0)]*(N + 1)
    for i, ai in enumerate(a):
        if ai == 0: continue
        for j in range(0, N + 1 - i):
            if b[j]: out[i + j] += ai*b[j]
    return out
def s_inv(a):          # 1/a, a[0] != 0
    out = [Fr(0)]*(N + 1); out[0] = 1/a[0]
    for n in range(1, N + 1):
        out[n] = -sum(a[k]*out[n - k] for k in range(1, n + 1))/a[0]
    return out
def s_log1(a):         # log(a), a[0] = 1:  (log a)' = a'/a
    da = [(i + 1)*a[i + 1] for i in range(N)] + [Fr(0)]
    q = s_mul(da, s_inv(a))
    return [Fr(0)] + [q[i - 1]/i for i in range(1, N + 1)]
def s_pow(a, e):       # a^e for integer e (a[0] = 1)
    out = [Fr(1)] + [Fr(0)]*N; base = a if e >= 0 else s_inv(a)
    for _ in range(abs(e)): out = s_mul(out, base)
    return out
def compose_series(f, g):   # f(g(x)), g[0] = 0
    out = [Fr(0)]*(N + 1); gp = [Fr(1)] + [Fr(0)]*N
    for i in range(1, N + 1):
        gp = s_mul(gp, g)
        if f[i]: out = [o + f[i]*v for o, v in zip(out, gp)]
    return out

def fatou_series(co):
    """co[i] = coefficient of x^i of f (co[1] = 1), Fractions. Returns (a1, rho, [c_1..c_M]) exactly."""
    a1 = co[2]
    den = [Fr(1)] + [co[j + 1]*(-1/a1)**j for j in range(1, N)] + [Fr(0)]
    P = s_inv(den)                                 # F(Z) = Z * P(u), u = 1/Z
    rho = P[2]
    logP = s_log1(P)
    # base = (P - 1)/u - rho log P  (series in u); then add c_j u^j (P^-j - 1); require = 1 + 0 u + ...
    base = [P[i + 1] for i in range(N)] + [Fr(0)]
    base = [b - rho*l for b, l in zip(base, logP)]
    cs = []
    tot = base[:]
    for m in range(1, M + 1):
        # order u^(m+1): tot[m+1] + c_m * (coefficient of u^(m+1) in u^m (P^-m - 1)) = 0  -> that coefficient is -m
        Pm = s_pow(P, -m)
        term = [Fr(0)]*(N + 1)
        for i in range(0, N + 1 - m):
            term[i + m] = Pm[i] - (1 if i == 0 else 0)
        cm = -tot[m + 1]/term[m + 1]
        cs.append(cm)
        tot = [t_ + cm*v for t_, v in zip(tot, term)]
    assert tot[0] == 1 and tot[1] == 0 and all(tot[i] == 0 for i in range(2, M + 2))
    return a1, rho, cs

def exp_series(scale=Fr(1)):
    return [scale**i/math.factorial(i) for i in range(N + 1)]
CO_E = [Fr(0)] + [Fr(1, math.factorial(i)) for i in range(1, N + 1)]

class Germ:
    def __init__(self, name, co, f, finv):
        self.name = name
        assert co[0] == 0 and co[1] == 1, name
        self.a1, self.rho, cs = fatou_series(co)
        self.resit_formula = 1 - co[3]/self.a1**2
        fr = lambda q: mp.mpf(q.numerator)/q.denominator
        self.c = [fr(cj) for cj in cs]
        self.A1 = fr(self.a1); self.rhof = fr(self.rho)
        self.f, self.finv = f, finv
    def F(self, Z):  return -1/(self.A1*self.f(-1/(self.A1*Z)))
    def Fi(self, Z): return -1/(self.A1*self.finv(-1/(self.A1*Z)))
    def phi0(self, Z, branch=0):
        L = mp.log(Z)
        if branch and mp.im(Z) < 0: L += 2j*mp.pi       # repelling petal: arg in (0, 2pi)
        return Z - self.rhof*L + sum(cj*Z**(-j) for j, cj in enumerate(self.c, start=1))
    def phi_att(self, Z, R=900):
        n = 0
        while abs(Z) < R or mp.re(Z) < 0.5*abs(Z):
            Z = self.F(Z); n += 1
            if n > 200000: raise RuntimeError("no convergence")
        return self.phi0(Z) - n
    def phi_rep(self, Z, R=900):
        n = 0
        while abs(Z) < R or mp.re(Z) > -0.5*abs(Z):
            Z = self.Fi(Z); n += 1
            if n > 200000: raise RuntimeError("no convergence")
        return self.phi0(Z, branch=1) + n

def horn(g, Y, sign=+1, npts=24, nmodes=5):
    """Fit h(W) - W = A0 + sum_{m=1..nmodes} A_{sign*m} q^m, q = exp(2 pi i sign W), from points Z = X + i*sign*Y."""
    Ws, Hs = [], []
    for k in range(npts):
        Z = mp.mpc(-1.5 + 3.0*k/npts, sign*Y)
        Ws.append(g.phi_rep(Z)); Hs.append(g.phi_att(Z))
    imW = sum(mp.im(W) for W in Ws)/len(Ws)
    qs = mp.exp(-2*mp.pi*sign*imW)                       # scale each mode so the fit is well conditioned
    rows = []; rhs = []
    for W, H in zip(Ws, Hs):
        q = mp.exp(2j*mp.pi*sign*W)/qs
        rows.append([q**m for m in range(0, nmodes + 1)]); rhs.append(H - W)
    sol, _ = mp.qr_solve(mp.matrix(rows), mp.matrix(rhs))
    resid = max(abs(sum(rows[i][m]*sol[m] for m in range(nmodes + 1)) - rhs[i]) for i in range(len(rows)))
    sol = [sol[m]/qs**m for m in range(nmodes + 1)]
    return sol, resid, imW

E = Germ("e^x - 1", CO_E, mp.expm1, mp.log1p)
cplx = lambda z: mp.mpc(z)

# ---------- 1. iterative residue ----------
rep("resit(e^x - 1) = 1 - a2/a1^2 = 1 - (1/6)/(1/4) = 1/3, and equals the Fatou log coefficient rho (= the running-coupling two-loop invariant)",
    E.rho == Fr(1, 3) and E.resit_formula == Fr(1, 3))
print("Fatou series c_j for e^x - 1 (Z = -2/x):", [str(c) for c in fatou_series(CO_E)[2][:6]])

# heights: formal E^s, coefficients exact polynomials in s
s = sp.symbols("s")
sys.path.insert(0, os.path.join(HERE, "..", "fractional-rungs"))
from regular_iteration import ITS
def formal_coeff(i):
    pts = list(range(0, i)); val = 0
    for j in pts:
        w = 1
        for m in pts:
            if m != j: w *= (s - m)/sp.Integer(j - m)
        val += sp.Rational(ITS[j][i].numerator, ITS[j][i].denominator)*w
    return sp.expand(val)
c2, c3 = formal_coeff(2), formal_coeff(3)
resit_s = sp.simplify(1 - c3/c2**2)
print("formal E^s = x + (", c2, ") x^2 + (", c3, ") x^3 + ...;  resit(E^s) =", resit_s)
rep("resit(E^s) = 1/(3s) at every height s: height times iterative residue is the constant 1/3", sp.simplify(resit_s - sp.Rational(1, 3)/s) == 0)

# corner dial: (e^{kx}-1)/k is conjugate to E by scaling, so resit = 1/3 for every k != 0
k, d = sp.symbols("k d", nonzero=True)
co_k = [sp.series((sp.exp(k*x) - 1)/k, x, 0, 4).removeO().coeff(x, i) for i in range(4)]
rk = sp.simplify(1 - co_k[3]/co_k[2]**2)
rep("corner dial: resit((e^{kx}-1)/k) = 1/3 for every k (conjugate to e^x - 1 by x -> kx; the dial is invisible locally)", rk == sp.Rational(1, 3))
# plane dial: group exponential of (a + b + k ab)/(1 - d ab) solves G' = 1 + kG + dG^2
t = sp.symbols("t")
Gser = [0, 1]
for n_ in range(2, 5):
    # coefficient recursion from G' = 1 + k G + d G^2
    Gp = sum(sp.Rational(1)*i*Gser[i]*t**(i - 1) for i in range(1, len(Gser)))
    Gt = sum(Gser[i]*t**i for i in range(len(Gser)))
    rhs_ = sp.expand(1 + k*Gt + d*Gt**2)
    Gser.append(sp.simplify(rhs_.coeff(t, n_ - 1)/n_))
rp = sp.simplify(1 - Gser[3]/Gser[2]**2)
kappa = d - k**2/4
print("plane dial: G = t + (", Gser[2], ") t^2 + (", sp.factor(Gser[3]), ") t^3;  resit =", sp.factor(rp))
rep("plane dial: resit(G_{k,d}) = 1/3 - 4d/(3k^2) = -(4/3) kappa/k^2 with kappa = d - k^2/4 (the Cayley-Klein number); zero exactly at kappa = 0",
    sp.simplify(rp - (sp.Rational(1, 3) - 4*d/(3*k**2))) == 0 and sp.simplify(rp + sp.Rational(4, 3)*kappa/k**2) == 0)

# ---------- 2. horn maps of e^x - 1 ----------
hp = {}
for Y in (2.0, 2.5, 3.0):
    A, r_, imW = horn(E, Y, +1)
    hp[Y] = (A, r_, imW)
    print(f"upper horn, Im Z = {Y}: Im W ~ {mp.nstr(imW, 4)}; A0 = {mp.nstr(A[0], 6)}, A1 = {mp.nstr(A[1], 12)}, A2 = {mp.nstr(A[2], 8)}, fit resid {mp.nstr(r_, 3)}")
A1s = [hp[Y][0][1] for Y in hp]
spread = max(abs(a - A1s[-1]) for a in A1s)/abs(A1s[-1])
rep("upper horn map of e^x - 1 is NOT a translation: A1 is nonzero and the same from three strips (Im Z = 2, 2.5, 3)",
    abs(A1s[-1]) > 1e-3 and spread < 1e-8, f"|A1| = {mp.nstr(abs(A1s[-1]), 12)}, spread {mp.nstr(spread, 2)}")
rep("upper constant A0+ = 0 in the canonical normalisation", abs(hp[3.0][0][0]) < 1e-12, mp.nstr(hp[3.0][0][0], 3))
hm = horn(E, 3.0, -1)
print(f"lower horn: A0- = {mp.nstr(hm[0][0], 12)}, A_-1 = {mp.nstr(hm[0][1], 12)}")
rep("lower constant A0- = 2 pi i * resit = 2 pi i / 3 (the known link between the two horn constants and the residue)",
    abs(hm[0][0] - 2j*mp.pi/3) < 1e-10, mp.nstr(hm[0][0] - 2j*mp.pi/3, 3))
inv = A1s[-1]*hm[0][1]
pred_m = abs(A1s[-1])*mp.exp(4*mp.pi**2*E.rhof)
rep("real map: the lower horn mirrors the upper once the 2 pi i/3 branch shift is undone: |A_-1| = |A1| e^(4 pi^2/3); the shift-free invariant A1*A_-1 is real",
    abs(abs(hm[0][1]) - pred_m) < 1e-10*pred_m and abs(mp.im(inv)) < 1e-15*abs(inv),
    f"A1*A_-1 = {mp.nstr(mp.re(inv), 15)}, sqrt = {mp.nstr(mp.sqrt(mp.re(inv)), 12)}")
A1E = A1s[-1]
res["E"] = {"resit": "1/3", "A1_upper": [str(mp.re(A1E)), str(mp.im(A1E))], "abs_A1": str(abs(A1E)),
            "A2_upper": str(abs(hp[3.0][0][2])), "A0_lower": str(hm[0][0]), "A1_times_Am1": str(inv)}

# ---------- 3. tie to height ----------
# (a) integer height 2: E o E. Expect |A2(E^2)| = |A1(E)|/2 and A1(E^2) = 0.
E2 = Germ("E o E", compose_series(CO_E, CO_E), lambda z: mp.expm1(mp.expm1(z)), lambda z: mp.log1p(mp.log1p(z)))
rep("resit(E^2) = 1/6 (= 1/(3*2))", E2.rho == Fr(1, 6))
A2h, r2, _ = horn(E2, 3.0, +1, npts=28, nmodes=6)
print(f"E^2 upper horn: A1 = {mp.nstr(A2h[1], 4)}, |A2| = {mp.nstr(abs(A2h[2]), 12)}, |A1(E)|/2 = {mp.nstr(abs(A1E)/2, 12)}")
rep("height 2: the horn map of E^2 keeps only even modes and |A2(E^2)| = |A1(E)|/2  (h_2(W) = h(2W - c)/2 + c)",
    abs(A2h[1]) < 1e-8*abs(A1E) and abs(abs(A2h[2]) - abs(A1E)/2) < 1e-6*abs(A1E))
res["E2"] = {"abs_A1": str(abs(A2h[1])), "abs_A2": str(abs(A2h[2]))}

# (b) fractional heights: the half-rung built from the attracting side and from the repelling side disagree
#     by  h(W+s) - h(W) - s = sum A_n e^{2 pi i n W}(e^{2 pi i n s} - 1)  ~  2 |A1| |sin(pi s)| e^{-2 pi Im W}
def psi_rep(W, Z0):
    return mp.findroot(lambda Z: E.phi_rep(Z) - W, Z0)
Zp = mp.mpc(0.3, 2.5)
W0 = E.phi_rep(Zp); Hp0 = E.phi_att(Zp)
mism = {}
for sv in [0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0, 1.25, 1.5, 2.0]:
    Zs = psi_rep(W0 + sv, Zp + sv)
    D = E.phi_att(Zs) - (Hp0 + sv)
    pred = sum(hp[3.0][0][m]*mp.exp(2j*mp.pi*m*W0)*(mp.exp(2j*mp.pi*m*sv) - 1) for m in range(1, 5))
    mism[sv] = (D, pred)
    print(f"  s = {sv:4}: attracting-vs-repelling half-rung gap |D| = {mp.nstr(abs(D), 8)}   predicted {mp.nstr(abs(pred), 8)}")
ok = all(abs(D - P) < 1e-6*abs(A1E)*abs(mp.exp(2j*mp.pi*W0)) + 1e-25 for D, P in mism.values())
ok0 = all(abs(mism[sv][0]) < 1e-15 for sv in (0.0, 1.0, 2.0))
rep("fractional heights: the two petals' E^s disagree by exactly the horn coefficients; the gap is ~ 2|A1||sin(pi s)| e^{-2 pi Im W}, zero only at whole heights",
    ok and ok0)
res["half_rung_gap"] = {str(k_): [str(abs(a)), str(abs(b))] for k_, (a, b) in mism.items()}
res["half_rung_gap_point"] = {"Z": str(Zp), "ImW": str(mp.im(W0))}

# ---------- 4. other germs with the same fading rate (all x + x^2/2 + ...) ----------
def plane_coeffs(kk, dd):
    G = [Fr(0), Fr(1)] + [Fr(0)]*(N - 1)
    for n_ in range(2, N + 1):
        G2 = s_mul(G, G)
        G[n_] = (kk*G[n_ - 1] + dd*G2[n_ - 1])/n_      # from G' = 1 + kG + dG^2
    return G
def plane(kk, dd):
    kk, dd = mp.mpf(kk), mp.mpf(dd)
    r = mp.sqrt(kk**2 - 4*dd)
    if abs(r) < mp.mpf(10)**(-30):
        f = lambda z: z/(1 - kk*z/2); finv = lambda y: y/(1 + kk*y/2)
    else:
        f = lambda z: 2*mp.sinh(r*z/2)/(r*mp.cosh(r*z/2) - kk*mp.sinh(r*z/2))
        # group logarithm: integral of dt/(1 + k t + d t^2), closed form
        def finv(y):
            return (2/r)*mp.atanh(r*y/(2 + kk*y))
    return f, finv
def cubic_inv(y):
    return mp.findroot(lambda z: z + z**2/2 + z**3/6 - y, y - y**2/2 + y**3/3)
half = Fr(1, 2)
co_mob = [Fr(0)] + [half**(i - 1) for i in range(1, N + 1)]
co_cut = [Fr(0), Fr(1), Fr(1, 2), Fr(1, 6)] + [Fr(0)]*(N - 3)
co_xe = [Fr(0)] + [half**(i - 1)/math.factorial(i - 1) for i in range(1, N + 1)]
co_log = [Fr(0)] + [Fr(1, i) for i in range(1, N + 1)]
co_quad = [Fr(0), Fr(1), Fr(1, 2)] + [Fr(0)]*(N - 2)
others = [
    ("Möbius x/(1 − x/2)", co_mob, lambda z: z/(1 - z/2), lambda z: z/(1 + z/2)),
    ("x + x²/2 + x³/6 (e^x − 1 cut at x³)", co_cut, lambda z: z + z**2/2 + z**3/6, cubic_inv),
    ("x·e^(x/2)", co_xe, lambda z: z*mp.exp(z/2), lambda z: 2*mp.lambertw(z/2)),
    ("−log(1 − x)", co_log, lambda z: -mp.log(1 - z), lambda z: -mp.expm1(-z)),
    ("x + x²/2", co_quad, lambda z: z + z**2/2, lambda z: mp.sqrt(1 + 2*z) - 1),
    ("plane dial k = 1, d = −1/4", plane_coeffs(Fr(1), Fr(-1, 4))) + plane(1, -0.25),
    ("plane dial k = 1, d = 1/4 (κ = 0)", plane_coeffs(Fr(1), Fr(1, 4))) + plane(1, 0.25),
]
table = [("e^x − 1 (our height-1 map)", Fr(1, 3), abs(A1E), mp.sqrt(abs(inv)))]
for name, co, f, finv in others:
    g = Germ(name, co, f, finv)
    try:
        A, r_, imW = horn(g, 3.0, +1)
        Ab, rb, _ = horn(g, 2.5, +1)
        Am, rm, _ = horn(g, 3.0, -1)
        prod = A[1]*Am[1]
        a1 = abs(A[1]); sym = mp.sqrt(abs(prod))
        stable = abs(A[1] - Ab[1])/max(abs(A[1]), mp.mpf(10)**-30)
        print(f"  {name}: A1 = {mp.nstr(A[1], 10)}, A_-1 = {mp.nstr(Am[1], 10)}, sqrt|A1 A_-1| = {mp.nstr(sym, 10)}, strip spread {mp.nstr(stable, 2)}, fit resid {mp.nstr(max(r_, rm), 2)}", flush=True)
    except Exception as ex:
        a1 = sym = mp.nan; print("   horn failed for", name, ex)
    table.append((name, g.rho, a1, sym))
print("germs that all fade like -2/n:")
for name, rr, a1, sym in table:
    print(f"   {name:42s} resit = {str(rr):6s}  |A1+| = {mp.nstr(a1, 10):14s} shift-free size sqrt|A1 A_-1| = {mp.nstr(sym, 10)}")
mob = [t_ for t_ in table if t_[0].startswith("Möbius")][0]
cut = [t_ for t_ in table if "cut" in t_[0]][0]
k0 = [t_ for t_ in table if "κ = 0" in t_[0]][0]
rep("Möbius x/(1 - x/2): resit 0 and trivial horn map (A1 = 0): it embeds in an analytic flow, unlike e^x - 1", mob[1] == 0 and mob[2] < 1e-12)
rep("plane dial at kappa = 0 is that Möbius map: every height analytic there", k0[1] == 0 and k0[2] < 1e-12)
rep("e^x - 1 cut at x^3 has the SAME fading rate and the SAME residue 1/3, yet a different horn map: only the horn map tells them apart",
    cut[1] == Fr(1, 3) and abs(cut[3] - table[0][3]) > 1e-3*table[0][3])
res["table"] = [[n_, str(r_), str(a_), str(b_)] for n_, r_, a_, b_ in table]

json.dump(res, open(os.path.join(HERE, "horn_invariants.json"), "w"), indent=1, default=str)
nb = sum(v for v in res.values() if isinstance(v, bool)); nt = sum(isinstance(v, bool) for v in res.values())
print(f"\n{nb} of {nt} checks pass")
