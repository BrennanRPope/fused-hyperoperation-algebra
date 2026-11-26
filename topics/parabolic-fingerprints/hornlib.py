"""Horn-map machinery copied unchanged from ../height-flow/horn_invariants.py (everything above its "E = Germ" line)."""
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

