"""Solitons on the rungs, and what the gluing number q does to them.

Every claim on the page is checked here. Output: checks_output.txt, results.json.
Conventions: KdV  u_t + 6 u u_x + u_xxx = 0,  u = 2 d^2/dx^2 log tau  (Hirota 1971).
Rung s law  B_s(a,b) = E^s(E^-s(a) + E^-s(b)),  E(x) = e^x - 1  (s=0 add, s=1 mult-add).
Weighted mult-add  a + b + k ab.   Theta  Th(w;q) = sum_n q^{n(n-1)/2} w^n.
"""
import sys, json, math, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, "/mnt/project-files/math-theory/fractional-rungs")
import numpy as np
import sympy as sp
import mpmath as mp
from regular_iteration import Es

RES = {}
OUT = []
def rep(name, ok, extra=""):
    line = ("PASS " if ok else "FAIL ") + name + (("   " + extra) if extra else "")
    print(line); OUT.append(line); RES[name] = {"pass": bool(ok), "detail": extra}
def note(s):
    print("   " + s); OUT.append("   " + s)

# ---------------------------------------------------------------------------
print("== 1. The three-body test: which rung has no three-body term ==")
s, a, b, c, x = sp.symbols("s a b c x")
ORDER = 3
def compose(f, g):
    P = sp.Poly(sp.expand(f.subs(x, g)), x)
    return sp.expand(sum(cf * x**mon[0] for mon, cf in P.terms() if mon[0] <= ORDER))
Ex = sum(x**i / sp.factorial(i) for i in range(1, ORDER + 1))
its = [x]
for n in range(1, ORDER + 1):
    its.append(sp.expand(compose(Ex, its[-1])))
Es_series = x
for i in range(2, ORDER + 1):
    pts = list(range(i))
    coeff = 0
    for j in pts:
        w = 1
        for m in pts:
            if m != j: w *= (s - m) / sp.Integer(j - m)
        coeff += its[j].coeff(x, i) * w
    Es_series += sp.factor(coeff) * x**i
def Esym(sv, arg):
    return Es_series.subs(s, sv).subs(x, arg)
t = sp.symbols("t")
def trunc(expr, vars_, deg):
    P = sp.Poly(sp.expand(expr), *vars_)
    return sp.expand(sum(cf * sp.prod([v**e for v, e in zip(vars_, mon)])
                         for mon, cf in P.terms() if sum(mon) <= deg))
Bs_ab = trunc(Esym(s, Esym(-s, a) + Esym(-s, b)), [a, b], 3)
note("B_s(a,b) to degree 3: " + str(sp.collect(Bs_ab, [a*b], sp.factor)))
Bs_abc = trunc(Esym(s, Esym(-s, Bs_ab) + Esym(-s, c)), [a, b, c], 3)
abc = sp.factor(sp.Poly(Bs_abc, a, b, c).coeff_monomial(a*b*c))
pair = sp.factor(sp.Poly(Bs_abc, a, b, c).coeff_monomial(a*b))
three_body = sp.factor(abc - pair**2)
note(f"rung s, three inputs: pair coefficient {pair}, triple-corner coefficient {abc}")
note(f"three-body term = triple corner - (pair coefficient)^2 = {three_body}")
rep("1a rung s: pair coefficient is s and triple corner is s^2 + s(s-1)/2",
    sp.simplify(pair - s) == 0 and sp.simplify(abc - (s**2 + s*(s-1)/2)) == 0)
rep("1b the three-body term s(s-1)/2 vanishes only at rungs 0 and 1",
    set(sp.solve(three_body, s)) == {0, 1})
k = sp.symbols("k")
W = lambda p, q_: p + q_ + k*p*q_
wabc = sp.Poly(sp.expand(W(W(a, b), c)), a, b, c).coeff_monomial(a*b*c)
rep("1c weighted mult-add a+b+k ab has triple corner k^2: no three-body term for any k",
    sp.simplify(wabc - k**2) == 0)
# the KdV 3-soliton tau has triple corner A12 A13 A23 (product of the pairs)
A12, A13, A23 = sp.symbols("A12 A13 A23")
tau3 = 1 + a + b + c + A12*a*b + A13*a*c + A23*b*c + A12*A13*A23*a*b*c
tau3k = sp.Poly(sp.expand(tau3.subs({A12: k, A13: k, A23: k})), a, b, c).coeff_monomial(a*b*c)
rep("1d KdV corners multiply over PAIRS: equal weights A give triple corner A^3, not weighted mult-add's A^2",
    sp.simplify(tau3k - k**3) == 0)
# corners of the cube for n inputs: rung s gives n! [y^n] E^s(y); solitons give A^C(n,2)
y = sp.symbols("y")
O5 = 6
E5 = sum(y**i / sp.factorial(i) for i in range(1, O5 + 1))
def comp5(f, g):
    P = sp.Poly(sp.expand(f.subs(y, g)), y)
    return sp.expand(sum(cf * y**m[0] for m, cf in P.terms() if m[0] <= O5))
it5 = [y]
for n in range(1, O5 + 1): it5.append(comp5(E5, it5[-1]))
def beta(n):
    pts = list(range(n)); cf = 0
    for j in pts:
        wl = 1
        for m in pts:
            if m != j: wl *= (s - m) / sp.Integer(j - m)
        cf += it5[j].coeff(y, n) * wl
    return sp.expand(cf)
corner_rung = {n: sp.factor(sp.factorial(n) * beta(n)) for n in range(2, O5 + 1)}
for n, cr in corner_rung.items():
    note(f"n={n}: rung-s corner {cr}   soliton corner A^{n*(n-1)//2}")
common = set([0, 1])
sols3 = set(sp.solve(corner_rung[3] - s**3, s))
sols4 = set(sp.solve(corner_rung[4] - s**6, s))
note(f"rung s = A and three-body corners agree at s in {sorted(sols3)}; four-body at s in {sorted(sols4, key=lambda z: complex(z).real)}")
rep("1e rung flow and soliton dial share 2-body corners (s = A) but part at 3 or 4 bodies unless s = A is 0 or 1",
    (sols3 & sols4 & set(sp.Reals.intersect(sp.Interval(-10, 10)).contains(z) and z for z in sols3)) <= common
    and 0 in sols3 and 1 in sols3 and sp.Rational(1, 2) in sols3 and sp.Rational(1, 2) not in sols4)
RES["corner_rung"] = {n: str(v) for n, v in corner_rung.items()}
RES["three_body_formula"] = str(three_body)

# ---------------------------------------------------------------------------
print("== 2. KdV: the corner rule is forced (Hirota bilinear residual) ==")
def tau_parts(ks, ds, X, T, corner_weight):
    """tau and its x,t derivatives as sums over cube corners S."""
    n = len(ks)
    out = {key: 0.0 for key in ["t", "x", "xx", "xxx", "xxxx", "T", "xT"]}
    for mask in range(1 << n):
        S = [i for i in range(n) if mask >> i & 1]
        kS = sum(ks[i] for i in S); wS = sum(-ks[i]**3 for i in S)
        eta = sum(ks[i]*X - ks[i]**3*T + ds[i] for i in S)
        cS = corner_weight(S)
        e = cS * np.exp(eta)
        out["t"] += e; out["x"] += kS*e; out["xx"] += kS**2*e
        out["xxx"] += kS**3*e; out["xxxx"] += kS**4*e
        out["T"] += wS*e; out["xT"] += kS*wS*e
    return out
def bilinear_residual(p):
    tau = p["t"]
    r = 2*(tau*p["xT"] - p["x"]*p["T"]) + 2*(tau*p["xxxx"] - 4*p["x"]*p["xxx"] + 3*p["xx"]**2)
    return r / tau**2
def Aij(ki, kj): return ((ki - kj)/(ki + kj))**2
ks = [0.8, 1.3, 1.9]; ds = [0.0, 1.0, -2.0]
Xg, Tg = np.meshgrid(np.linspace(-15, 15, 241), np.linspace(-3, 3, 61))
def cw_kdv(S, delta=0.0, pair_scale=1.0):
    w = 1.0
    for i in S:
        for j in S:
            if i < j: w *= Aij(ks[i], ks[j]) * pair_scale
    if len(S) == 3: w += delta
    return w
r0 = np.max(np.abs(bilinear_residual(tau_parts(ks, ds, Xg, Tg, cw_kdv))))
rep("2a three KdV solitons with corner weight = product of pair weights solve KdV", r0 < 1e-10,
    f"max bilinear residual {r0:.1e}")
for delta in [1e-3, 1e-2, 1e-1]:
    r = np.max(np.abs(bilinear_residual(tau_parts(ks, ds, Xg, Tg, lambda S: cw_kdv(S, delta)))))
    note(f"add a three-body term {delta:g} to the triple corner: residual {r:.2e}")
rd = np.max(np.abs(bilinear_residual(tau_parts(ks, ds, Xg, Tg, lambda S: cw_kdv(S, 1e-2)))))
rep("2b any three-body term breaks KdV (residual grows with it)", rd > 1e-6, f"delta=0.01 gives {rd:.1e}")
rs = np.max(np.abs(bilinear_residual(tau_parts(ks[:2], ds[:2], Xg, Tg, lambda S: cw_kdv(S, 0, 1.05)))))
rep("2c the pair weight is pinned too: A = ((k1-k2)/(k1+k2))^2, a 5% change breaks it", rs > 1e-6,
    f"residual {rs:.1e}")
kk = np.linspace(0.05, 5, 400)
AA = Aij(1.0, kk)
rep("2d every KdV collision weight lies in [0,1): 0 when speeds match, -> 1 when one is far taller",
    AA.min() >= 0 and AA.max() < 1, f"range {AA.min():.3f} .. {AA.max():.3f}")
RES["kdv_ks"] = ks

# ---------------------------------------------------------------------------
print("== 3. Tropical picture: the soliton paths are the creases of a max ==")
k1, k2 = 1.0, 1.6
A = Aij(k1, k2); phase = math.log(A)
def corners(X, T):
    e1 = k1*X - k1**3*T; e2 = k2*X - k2**3*T
    return np.stack([0*X, e1, e2, e1 + e2 + phase])
X3, T3 = np.meshgrid(np.linspace(-25, 25, 801), np.linspace(-8, 8, 321))
H = corners(X3, T3)
logtau = np.log(np.exp(H).sum(0))
gap = logtau - H.max(0)
rep("3a log tau = log-sum-exp of the four corner heights; it differs from their max by at most log 4",
    gap.min() >= -1e-12 and gap.max() <= math.log(4) + 1e-12, f"gap range {gap.min():.2e} .. {gap.max():.3f}")
# fast soliton (k2) peak position far before and far after collision vs tropical lines
def peak_x(T, which):
    xs = np.linspace(-80, 80, 160001)
    e1 = k1*xs - k1**3*T; e2 = k2*xs - k2**3*T
    tau = 1 + np.exp(e1) + np.exp(e2) + A*np.exp(e1 + e2)
    u = 2*np.gradient(np.gradient(np.log(tau), xs), xs)
    # the k2 soliton is the taller one: height k2^2/2
    return xs[np.argmax(u)]
xb = peak_x(-25.0, 2); xa = peak_x(25.0, 2)
# tropical prediction: before (slow soliton far right) the fast crease is 0 = e2 -> x = k2^2 T,
# after: crease e1 = e1+e2+phase -> x = k2^2 T - phase/k2
shift_obs = (xa - k2**2*25.0) - (xb - k2**2*(-25.0))
rep("3b the fast soliton's jump equals the tropical crease offset -log(A)/k2",
    abs(shift_obs - (-phase/k2)) < 2e-3, f"observed {shift_obs:.4f}, predicted {-phase/k2:.4f}")
B_bottom = Es(-4, Es(4, np.array([0.4])) + Es(4, np.array([0.7])))
note(f"rung sheet bottom check: B_-4(0.4, 0.7) = {B_bottom[0]:.6f}  (max = 0.7)")
rep("3c far down the rungs (s = -4) the rung law is already nearly the max", abs(B_bottom[0] - 0.7) < 0.01)
RES["tropical"] = {"k1": k1, "k2": k2, "A": A, "phase": phase, "shift": shift_obs}

# ---------------------------------------------------------------------------
print("== 4. The gluing number q: cnoidal waves ==")
mp.mp.dps = 40
def Th(w, q, N=None):
    if N is None:
        N = int(math.sqrt(2*120/max(-math.log(float(q)), 1e-9))) + 6
    return mp.fsum(mp.power(q, n*(n - 1)//2) * mp.power(w, n) for n in range(-N, N + 1))
def triple(w, q, M=400):
    out = mp.mpf(1)
    for m in range(1, M):
        out *= (1 - q**m) * (1 + w*q**(m - 1)) * (1 + q**m / w)
    return out
errs = []
for q in [mp.mpf("0.05"), mp.mpf("0.3"), mp.mpf("0.6")]:
    for w in [mp.mpf("0.4"), mp.mpf("1.7"), mp.mpf("-0.37")]:
        errs.append(abs(Th(w, q) - triple(w, q)) / abs(Th(w, q)))
rep("4a Jacobi triple product: Th(w;q) = (q;q) (1+w) prod_m (1+w q^m)(1+q^m/w)", max(errs) < 1e-25,
    f"max rel err {float(max(errs)):.1e}")
# so Th/(q;q) - 1 = mult-add of the whole q-orbit {w q^m} U {q^m / w}
q = mp.mpf("0.3"); w = mp.mpf("0.8")
orbit = [w*q**m for m in range(0, 300)] + [q**m/w for m in range(1, 300)]
ma = mp.mpf(0)
for y in orbit: ma = ma + y + ma*y
qq = mp.mpf(1)
for m in range(1, 400): qq *= (1 - q**m)
e = abs((Th(w, q)/qq - 1) - ma)
rep("4b Th/(q;q)_inf - 1 is the mult-add of every point of w's q-orbit", e < 1e-25, f"err {float(e):.1e}")
def g_theta(eta, q):  # 2 d^2/deta^2 log Th(e^eta; q)
    return 2*mp.diff(lambda z: mp.log(Th(mp.e**z, q)), eta, 2)
def g_sum(eta, q, M=200):
    L = -mp.log(q)
    return mp.fsum(mp.mpf(1)/2 * mp.sech((eta - m*L)/2)**2 for m in range(-M, M + 1))
errs = []
for qv in ["0.01", "0.2", "0.5", "0.8"]:
    qm = mp.mpf(qv)
    for eta in [-0.7, 0.0, 0.9, 2.3]:
        errs.append(abs(g_theta(eta, qm) - g_sum(eta, qm)))
rep("4c cnoidal profile = exact sum of copies of the q=0 soliton, spaced by -log q (Toda; Boyd 1984)",
    max(errs) < 1e-15, f"max err {float(max(errs)):.1e}")
# KdV: g'' + 3 g^2 + alpha g + beta = 0 (travelling wave, any constant offset)
kdv_fit = {}
for qv in ["1e-10", "0.02", "0.2", "0.5", "0.8"]:
    qm = mp.mpf(qv); L = float(-mp.log(qm))
    etas = [mp.mpf(L)*j/9 - mp.mpf(L)/2 for j in range(10)]
    G = [g_sum(e_, qm) for e_ in etas]
    G2 = [mp.diff(lambda z: g_sum(z, qm), e_, 2) for e_ in etas]
    Mtx = mp.matrix([[G[i], 1] for i in range(10)])
    rhs = mp.matrix([-(G2[i] + 3*G[i]**2) for i in range(10)])
    sol = mp.lu_solve(Mtx.T*Mtx, Mtx.T*rhs)
    resid = max(abs(G2[i] + 3*G[i]**2 + sol[0]*G[i] + sol[1]) for i in range(10))
    kdv_fit[qv] = {"alpha": float(sol[0]), "beta": float(sol[1]), "resid": float(resid)}
    note(f"q={qv}: alpha={float(sol[0]):+.6f} beta={float(sol[1]):+.3e} residual {float(resid):.1e}")
rep("4d every q gives an exact KdV travelling wave (constant offset and speed set by q)",
    max(v["resid"] for v in kdv_fit.values()) < 1e-15)
rep("4e at q -> 0 the offset and speed shift vanish: alpha -> -1, the lone soliton's speed",
    abs(kdv_fit["1e-10"]["alpha"] + 1) < 1e-3, f"alpha(1e-10) = {kdv_fit['1e-10']['alpha']:.6f}")
RES["kdv_fit"] = kdv_fit
# harmonic ratio c2/c1 = sech(2 pi^2 / L)
hr = {}
for qv in [1e-6, 1e-3, 0.01, 0.1, 0.3]:
    L = -math.log(qv); Npts = 4096
    etas = np.arange(Npts) * L / Npts
    M = 60
    gvals = sum(0.5/np.cosh((etas - m*L)/2)**2 for m in range(-M, M + 1))
    F = np.fft.rfft(gvals)
    r = abs(F[2])/abs(F[1]); pred = 1/math.cosh(2*math.pi**2/L)
    hr[qv] = (r, pred)
ok = all(abs(r - p) < 1e-7 + 1e-6*p for r, p in hr.values())
rep("4f harmonic ratio c2/c1 = sech(2 pi^2 / (-log q)): 1 at q=0 (all harmonics: a soliton), "
    "-> 0 at q=1 (one sinusoid: a linear wave)", ok,
    ", ".join(f"q={k_}: {v[0]:.3g}" for k_, v in hr.items()))
RES["harmonic_ratio"] = {str(k_): v[0] for k_, v in hr.items()}
# ripple: as q->1 the wave is mean + tiny cosine, log Th ~ linear superposition
for qv in [0.9, 0.97]:
    L = -math.log(qv); etas = np.linspace(0, L, 2001)
    gvals = sum(0.5/np.cosh((etas - m*L)/2)**2 for m in range(-400, 401))
    note(f"q={qv}: mean {gvals.mean():.4f}, ripple height {gvals.max()-gvals.min():.2e}")
# gluing is a rung-1 step: w q - 1 = (w-1) (+) (q-1),  and q - 1 = E(log q)
wq = 0.73; qv = 0.21
lhs = wq*qv - 1; rhs_ = (wq - 1) + (qv - 1) + (wq - 1)*(qv - 1)
rep("4g gluing w -> w q is mult-add by q-1 in the 1+x picture, and q-1 = E(log q)",
    abs(lhs - rhs_) < 1e-15 and abs(np.expm1(np.log(qv)) - (qv - 1)) < 1e-15)

# ---------------------------------------------------------------------------
print("== 5. Two waves: the soliton cube is the q=0 face of the theta lattice ==")
def Th2(w1, w2, q1, q2, Aw, N=12):
    tot = mp.mpf(0)
    for n1 in range(-N, N + 1):
        for n2 in range(-N, N + 1):
            tot += (q1**(n1*(n1 - 1)//2) * q2**(n2*(n2 - 1)//2) * mp.power(Aw, n1*n2)
                    * mp.power(w1, n1) * mp.power(w2, n2))
    return tot
w1, w2, Aw = mp.mpf("0.6"), mp.mpf("1.4"), mp.mpf(Aij(1.0, 1.6))
cube = 1 + w1 + w2 + Aw*w1*w2
face = []
for qv in ["1e-3", "1e-4", "1e-5", "1e-6"]:
    qm = mp.mpf(qv)
    face.append((qv, float(abs(Th2(w1, w2, qm, qm, Aw) - cube) / qm)))
note("error / q as q -> 0: " + ", ".join(f"{a_}: {b_:.4f}" for a_, b_ in face))
rep("5a genus-2 theta -> 1 + w1 + w2 + A w1 w2 (the weighted mult-add cube), error proportional to q",
    abs(face[-1][1] - face[-2][1]) / face[-1][1] < 0.01)
rep("5b with A = 1 the two-wave theta splits into two one-wave thetas (no interaction)",
    abs(Th2(w1, w2, mp.mpf("0.3"), mp.mpf("0.4"), 1) - Th(w1, mp.mpf("0.3"))*Th(w2, mp.mpf("0.4"))) < 1e-25)
RES["face"] = face

# ---------------------------------------------------------------------------
print("== 6. The rung flow on a soliton: p_t = 2 d^2/deta^2 E^-t(e^eta) ==")
etas = np.linspace(-14, 22, 1801); h = etas[1] - etas[0]
flow = {}
tgrid = [round(v, 2) for v in np.arange(0.0, 2.51, 0.05)]
for tv in tgrid:
    F = Es(-tv, np.exp(etas)) if tv != 0 else np.exp(etas)
    p = 2*np.gradient(np.gradient(F, h), h)
    p2 = np.gradient(np.gradient(p, h), h)
    # solitary-wave test on eta in [-10, 10]: p'' + 3 p^2 - V p = C
    m_ = (etas > -10) & (etas < 10)
    Mt = np.stack([p[m_], np.ones(m_.sum())], 1)
    rhs = -(p2[m_] + 3*p[m_]**2)
    sol, *_ = np.linalg.lstsq(Mt, rhs, rcond=None)
    res = np.max(np.abs(Mt @ sol - rhs)) / max(np.max(np.abs(p2[m_])), 1e-12)
    i10 = np.argmin(abs(etas - 10)); i20 = np.argmin(abs(etas - 20))
    flow[tv] = {"V": float(-sol[0]), "resid": float(res), "p10": float(p[i10]), "p20": float(p[i20])}
best = min(flow, key=lambda tv: flow[tv]["resid"])
for tv in [0.0, 0.5, 0.9, 1.0, 1.1, 1.5, 2.0]:
    f = flow[tv]
    note(f"t={tv:.2f}: KdV mismatch {f['resid']:.1e}, tail p(10)={f['p10']:+.2e}, p(20)={f['p20']:+.2e}")
rep("6a only a whole rung step t = 1 turns the exponential into a KdV soliton", best == 1.0 and
    flow[1.0]["resid"] < 1e-3 and min(flow[0.9]["resid"], flow[1.1]["resid"]) > 10*flow[1.0]["resid"],
    f"best t = {best}, mismatch {flow[1.0]['resid']:.1e} vs {flow[0.9]['resid']:.1e} at 0.9")
rep("6b at t = 1 the wave speed is the lone soliton's (V = 1)", abs(flow[1.0]["V"] - 1) < 1e-3,
    f"V = {flow[1.0]['V']:.5f}")
rep("6c below one rung (t<1) the right tail grows: no lump", flow[0.5]["p20"] > flow[0.5]["p10"] > 0)
rep("6d above one rung (t>1) the tail turns negative and decays only slowly (algebraic)",
    flow[2.0]["p20"] < 0 and abs(flow[2.0]["p20"]) > 1e-4)
RES["flow"] = {str(k_): v for k_, v in flow.items()}

npass = sum(1 for v in RES.values() if isinstance(v, dict) and "pass" in v and v["pass"])
nall = sum(1 for v in RES.values() if isinstance(v, dict) and "pass" in v)
print(f"\n{npass} of {nall} checks pass"); OUT.append(f"\n{npass} of {nall} checks pass")
open("/mnt/project-files/math-theory/solitons-rungs/checks_output.txt", "w").write("\n".join(OUT) + "\n")
json.dump(RES, open("/mnt/project-files/math-theory/solitons-rungs/results.json", "w"), indent=1, default=str)
