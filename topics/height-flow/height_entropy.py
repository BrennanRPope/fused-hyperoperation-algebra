"""Entropies built from the hypernatural logarithm, tested.

ln_s(X) = E^-s(X - 1) = E^(1-s)(log X)      (E(x) = e^x - 1, E^r its regular iteration)
So ln_s is a "group logarithm" in Tempesta's sense with G = E^r, r = 1 - s, and the
formal group law G(G^-1 x + G^-1 y) = B_r(x, y), the rung-r operation.

Forms tested (r = 1 - s throughout):
  T_s(p)       = sum p_i ln_s(1/p_i)                 trace form (Tempesta's universal-group form)
  H_{s,a}(p)   = ln_s(D_a(p)) = E^r(R_a(p))          hypernatural log of the Hill number
                                                     (= Jensen-Tempesta's non-trace form, eq. 14)
  Z_{s,a}(p)   = ln_s(sum p^a)/(1-a)                 Tempesta's Z-entropy with this G
  Esc_{s,q}(p) = sum P_i ln_s(1/p_i), P = p^q/sum p^q   escort form
Writes height_entropy.json; prints PASS/FAIL lines.
"""
import sys, os, json, math
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "fractional-rungs"))
from regular_iteration import Es, B, coeffs_at   # checked 16/16 in that folder

res = {}
def rep(name, ok, extra=""):
    print(("PASS " if ok else "FAIL ") + name + (("   " + extra) if extra else ""))
    res[name] = bool(ok)

def G(r, t):
    t = np.asarray(t, float)
    if r == 0: return t.copy()
    if r == 1: return np.expm1(t)
    if r == -1: return np.log1p(t)
    return Es(r, t)

def ln_s(s, X): return G(1 - s, np.log(np.asarray(X, float)))

def renyi(p, a):
    p = p[p > 0]
    if abs(a - 1) < 1e-12: return float(-(p*np.log(p)).sum())
    return float(np.log((p**a).sum())/(1 - a))

def T(s, p):
    p = p[p > 0]; return float((p*ln_s(s, 1/p)).sum())
def H(s, p, a=1.0): return float(G(1 - s, np.array([renyi(p, a)]))[0])
def Z(s, p, a):
    p = p[p > 0]; return float(ln_s(s, np.array([(p**a).sum()]))[0]/(1 - a))
def Esc(s, p, q):
    p = p[p > 0]; P = p**q/(p**q).sum(); return float((P*ln_s(s, 1/p)).sum())

rng = np.random.default_rng(7)
def rand_p(n):
    x = rng.exponential(1, n); return x/x.sum()

S_LIST = [-0.5, 0.25, 0.5, 0.75, 1.5, 2.0]

# ---------- 1. the logarithm itself ----------
Xs = np.array([0.3, 0.8, 1.7, 3.0, 9.0])
worst = 0
for s in S_LIST:
    worst = max(worst, np.max(np.abs(ln_s(s, Xs) - Es(-s, Xs - 1))))
rep("ln_s(X) = E^(1-s)(log X) equals E^-s(X-1)", worst < 1e-8, f"worst {worst:.1e}")
worst = 0
for s in S_LIST:
    X, Y = np.array([1.3, 2.2, 4.0]), np.array([1.9, 1.1, 2.5])
    worst = max(worst, np.max(np.abs(ln_s(s, X*Y) - B(1 - s, ln_s(s, X), ln_s(s, Y)))))
rep("ln_s(XY) = B_(1-s)(ln_s X, ln_s Y)  (group logarithm with law B_(1-s))", worst < 1e-8, f"worst {worst:.1e}")

# ---------- 2. landmarks: which heights are already known entropies ----------
p = rand_p(7)
sh = renyi(p, 1)
rep("T_1 = Shannon", abs(T(1, p) - sh) < 1e-12)
rep("T_0 = (support size) - 1  = Tsallis q = 0", abs(T(0, p) - (len(p) - 1)) < 1e-10)
ok = all(abs(H(1, p, a) - renyi(p, a)) < 1e-12 for a in (0.5, 2, 3))
rep("H_{1,a} = Renyi R_a", ok)
ok = all(abs(H(0, p, a) - (math.exp(renyi(p, a)) - 1)) < 1e-10 for a in (0.5, 1, 2))
sm = lambda p, a, b: (math.exp((1 - b)*renyi(p, a)) - 1)/(1 - b)
ok2 = all(abs(H(0, p, a) - sm(p, a, 0.0)) < 1e-10 for a in (0.5, 1, 2))
rep("H_{0,a} = D_a - 1 = Sharma-Mittal with beta = 0", ok and ok2)
# Tempesta-Jensen (arXiv 1912.10907, eq. 21-22) "rapidly growing" W = exp(k2 e^{k3 N} - k1): Z = k1 ln(ln(sum p^a)/(k1(1-a)) + 1)
tj = lambda p, a, k1: k1*math.log(math.log((p**a).sum())/(k1*(1 - a)) + 1)
ok = all(abs(H(2, p, a) - tj(p, a, 1.0)) < 1e-10 for a in (0.5, 2, 3))
rep("H_{2,a} = log(1 + R_a) = Tempesta-Jensen 2019 'rapidly growing' entropy with k1 = 1 (known)", ok)

# ---------- 3. composition for independent systems ----------
comp = {}
for form in ("H", "Z", "T", "Esc"):
    comp[form] = {}
for s in [-0.5, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 0.0]:
    eH = eZ = eT = eE = 0.0; eTu = 0.0
    for trial in range(20):
        pa, pb = rand_p(rng.integers(2, 6)), rand_p(rng.integers(2, 6))
        pab = np.outer(pa, pb).ravel()
        r = 1 - s
        for a in (0.5, 1.0, 2.0):
            lhs = H(s, pab, a); rhs = float(B(r, H(s, pa, a), H(s, pb, a)))
            eH = max(eH, abs(lhs - rhs)/max(1, abs(lhs)))
        for a in (0.5, 2.0):
            za, zb = Z(s, pa, a), Z(s, pb, a)
            rhs = float(B(r, (1 - a)*za, (1 - a)*zb))/(1 - a)
            eZ = max(eZ, abs(Z(s, pab, a) - rhs)/max(1, abs(rhs)))
        lhs = T(s, pab); rhs = float(B(r, T(s, pa), T(s, pb)))
        eT = max(eT, abs(lhs - rhs)/max(1, abs(lhs)))
        lhs = Esc(s, pab, 2.0); rhs = float(B(r, Esc(s, pa, 2.0), Esc(s, pb, 2.0)))
        eE = max(eE, abs(lhs - rhs)/max(1, abs(lhs)))
        na, nb = rng.integers(2, 6), rng.integers(2, 6)
        ua, ub = np.ones(na)/na, np.ones(nb)/nb
        eTu = max(eTu, abs(T(s, np.outer(ua, ub).ravel()) - float(B(r, T(s, ua), T(s, ub)))))
    comp["H"][s] = eH; comp["Z"][s] = eZ; comp["T"][s] = eT; comp["Esc"][s] = eE
    comp.setdefault("T_uniform", {})[s] = eTu
print("worst relative composition error by height:")
for f in comp: print("  ", f, {k: f"{v:.1e}" for k, v in comp[f].items()})
rep("Hill form H_{s,a} composes EXACTLY by B_(1-s), all distributions, a = 0.5, 1, 2, every height",
    max(comp["H"].values()) < 1e-7)
rep("Tempesta Z form composes exactly by its own law (Jensen-Tempesta eq. 23)", max(comp["Z"].values()) < 1e-7)
rep("trace form T_s composes exactly on uniform distributions (weak composability)", max(comp["T_uniform"].values()) < 1e-7)
fracT = [comp["T"][s] for s in comp["T"] if s not in (0.0, 1.0)]
rep("trace form T_s FAILS exact composition off uniform at every height except s = 0, 1 (Tempesta: only Tsallis is strongly composable)",
    min(fracT) > 1e-3 and comp["T"][0.0] < 1e-8 and comp["T"][1.0] < 1e-8,
    f"smallest failure {min(fracT):.2e}")
rep("escort form composes exactly only at s = 0 and 1, fails elsewhere",
    min(comp["Esc"][s] for s in comp["Esc"] if s not in (0.0, 1.0)) > 1e-3 and comp["Esc"][0.0] < 1e-8 and comp["Esc"][1.0] < 1e-8)
res["composition_errors"] = {f: {str(k): v for k, v in d.items()} for f, d in comp.items()}

# ---------- 4. what only this family has: heights compose ----------
worst = 0
for (s, t) in [(0.5, 0.25), (1.0, 0.5), (1.5, -0.5), (0.75, 0.75), (2.0, 1.0)]:
    for a in (0.5, 1.0, 2.0):
        q = rand_p(6)
        worst = max(worst, abs(H(s - t, q, a) - float(G(t, np.array([H(s, q, a)]))[0])))
rep("height shift: H_{s-t} = E^t(H_s)  (one flow; Tsallis/Sharma-Mittal dial has no such law)", worst < 1e-7, f"worst {worst:.1e}")
# Sharma-Mittal check: G_k o G_j != G_{k+j}
Gk = lambda k, t: np.expm1(k*t)/k
d = abs(Gk(0.5, Gk(0.5, 1.0)) - Gk(1.0, 1.0))
rep("dial family does NOT compose: G_0.5(G_0.5(1)) != G_1(1) for (e^{kt}-1)/k", d > 1e-2, f"gap {d:.3f}")
worst = 0
for s in [-0.5, 0.0, 0.5, 1.0, 1.5]:
    for W in [2.0, 3.5, 10.0, 40.0]:
        worst = max(worst, abs(float(ln_s(s + 1, np.array([math.exp(W - 1)]))[0]) - float(ln_s(s, np.array([W]))[0])))
rep("power-set step: ln_(s+1)(e^(W-1)) = ln_s(W)  (one height up undoes one 'exponentiate the state count')", worst < 1e-7, f"worst {worst:.1e}")

# ---------- 5. comparison with Sharma-Mittal and Tsallis at third order ----------
fit = {}
for r in [0.25, 0.5, 0.75]:
    Rs = np.array([0.02, 0.04, 0.06])
    diff = G(r, Rs) - Gk(r, Rs)
    c3 = diff/Rs**3
    fit[r] = (float(c3[0]), r*(r - 1)/12)
okf = all(abs(v[0] - v[1]) < 0.01*max(1, abs(v[1])) + 2e-3 for v in fit.values())
rep("H_s minus Sharma-Mittal (matched k = r) = r(r-1)/12 * R^3 + ...: agree to 2nd order, part at 3rd", okf,
    json.dumps({str(k): [round(a, 4), round(b, 4)] for k, (a, b) in fit.items()}))
cs = coeffs_at(0.5)
rep("E^r has a quadratic term (r/2): not Kaniadakis (sinh, odd) unless r = 0", abs(cs[2] - 0.25) < 1e-12)
# Non-analytic at fractional height (Baker): coefficient growth of formal E^(1/2)
# (computed in ../fractional-rungs/regular_iteration_output.txt, 'formal half-rung series diverges')
txt = open(os.path.join(HERE, "..", "fractional-rungs", "regular_iteration_output.txt")).read()
rep("fractional-height G is not analytic at 0 (Baker; reused check), so it is none of Tempesta's listed entire G's",
    "PASS formal half-rung series diverges" in txt)

# ---------- 6. Shannon-Khinchin-type axioms ----------
# concavity: entropies only use X >= 1 (X = 1/p_i or a Hill number), so test there.
def concave_on(f, x):
    y = f(x)
    mid = y[1:-1] - 0.5*(y[:-2] + y[2:])          # >= 0 for concave on an evenly spaced grid
    return bool(np.all(mid >= -1e-11*np.maximum(1, np.abs(y[1:-1]))))
Xg = np.linspace(1.0, 400.0, 4001)
pg = np.linspace(1e-6, 1.0, 4001)
conc = {}
for s in [-1.0, -0.5, 0.0, 0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0]:
    lncon = concave_on(lambda X: ln_s(s, X), Xg)
    gcon = concave_on(lambda q: q*ln_s(s, 1/q), pg)
    bend = float((ln_s(s, np.array([1 + 1e-3]))[0] + ln_s(s, np.array([1 - 1e-3]))[0])/1e-6)
    g0 = float(1e-300*ln_s(s, np.array([1e300]))[0]) if s > 0 else (1.0 if s == 0 else float("inf"))
    conc[s] = {"ln_s concave on X>=1": lncon, "trace g concave": gcon, "bend at 1": round(bend, 4), "g(0)": g0}
print("concavity / expansibility table:")
for s, v in conc.items(): print(f"   s={s:5}: ", v)
rep("ln_s concave on X >= 1 for every s >= 0 (Tempesta's group-log condition where entropy uses it); bend at 1 is -s, so convex for s < 0",
    all(conc[s]["ln_s concave on X>=1"] for s in conc if s >= 0) and not any(conc[s]["ln_s concave on X>=1"] for s in conc if s < 0)
    and all(abs(conc[s]["bend at 1"] + s) < 1e-3 for s in conc))
rep("trace form g(p) = p ln_s(1/p) concave on [0,1] for every s >= 0 tested (0 to 3), not for s < 0",
    all(conc[s]["trace g concave"] for s in conc if s >= 0) and not any(conc[s]["trace g concave"] for s in conc if s < 0))
rep("expansibility g(0) = 0 for s > 0; fails at s = 0 (counts W-1) and blows up for s < 0",
    all(conc[s]["g(0)"] < 1e-6 for s in conc if s > 0))
res["concavity"] = {str(k): v for k, v in conc.items()}

# ---------- 7. Lesche stability (numerical evidence on two standard families, not a proof) ----------
# For huge W we never build arrays: entropies of the two families are written in closed form, and
# ratios of G use the exact identity log(1 + E^r(x)) = E^(r-1)(x).
def logG1(r, x): return float(G(r - 1, np.array([x]))[0])        # log(1 + G(x))
def ratio_G(r, a, b):   # G(a)/G(b) for large arguments, a, b > 0
    ga = float(G(r, np.array([a]))[0]); gb = float(G(r, np.array([b]))[0])
    if np.isfinite(ga) and np.isfinite(gb) and gb < 1e300: return ga/gb
    return math.exp(logG1(r, a) - logG1(r, b))
def lesche(form, s, delta, Ls):
    r = 1 - s; out = []
    for L in Ls:
        if form == "T":
            # A: point mass vs smeared (1-d+d/W, d/W, ...): S(p) = 0, S(p') ~ d*G(L + ln(1/d)) + top term
            top = 1 - delta
            gtop = top*float(G(r, np.array([-math.log(top)]))[0])
            dA = (gtop + delta*ratio_G(r, L - math.log(delta), L)*float(G(r, np.array([L]))[0])) if L < 600 else None
            rA = (gtop/float(G(r, np.array([L]))[0]) if L < 600 else 0.0) + delta*ratio_G(r, L - math.log(delta), L)
            # B: uniform vs uniform + spike d: rest of mass (1-d) sits at ln(1/p) = L - ln(1-d)
            rB = abs((1 - delta)*ratio_G(r, L - math.log(1 - delta), L) - 1)
            out.append(max(rA, rB))
        else:  # Hill form, a = 1: E^r(Shannon)
            hA = -(1 - delta)*math.log(1 - delta) + delta*(L - math.log(delta))     # W -> infinity form
            hB = -delta*math.log(delta) - (1 - delta)*(math.log(1 - delta) - L)     # leading terms
            rA = ratio_G(r, hA, L)
            rB = abs(1 - ratio_G(r, hB, L))
            out.append(max(rA, rB))
    return out
Ls = np.geomspace(3, 1e5, 40)
les = {}
DELTAS = (1e-2, 1e-3, 1e-4, 1e-5)
for form in ("T", "H"):
    for s in [0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0]:
        les[(form, s)] = {d: lesche(form, s, d, Ls) for d in DELTAS}
print("Lesche: sup over ln W in [3, 1e5] of |dS|/S_max for delta = 1e-2, 1e-3, 1e-4, 1e-5, and the value at ln W = 1e5 (delta = 1e-2):")
stab = {}
for (form, s), v in les.items():
    sups = [max(v[d]) for d in DELTAS]; end2 = v[1e-2][-1]
    stab[(form, s)] = (sups, end2)
    print(f"   {form} s={s:4}: " + " ".join(f"{x:.4f}" for x in sups) + f"   (at ln W = 1e5: {end2:.3f})")
rep("trace form: the sup shrinks steadily as delta -> 0 for every s > 0 (evidence of Lesche stability)",
    all(stab[("T", s)][0][3] < 0.35*stab[("T", s)][0][0] and all(np.diff(stab[("T", s)][0]) < 0) for s in [0.25, 0.5, 0.75, 1.0, 1.5, 2.0]))
rep("trace form unstable at s = 0 (Tsallis q = 0 edge): sup stays at 1", min(stab[("T", 0.0)][0]) > 0.5)
# Hill form: the analytic reason (ours): for s < 1, G = E^r grows faster than any power, so G((1-d)L)/G(L) -> 0
# (family B, ratio -> 1); for s > 1, G grows slower than any power, so G(d L)/G(L) -> 1 (family A, ratio -> 1).
# Only s = 1 (Shannon) escapes both. The climb is extremely slow, so we follow it out to ln W = 1e12.
LsH = np.geomspace(3, 1e12, 45)
lesH = {s: lesche("H", s, 1e-2, LsH) for s in [0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0]}
print("Hill form, delta = 1e-2, |dS|/S_max at ln W = 1e2, 1e4, 1e6, 1e8, 1e10, 1e12:")
idx = [int(np.argmin(np.abs(LsH - 10.0**k))) for k in (2, 4, 6, 8, 10, 12)]
for s, v in lesH.items(): print(f"   s={s:4}: " + " ".join(f"{v[i]:.3f}" for i in idx))
climb = {s: lesH[s][idx[-1]] - lesH[s][idx[1]] for s in lesH}
rep("Hill form (a = 1) is Lesche-unstable at every height except s = 1: its error ratio climbs with W at fixed delta for s != 1, falls at s = 1",
    all(climb[s] > 0 for s in lesH if s not in (0.0, 1.0)) and lesH[0.0][idx[-1]] > 0.99 and climb[1.0] < 0,
    json.dumps({str(k): round(v, 3) for k, v in climb.items()}))
res["lesche_hill_long"] = {str(k): v for k, v in lesH.items()}; res["lesche_hill_L"] = LsH.tolist()
res["lesche"] = {f"{f}|{s}": {str(d): v for d, v in dd.items()} for (f, s), dd in les.items()}
res["lesche_L"] = Ls.tolist()

# ---------- 8. maximum-entropy distribution of the trace form ----------
def Gp(r, t, h=1e-5):
    return (G(r, t + h) - G(r, np.maximum(t - h, 0) if np.ndim(t) else max(t - h, 0)))/(t + h - np.maximum(t - h, 0))
Lgrid = np.linspace(0, 80, 8001)
maxent = {}
for s in [0.25, 0.5, 0.75, 1.0, 1.5]:
    r = 1 - s
    Gv = G(r, Lgrid); Gd = np.gradient(Gv, Lgrid)
    Lam = Gv - Gd
    mono = bool(np.all(np.diff(Lam[1:]) > 0))
    eps = np.linspace(0, 400, 40001)
    def probs(a):
        target = a + eps
        Lval = np.interp(target, Lam, Lgrid, left=np.nan, right=np.nan)
        return np.exp(-Lval)
    lo, hi = Lam[1] - 1.0, Lam[1] + 50
    # p_i = exp(-Lam^-1(a + eps_i)); need a >= Lam(0); sum decreasing in a
    lo = Lam[0]
    for _ in range(200):
        mid = (lo + hi)/2
        tot = np.nansum(probs(mid))*(eps[1] - eps[0])
        if tot > 1: lo = mid
        else: hi = mid
    pe = probs((lo + hi)/2)
    maxent[s] = {"monotone": mono, "eps": eps[::100].tolist(), "p": np.nan_to_num(pe[::100]).tolist()}
    # local maximality: perturb within the constraints (density version, discretised)
rep("maxent condition G(L) - G'(L) = a + b*eps is solvable (left side increasing) for s = 0.25 ... 1.5",
    all(v["monotone"] for v in maxent.values()))
# tail class: -ln p vs E^{-r}(eps)
tails = {}
for s in [0.25, 0.5, 0.75]:
    e = np.array(maxent[s]["eps"]); pp = np.array(maxent[s]["p"])
    m = (e > 50) & (pp > 0)
    mlp = -np.log(pp[m])
    slope = np.polyfit(np.log(e[m]), np.log(mlp), 1)[0]   # 1 = exponential, 0 = power law
    tails[s] = float(slope)
print("maxent tail: d log(-log p)/d log eps on eps in [50, 400] (1 = exponential, c = stretched exp, -> 0 power law):",
      {k: round(v, 3) for k, v in tails.items()})
rep("maxent tails at s in (0,1) lie between power law and exponential (local slope strictly between 0 and 1)",
    all(0.02 < v < 0.98 for v in tails.values()))
res["maxent"] = {str(k): v for k, v in maxent.items()}
res["maxent_tail_slopes"] = {str(k): v for k, v in tails.items()}

# ---------- 9. extensivity: which state spaces make S_s proportional to N ----------
ext = {}
for s in [-0.5, 0.0, 0.5, 1.0, 1.5]:
    r = 1 - s
    Ns = np.array([2.0, 5.0, 10.0, 20.0])
    lnW = G(-r, Ns)                     # ln W(N) = E^{-r}(N)
    S = G(r, lnW)                       # uniform: S = ln_s(W) = E^r(ln W)
    ext[s] = float(np.max(np.abs(S - Ns)))
rep("uniform S_s = N exactly on state spaces W(N) = exp(E^(s-1)(N)) (Jensen-Tempesta procedure gives G = E^(1-s))",
    max(ext.values()) < 1e-6, json.dumps({str(k): f"{v:.1e}" for k, v in ext.items()}))
bf = []
for N in [4, 8, 16, 32, 48]:
    lnW = (2.0**N)*math.log(2)          # Boolean functions of N bits: W = 2^(2^N)
    bf.append((N, math.log1p(lnW)/N, lnW/N))
print("Boolean functions on N bits: S_2/N and S_1/N:", [(n, round(a, 4), f"{b:.2e}") for n, a, b in bf])
rep("Boolean functions on N bits (W = 2^(2^N)): S_2/N -> ln 2 (extensive) while Shannon S_1/N blows up",
    abs(bf[-1][1] - math.log(2)) < 0.02 and bf[-1][2] > 1e10)

# ---------- 10. the entropy family as a running coupling (joins the running-coupling handoff) ----------
from fractions import Fraction as Fr
RC = json.load(open(os.path.join(HERE, "..", "height-running-coupling", "running_coupling.json")))
VC = {int(k): float(Fr(v)) for k, v in RC["v_exact"].items()}
def v_small(y):
    return sum(c*y**k for k, c in VC.items() if k <= 8)
def v_of(x, n=4000):
    """Julia's iterative logarithm v(x) = d/dt E^t(x) at t = 0, via v(E(y)) = E'(y) v(y) (independent of Es)."""
    y = x; logprod = 0.0
    for _ in range(n): y = math.log1p(y)
    z = y; vy = v_small(z)
    for _ in range(n):
        logprod += z            # E'(z) = e^z
        z = math.expm1(z)
    return vy*math.exp(logprod)
def Phi(x, n=200000):
    """Abel function (Fatou coordinate) of E on x > 0: Phi(E(x)) = Phi(x) + 1, normalised by Phi(y) ~ -2/y + (1/3) ln y."""
    y = x
    for _ in range(n): y = math.log1p(y)
    return -2/y + math.log(y)/3 + y/36 + n      # y/36: next term of the Fatou expansion (checked below)
# check the Fatou normalisation term numerically: Phi(E(x)) - Phi(x) = 1
e_phi = max(abs(Phi(math.expm1(x)) - Phi(x) - 1) for x in (0.05, 0.3, 1.0))
rep("Fatou coordinate Phi(E(x)) = Phi(x) + 1 (rung ruler, used as the RG-invariant scale)", e_phi < 1e-6, f"err {e_phi:.1e}")
# dH/ds = -v(H): moving along height is a running coupling with beta = v
worst = 0
for R in (0.4, 1.3, 2.5):
    for s in (0.5, 1.5, 2.5):
        h = 1e-4
        dH = (float(G(1 - s - h, np.array([R]))[0]) - float(G(1 - s + h, np.array([R]))[0]))/(2*h)   # d/ds H = -d/dr
        Hs = float(G(1 - s, np.array([R]))[0])
        worst = max(worst, abs(dH + v_of(Hs))/v_of(Hs))
rep("dH_s/ds = -beta(H_s) with beta = Julia's iterative log v(x) = x^2/2 - x^3/12 + ... (the entropy family IS a running-coupling trajectory)",
    worst < 1e-5, f"worst rel {worst:.1e}")
# invariant: Phi(H_s) + s is the same at every height; it is fixed by the Renyi (or Shannon) value alone
inv = {}
for R in (0.4, 1.3, 2.5):
    vals = [Phi(float(G(1 - s, np.array([R]))[0])) + s for s in (0.5, 1.0, 2.0, 4.0, 8.0)]
    inv[R] = max(vals) - min(vals)
rep("RG invariant: Phi(H_s) + s is constant along the family (the 'Lambda' of an entropy is Phi(R) + 1)", max(inv.values()) < 1e-5,
    json.dumps({str(k): f"{v:.1e}" for k, v in inv.items()}))
# universality above Shannon: H_s -> 2/s for every starting entropy, two-loop correction from rho = 1/3
uni = {}
for R in (0.4, 1.3, 2.5, 10.0):
    S = 60.0
    Hs = float(G(1 - S, np.array([R]))[0])
    c = Phi(R) + 1                           # invariant
    t = S - c                                # Phi(H) = c - S  =>  H solves -2/H + (1/3) ln H = c - S
    one = 2/t
    two = 2/(t + (1/3)*math.log(2/t))        # one step of the two-loop correction
    uni[R] = (Hs, one, two)
rep("far above Shannon every entropy runs to 0 like 2/s (one loop), and the two-loop log term closes most of the gap",
    all(abs(v[2] - v[0]) < 0.1*abs(v[1] - v[0]) for v in uni.values()),
    json.dumps({str(k): [f"{a:.6f}", f"{b:.6f}", f"{c:.6f}"] for k, (a, b, c) in uni.items()}))
res["rg_invariant_spread"] = {str(k): v for k, v in inv.items()}

json.dump(res, open(os.path.join(HERE, "height_entropy.json"), "w"), indent=1, default=float)
nb = sum(v for v in res.values() if isinstance(v, bool)); nt = sum(isinstance(v, bool) for v in res.values())
print(f"\n{nb} of {nt} checks pass")
