"""e and the natural log at every rung height.

E^s = regular iteration of E(x) = e^x - 1 at its parabolic fixed point 0 (Szekeres, Ecalle);
computed by ../regular_iteration.py (16/16 checks).  B_s(a,b) = E^s(E^-s(a) + E^-s(b)).

Claim (ours): building e at height s the way e is built from the n-cube,
    e_s := 1 + lim_n  (1/n combined n times with B_s)   [the 1 is the cube's empty corner]
equals 1 + E^s(1).
Reason: B_s is addition seen through E^s, so the n-fold B_s of x is E^s(n E^-s(x)) exactly,
and n E^-s(1/n) = 1 - s/(2n) + O(1/n^2) -> 1 because E^-s(x) = x - s x^2/2 + ...
So the error should be about -(s/(2n)) (E^s)'(1).
"""
import sys, os, math, json
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
import numpy as np
from regular_iteration import Es, B, ITS, ORDER

out = {}
def rep(name, ok, extra=""):
    print(("PASS " if ok else "FAIL ") + name + (("   " + extra) if extra else ""))
    out.setdefault("checks", {})[name] = bool(ok)

def e_s(s):
    return 1.0 + float(Es(s, np.array([1.0]))[0])

def dEs_dx(s, x, h=1e-5):
    return float((Es(s, np.array([x + h])) - Es(s, np.array([x - h])))[0] / (2*h))

# ---------- 1. literal n-fold build ----------
S = [-2.0, -1.0, -0.5, 0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0]
NS = [10, 30, 100, 300, 1000]
build = {}
print("n-fold build of 1/n with B_s, literally combined one at a time")
for s in S:
    target = e_s(s)
    row = {}
    for n in NS:
        y = np.array([1.0/n]); u = np.array([1.0/n])
        for _ in range(n - 1):
            y = B(s, y, u)
        row[n] = 1.0 + float(y[0])
    pred = -(s/2) * dEs_dx(s, 1.0)  # predicted n*(error)
    build[s] = dict(target=target, values=row, pred_n_err=pred)
    errs = {n: row[n] - target for n in NS}
    print(f" s={s:+.2f}: 1+E^s(1)={target:.10f}  n=1000 build={row[1000]:.10f}  "
          f"n*err at n=100,1000: {100*errs[100]:+.5f}, {1000*errs[1000]:+.5f}  predicted {pred:+.5f}")
ok1 = all(abs(build[s]["values"][1000] - build[s]["target"]) < 2e-3*max(1, abs(build[s]["pred_n_err"])) for s in S)
rep("n-fold build converges to 1 + E^s(1) at all 10 heights (n=1000)", ok1)
ok1b = all(abs(1000*(build[s]["values"][1000]-build[s]["target"]) - build[s]["pred_n_err"]) < 0.01*max(1, abs(build[s]["pred_n_err"])) for s in S)
rep("error matches the predicted -(s/2n)(E^s)'(1) to 1%", ok1b)

# ---------- known heights ----------
rep("e_0 = 2", abs(e_s(0) - 2) < 1e-10)
rep("e_1 = e", abs(e_s(1) - math.e) < 1e-9)
rep("e_2 = e^(e-1)", abs(e_s(2) - math.exp(math.e - 1)) < 1e-8)
rep("e_-1 = 1 + log 2", abs(e_s(-1) - (1 + math.log(2))) < 1e-10)
rep("e_-2 = 1 + log(1 + log 2)", abs(e_s(-2) - (1 + math.log1p(math.log(2)))) < 1e-10)
rep("e_3 = e^(e^(e-1)-1)", abs(e_s(3) - math.exp(math.exp(math.e - 1) - 1)) < 1e-6)
table = {s: e_s(s) for s in [-3, -2, -1, -0.5, 0, 1/3, 0.5, 2/3, 1, 1.5, 2, 2.5, 3]}
print("e_s table:", json.dumps({f"{k:.4g}": round(v, 10) for k, v in table.items()}))
out["table"] = {f"{k:.6g}": v for k, v in table.items()}
# half rung twice
h = e_s(0.5) - 1
rep("E^1/2 applied to (e_1/2 - 1) gives e - 1", abs(float(Es(0.5, np.array([h]))[0]) - (math.e - 1)) < 1e-8)

# ---------- 2. closed forms ----------
# (a) Levy's one-line limit:  E^s(x) = lim E^n( y/(1 - s y/2) ),  y = L^n(x)   (x > 0)
def levy(s, x, n):
    y = x
    for _ in range(n): y = math.log1p(y)
    y = y/(1 - s*y/2)
    for _ in range(n): y = math.expm1(y)
    return y
lev = {}
for n in [10, 100, 1000, 10000, 100000]:
    lev[n] = 1 + levy(0.5, 1.0, n)
print("Levy one-line formula for e_1/2:", {n: round(v, 10) for n, v in lev.items()}, " target", round(e_s(.5), 10))
out["levy_e_half"] = lev
rep("Levy's formula converges to e_1/2 (n=1e5 within 1e-4)", abs(lev[100000] - e_s(.5)) < 1e-4)
# (a') with the log term of the Abel function: alpha(y) = -2/y + (1/3) log y  (ours to derive, standard method)
def alpha0(y): return -2/y + math.log(y)/3
def alpha0_inv(w, guess):
    y = guess
    for _ in range(60):
        f = alpha0(y) - w; fp = 2/y**2 + 1/(3*y)
        y -= f/fp
    return y
def levy2(s, x, n):
    y = x
    for _ in range(n): y = math.log1p(y)
    y = alpha0_inv(alpha0(y) + s, y)
    for _ in range(n): y = math.expm1(y)
    return y
lev2 = {n: 1 + levy2(0.5, 1.0, n) for n in [10, 100, 1000, 10000]}
print("Levy with log correction:", {n: round(v, 12) for n, v in lev2.items()})
out["levy2_e_half"] = lev2
rep("Levy + log term converges faster (n=1e4 within 1e-7)", abs(lev2[10000] - e_s(.5)) < 1e-7)

# (b) e_s is an entire function of s (inverse Fatou coordinate of an entire map extends to C)
#     -> Taylor series in s.  Compute coefficients by Cauchy integral with complex s.
def coeffs_c(s):
    cs = [0j]*(ORDER + 1)
    for i in range(1, ORDER + 1):
        pts = range(0, i); val = 0j
        for j in pts:
            w = 1+0j
            for m in pts:
                if m != j: w *= (s - m)/(j - m)
            val += float(ITS[j][i])*w
        cs[i] = val
    return cs
def Es_c(s, x, n=600):  # x > 0 only
    y = x
    for _ in range(n): y = math.log1p(y)
    cs = coeffs_c(s); p = 0j
    for i in range(ORDER, 0, -1): p = (p + cs[i])*y
    y = p
    for _ in range(n): y = np.expm1(y)
    return y
M = 64; R = 1.0
pts = [R*np.exp(2j*np.pi*k/M) for k in range(M)]
vals = [1 + Es_c(z, 1.0) for z in pts]
taylor = [np.real(sum(v*np.exp(-2j*np.pi*k*m/M) for k, v in enumerate(vals))/M)/R**m for m in range(20)]
print("Taylor coefficients of e_s in s (about s=0):")
for m, c in enumerate(taylor[:14]): print(f"  s^{m}: {c:+.12e}")
out["taylor"] = taylor
def taylor_eval(s, K=20): return sum(taylor[m]*s**m for m in range(K))
tt = {s: (taylor_eval(s), e_s(s)) for s in [-1, -0.5, 0.5, 1]}
print("series vs direct:", {s: (round(a, 10), round(b, 10)) for s, (a, b) in tt.items()})
rep("power series in s reproduces e_s on [-1, 1] to 1e-8", all(abs(a - b) < 1e-8 for a, b in tt.values()))
rep("first Taylor coefficient matches the numeric d/ds of e_s (Julia iterative logarithm at 1)",
    abs(taylor[1] - (e_s(1e-4) - e_s(-1e-4))/2e-4) < 1e-6, f"c1 = {taylor[1]:.10f}")
# Julia/iterative logarithm lambda(x) = d/ds E^s(x) at s=0 satisfies lambda(E(x)) = E'(x) lambda(x)
lam = lambda x: (float(Es(1e-4, np.array([x]))[0]) - float(Es(-1e-4, np.array([x]))[0]))/2e-4
rep("Julia equation lambda(e^x-1) = e^x lambda(x) at x=1", abs(lam(math.e - 1) - math.e*lam(1.0)) < 1e-5,
    f"{lam(math.e-1):.8f} vs {math.e*lam(1.0):.8f}")
# far below: e_s -> 1 like 1 + 2/|s|
far = {s: (e_s(s) - 1)*abs(s) for s in [-10, -30, -100, -300]}
print("|s|(e_s - 1) for s = -10..-300:", {k: round(v, 5) for k, v in far.items()}, " (-> 2, log correction)")
out["far_below"] = far

# ---------- 3. natural log at each height ----------
def ln_s(s, X): return Es(-s, np.asarray(X, float) - 1)
def exp_s(s, t): return 1 + Es(s, np.asarray(t, float))
def otimes(s, X, Y): return 1 + B(s, np.asarray(X, float) - 1, np.asarray(Y, float) - 1)
worst = dict(inv=0, one=0, base=0, slope=0, curv=0, hom=0, ln1=0, ln2=0)
rng = np.random.default_rng(7)
for s in [-1.0, -0.5, 0.25, 0.5, 1.0, 1.5, 2.0]:
    t = rng.uniform(-0.5, 1.2, 6)
    worst["inv"] = max(worst["inv"], np.max(np.abs(ln_s(s, exp_s(s, t)) - t)))
    worst["one"] = max(worst["one"], abs(float(ln_s(s, [1.0])[0])))
    worst["base"] = max(worst["base"], abs(float(ln_s(s, [e_s(s)])[0]) - 1))
    hh = 1e-4
    d1 = float((ln_s(s, [1 + hh]) - ln_s(s, [1 - hh]))[0])/(2*hh)
    d2 = float((ln_s(s, [1 + 1e-3]) - 2*ln_s(s, [1.0]) + ln_s(s, [1 - 1e-3]))[0])/1e-6
    worst["slope"] = max(worst["slope"], abs(d1 - 1))
    worst["curv"] = max(worst["curv"], abs(d2 + s))
    X = rng.uniform(1.05, 2.5, 6); Y = rng.uniform(1.05, 2.5, 6)
    worst["hom"] = max(worst["hom"], np.max(np.abs(ln_s(s, otimes(s, X, Y)) - ln_s(s, X) - ln_s(s, Y))))
X = np.array([0.5, 1.7, 3.0, 9.0])
worst["ln1"] = np.max(np.abs(ln_s(1, X) - np.log(X)))
worst["ln2"] = np.max(np.abs(ln_s(2, X) - np.log1p(np.log(X))))
worst["prod1"] = float(np.max(np.abs(otimes(1, X, X[::-1]) - X*X[::-1])))
for k, v in worst.items(): print(f"  ln_s worst {k}: {float(v):.2e}")
rep("ln_s(exp_s(t)) = t", worst["inv"] < 1e-8)
rep("ln_s(1) = 0 at every height", worst["one"] < 1e-10)
rep("ln_s(e_s) = 1 at every height", worst["base"] < 1e-8)
rep("slope of ln_s at 1 is 1 at every height", worst["slope"] < 1e-6)
rep("curvature of ln_s at 1 is -s (the height itself)", worst["curv"] < 1e-3)
rep("ln_s turns X (x)_s Y into addition", worst["hom"] < 1e-8)
rep("ln_1 = natural log, (x)_1 = ordinary product", worst["ln1"] < 1e-9 and worst["prod1"] < 1e-8)
rep("ln_2(X) = log(1 + log X)", worst["ln2"] < 1e-9)
# Tsallis comparison: ln_q(1+u) = u - q u^2/2 + q(q+1)/6 u^3 ;  ln_s(1+u) = u - s u^2/2 + (s^2/4 + s/12) u^3
u = 1e-2
for s in [0.5, 2.0]:
    ours = float(ln_s(s, [1 + u])[0]); ts = ((1+u)**(1-s) - 1)/(1-s)
    c3o = (ours - (u - s*u*u/2))/u**3; c3t = (ts - (u - s*u*u/2))/u**3
    print(f"  s=q={s}: cubic coeff ours {c3o:.4f} (s^2/4+s/12={s*s/4+s/12:.4f}), Tsallis {c3t:.4f} (q(q+1)/6={s*(s+1)/6:.4f})")
out["worst_ln"] = {k: float(v) for k, v in worst.items()}
out["build"] = {str(s): dict(target=v["target"], values={str(n): x for n, x in v["values"].items()}, pred_n_err=v["pred_n_err"]) for s, v in build.items()}
json.dump(out, open(os.path.join(HERE, "e_every_height.json"), "w"), indent=1)
c = out["checks"]; print(f"\n{sum(c.values())} of {len(c)} checks pass")

# ---------- 4. two structural identities (ours, exact by construction; checked numerically) ----------
# (i) heights compose:  ln_{s+t}(X) = ln_s(1 + ln_t(X))   and  ln_s(e_{s+t}) = e_t - 1
w_g = 0; w_e = 0
for s, t in [(0.5, 0.5), (0.3, 0.45), (-0.5, 1.25), (1.0, -0.4), (0.25, 1.5)]:
    X = np.array([1.3, 2.0, 3.5])
    w_g = max(w_g, float(np.max(np.abs(ln_s(s + t, X) - ln_s(s, 1 + ln_s(t, X))))))
    w_e = max(w_e, abs(float(ln_s(s, [e_s(s + t)])[0]) - (e_s(t) - 1)))
rep("heights compose: ln_(s+t) = ln_s(1 + ln_t)", w_g < 1e-8, f"worst {w_g:.1e}")
rep("ln_s(e_(s+t)) = e_t - 1", w_e < 1e-8, f"worst {w_e:.1e}")
# (ii) ln_s turns ORDINARY multiplication into the rung (1 - s) operation:  ln_s(XY) = B_(1-s)(ln_s X, ln_s Y)
w_m = 0
for s in [-0.5, 0.0, 0.25, 0.5, 0.75, 1.0, 1.5]:
    X = np.array([1.2, 1.8, 2.5]); Y = np.array([1.5, 1.1, 2.2])
    w_m = max(w_m, float(np.max(np.abs(ln_s(s, X*Y) - B(1 - s, ln_s(s, X), ln_s(s, Y))))))
rep("ln_s(XY) = B_(1-s)(ln_s X, ln_s Y)  (s=1: ln turns x into +;  s=0: X-1 turns x into a+b+ab)", w_m < 1e-8, f"worst {w_m:.1e}")

# ---------- 5. the height flow as an ODE: Julia's iterative logarithm (beta function reading) ----------
from fractions import Fraction as Fr
def dcoef(i):  # d/ds at s=0 of the Lagrange polynomial through n=0..i-1 for coefficient i
    pts = list(range(i)); tot = Fr(0)
    for j in pts:
        # derivative at 0 of prod_{m!=j} (s-m)/(j-m)
        others = [m for m in pts if m != j]; den = Fr(1)
        for m in others: den *= (j - m)
        d = Fr(0)
        for a in others:
            p = Fr(1)
            for m in others:
                if m != a: p *= (0 - m)
            d += p
        tot += ITS[j][i]*d/den
    return tot
lam_c = [dcoef(i) for i in range(1, 9)]
print("Julia iterative logarithm lambda(x) = dE^s/ds at s=0, series:", ", ".join(f"x^{i+1}: {c}" for i, c in enumerate(lam_c)))
out["julia_series"] = [str(c) for c in lam_c]
rep("lambda(x) = x^2/2 - x^3/12 + x^4/48 - ...", lam_c[1] == Fr(1, 2) and lam_c[2] == Fr(-1, 12) and lam_c[3] == Fr(1, 48))
json.dump(out, open(os.path.join(HERE, "e_every_height.json"), "w"), indent=1)
c = out["checks"]; print(f"\n{sum(c.values())} of {len(c)} checks pass (final)")
