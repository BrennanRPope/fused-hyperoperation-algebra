"""The s-ladder at real rung height, computed numerically.

E(x) = e^x - 1 fixes 0 with slope 1 (a parabolic fixed point), so its real iterates E^s
have one standard smooth choice: regular iteration (Szekeres 1958, Ecalle). We compute
it the classical way: walk toward 0 with n steps, apply the formal series of E^s there
(where it is extremely accurate), walk back out.

  x > 0: E^s(x) = E^n( P_s( L^n(x) ) )      L(x) = log(1+x) pulls positives to 0
  x < 0: E^s(x) = L^n( P_s( E^n(x) ) )      E pulls negatives to 0

B_s(a,b) = E^s( E^-s(a) + E^-s(b) ):  s = 0 addition, s = 1 circle operation, s = -1 log(e^a+e^b-1).
"""
from fractions import Fraction as Fr
import math, json, sys
import numpy as np

ORDER = 14
NSTEP = 600

# ---- exact formal coefficients of E^s as polynomials in s (Lagrange from integer iterates) ----
def compose(f, g, order):
    """coefficient lists (index = power, f[0]=g[0]=0) -> f(g(x)) truncated"""
    out = [Fr(0)]*(order + 1)
    gp = [Fr(0)]*(order + 1); gp[0] = Fr(1)  # g^0
    for i in range(1, order + 1):
        # gp <- gp * g
        new = [Fr(0)]*(order + 1)
        for p, cp in enumerate(gp):
            if cp == 0: continue
            for q in range(1, order + 1 - p):
                if g[q]: new[p + q] += cp*g[q]
        gp = new
        if f[i]:
            for j in range(order + 1):
                out[j] += f[i]*gp[j]
    return out

def formal_iterates(order):
    E = [Fr(0)] + [Fr(1, math.factorial(i)) for i in range(1, order + 1)]
    its = [[Fr(0), Fr(1)] + [Fr(0)]*(order - 1)]
    for n in range(1, order + 1):
        its.append(compose(E, its[-1], order))
    return its  # its[n][i] = coeff of x^i in E^n, degree <= i-1 in n

ITS = formal_iterates(ORDER)

def coeffs_at(s, order=ORDER):
    """coefficient i of E^s = Lagrange interpolation through n = 0..i-1 (exact polynomial of degree i-1)"""
    cs = [0.0]*(order + 1)
    for i in range(1, order + 1):
        pts = list(range(0, i))
        val = 0.0
        for j in pts:
            w = 1.0
            for m in pts:
                if m != j: w *= (s - m)/(j - m)
            val += float(ITS[j][i])*w
        cs[i] = val
    return cs

def P(s, y):
    cs = coeffs_at(s)
    out = np.zeros_like(y)
    for i in range(ORDER, 0, -1):
        out = (out + cs[i])*y
    return out

def Es(s, x, n=NSTEP):
    x = np.asarray(x, dtype=float)
    out = np.full_like(x, np.nan)
    pos = x > 0; neg = (x < 0) & (x > -np.inf); zero = x == 0
    out[zero] = 0.0
    if pos.any():
        y = x[pos].copy()
        for _ in range(n): y = np.log1p(y)
        y = P(s, y)
        with np.errstate(over='ignore', invalid='ignore'):
            for _ in range(n): y = np.expm1(y)
        out[pos] = y
    if neg.any():
        y = x[neg].copy()
        for _ in range(n): y = np.expm1(y)
        y = P(s, y)
        with np.errstate(invalid='ignore'):
            for _ in range(n):
                y = np.where(y > -1, np.log1p(np.maximum(y, -1 + 1e-300)), np.nan)
        out[neg] = y
    return out

def B(s, a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    return Es(s, Es(-s, a) + Es(-s, b))

if __name__ == "__main__":
    res = {}
    def rep(name, ok, extra=""):
        print(("PASS " if ok else "FAIL ") + name + (("   " + extra) if extra else ""))
        res[name] = bool(ok)

    xs = np.array([-0.9, -0.5, -0.1, 0.05, 0.3, 1.0, 2.0, 3.0])
    e1 = np.max(np.abs(Es(1, xs) - np.expm1(xs)))
    rep("E^1 = e^x - 1 on both sides of 0", e1 < 1e-9, f"max err {e1:.1e}")
    e2 = np.max(np.abs(Es(-1, xs) - np.log1p(xs)))
    rep("E^-1 = log(1+x)", e2 < 1e-9, f"max err {e2:.1e}")
    h = Es(0.5, Es(0.5, xs))
    e3 = np.max(np.abs(h - np.expm1(xs)))
    rep("half rung twice = one rung: E^1/2(E^1/2(x)) = e^x - 1", e3 < 1e-8, f"max err {e3:.1e}")
    e4 = np.max(np.abs(Es(0.3, Es(0.45, xs)) - Es(0.75, xs)))
    rep("E^0.3 o E^0.45 = E^0.75", e4 < 1e-8, f"max err {e4:.1e}")
    e5 = np.max(np.abs(Es(1/3, Es(1/3, Es(1/3, xs))) - np.expm1(xs)))
    rep("third rung three times = one rung", e5 < 1e-8, f"max err {e5:.1e}")
    # independence from the step count (the construction converges)
    e6 = np.max(np.abs(Es(0.5, xs, n=400) - Es(0.5, xs, n=1600)))
    rep("result independent of step count (400 vs 1600)", e6 < 1e-8, f"diff {e6:.1e}")

    # s-ladder values at whole rungs
    rep("B_0(1.3, 0.7) = 2.0", abs(B(0, 1.3, .7) - 2.0) < 1e-9)
    rep("B_1(1.3, 0.7) = the circle operation = 2.91", abs(B(1, 1.3, .7) - (1.3 + .7 + .91)) < 1e-9)
    rep("B_-1(1.3, 0.7) = log(e^1.3 + e^0.7 - 1)", abs(B(-1, 1.3, .7) - math.log(math.exp(1.3) + math.exp(.7) - 1)) < 1e-9)
    b2 = (1+1.3)*(1+.7)*math.exp(math.log1p(1.3)*math.log1p(.7)) - 1
    rep("B_2(1.3, 0.7) = (1+a)(1+b)e^{L(a)L(b)} - 1", abs(B(2, 1.3, .7) - b2) < 1e-8)

    # group laws at fractional heights
    worst_assoc = 0; worst_comm = 0; worst_id = 0
    rng = np.random.default_rng(1)
    for s in [-0.5, 0.25, 0.5, 0.75, 1.5]:
        p = rng.uniform(0.05, 1.5, 6); q = rng.uniform(0.05, 1.5, 6); r = rng.uniform(0.05, 1.5, 6)
        worst_assoc = max(worst_assoc, np.nanmax(np.abs(B(s, B(s, p, q), r) - B(s, p, B(s, q, r)))))
        worst_comm = max(worst_comm, np.nanmax(np.abs(B(s, p, q) - B(s, q, p))))
        worst_id = max(worst_id, np.nanmax(np.abs(B(s, p, 0*p) - p)))
    rep("B_s associative at s = -0.5, 0.25, 0.5, 0.75, 1.5", worst_assoc < 1e-8, f"worst {worst_assoc:.1e}")
    rep("B_s commutative", worst_comm < 1e-9, f"worst {worst_comm:.1e}")
    rep("0 is the identity at every height", worst_id < 1e-9, f"worst {worst_id:.1e}")

    # corner weight = s
    hh = 1e-3
    cw = {s: float((B(s, hh, hh) - 2*hh)/hh**2) for s in [-1, -0.5, 0, 0.25, 0.5, 1, 1.5, 2]}
    okc = all(abs(cw[s] - s) < 5e-3*max(1, abs(s)) + 2e-3 for s in cw)
    rep("corner weight (B_s(h,h) - 2h)/h^2 -> s", okc, json.dumps({str(k): round(v, 4) for k, v in cw.items()}))

    # the smooth curve through all rungs: B_s(1,1)
    ss = [-1, -0.5, 0, 0.5, 1, 1.5, 2]
    curve = {str(s): float(B(s, 1.0, 1.0)) for s in ss}
    print("B_s(1,1):", json.dumps({k: round(v, 6) for k, v in curve.items()}))
    res["B_s(1,1)"] = curve
    # half-iterate value
    print("E^(1/2)(1) =", float(Es(0.5, np.array([1.0]))[0]), " E^(1/2)(-0.5) =", float(Es(0.5, np.array([-0.5]))[0]))

    # smooth across 0: numeric derivatives from both sides match the formal coefficients
    cs = coeffs_at(0.5)
    for side in (1, -1):
        hs = side*np.array([1e-2, 2e-2, 3e-2])
        vals = Es(0.5, hs)
        # fit quadratic coefficient from values
        q2 = (vals - hs)/hs**2
        print(f"side {side:+d}: (E^1/2(x)-x)/x^2 at |x|=.01,.02,.03 ->", np.round(q2, 5), " formal c2 =", cs[2], " c3 =", round(cs[3], 6))
    # Baker: formal half-iterate diverges (coefficients grow like factorials)
    big = formal_iterates(26)
    def coeff_half(i):
        pts = list(range(0, i)); val = Fr(0)
        for j in pts:
            w = Fr(1)
            for m in pts:
                if m != j: w *= (Fr(1, 2) - m)/Fr(j - m)
            val += big[j][i]*w
        return val
    hc = [coeff_half(i) for i in range(1, 27)]
    mags = [abs(float(c)) for c in hc]
    win = [max(mags[i:i+5]) for i in range(5, 26, 5)]  # window maxima, n = 6-10, 11-15, ...
    roots = [w**(1/(6+5*j)) for j, w in enumerate(win)]
    print("half-rung formal coefficients: window max |c_n| for n=6-10,11-15,16-20,21-25:", ["%.2e" % w for w in win[:4]])
    print("  rough n-th roots:", [round(r, 3) for r in roots[:4]])
    rep("formal half-rung series diverges (n-th roots of coefficients keep growing; Baker 1958)",
        all(roots[i+1] > roots[i] for i in range(3)), " -> ".join("%.3f" % r for r in roots[:4]))

    # tropical / zoom limits of the negative rung
    for H in [1.0, 0.1, 0.01]:
        z = H*math.log(math.exp(2/H) + math.exp(1/H) - 1)
        print(f"zoom out B_-1: h={H}: h B_-1(2/h, 1/h) = {z:.6f}  (max = 2)")
    for H in [10, 100, 1000]:
        z = H*math.log(math.exp(2/H) + math.exp(1/H) - 1)
        print(f"zoom in  B_-1: h={H}: h B_-1(2/h, 1/h) = {z:.6f}  (sum = 3)")
    # lower rungs on positives approach max
    low = {str(s): float(B(s, 2.0, 1.0)) for s in [-1, -2, -3]}
    print("B_s(2,1) for s = -1..-3:", json.dumps({k: round(v, 6) for k, v in low.items()}))
    rep("lower rungs approach max on positives: B_-3(2,1) within 1e-6 of 2 (lower still overflows floats)", abs(low["-3"] - 2) < 1e-6)
    lowneg = {str(s): float(B(s, -0.3, -0.2)) for s in [-0.5, -1, -2]}
    print("B_s(-0.3,-0.2):", lowneg, " (nan = past the wall)")

    json.dump(res, open(sys.argv[1] if len(sys.argv) > 1 else "/dev/null", "w"), indent=1)
    print(f"\n{sum(v for v in res.values() if isinstance(v, bool))} of {sum(isinstance(v, bool) for v in res.values())} numeric checks pass")
