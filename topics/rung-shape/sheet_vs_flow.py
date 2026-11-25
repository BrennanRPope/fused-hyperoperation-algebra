"""Is the rung sheet exactly the "equals 1" slice of the flow surface z = E^s(x)?

Flow surface: points (s, x, E^s(x)).  Rung sheet: points (a, b, s) with B_s(a,b) = 1,
B_s(a,b) = E^s(E^-s a + E^-s b).  Run: python3 sheet_vs_flow.py
"""
import sys, os, math, warnings, json
import numpy as np
from scipy.optimize import brentq
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "fractional-rungs"))
from regular_iteration import Es, B
res = {}
def rep(name, ok, extra=""):
    print(("PASS " if ok else "FAIL ") + name + (("   " + extra) if extra else "")); res[name] = bool(ok)
f = lambda s, x: float(Es(s, np.array([float(x)]))[0])

# 1. the flow surface's own "= 1" set is one curve: x = E^-s(1), the path of 1 run backward
S = np.linspace(-3, 4, 29)
w = np.array([f(-s, 1.0) for s in S])
err1 = max(abs(f(s, wi) - 1) for s, wi in zip(S, w))
mono = all(f(s, wi*0.999) < 1 < f(s, wi*1.001) for s, wi in zip(S, w))
rep("1 flow surface level {E^s(x) = 1} is the single curve x = E^-s(1) (one x per height)", err1 < 1e-9 and mono, f"err {err1:.1e}")
print("   that curve w(s) = E^-s(1):", {float(round(s, 2)): round(float(t), 4) for s, t in list(zip(S, w))[::4]})

# 2. the sheet is the flow applied to straight lines: on the sheet, E^-s(a) + E^-s(b) = E^-s(1)
rng = np.random.default_rng(7)
heights = np.concatenate([[-2, -1, 0, 1, 2, 3], rng.uniform(-3, 4, 20)])
worst = 0
for s in heights:
    a = rng.uniform(0.02, 0.98, 25)
    b = Es(s, f(-s, 1.0) - Es(-s, a))
    worst = max(worst, np.max(np.abs(Es(-s, a) + Es(-s, b) - f(-s, 1.0))))
    worst = max(worst, np.max(np.abs(B(s, a, b) - 1)))
rep("2 at every height tested (6 whole, 20 random), each sheet point maps under E^-s onto the straight line u + v = E^-s(1) (rounding grows with E^-s(1), 96 at s = -3)", worst < 1e-8, f"worst {worst:.1e}")

# 3. where the flow's own path shows up on the sheet: the diagonal point is E^s(w/2)
dg = max(abs(B(s, f(s, f(-s, 1.0)/2), f(s, f(-s, 1.0)/2)) - 1) for s in heights)
rep("3 the sheet's diagonal point at height s is E^s(½ E^-s(1)): the flow carrying the midpoint", dg < 1e-9, f"err {dg:.1e}")

# 4. the rung operation at all levels gives back the flow (normalised by slope 1 at 0)
def recover_lens(s, a, k=14):
    """halve a k times using only the rung-s operation, then scale back up: 2^k * (a/2^k in rung s)"""
    h = a
    for _ in range(k):
        h = brentq(lambda t: float(B(s, np.array([t]), np.array([t]))[0]) - h, 0, h)
    return h*2**k
rec = {}
for s in [0.5, 1.0, -0.7]:
    for a in [0.4, 0.9]:
        rec[(s, a)] = (recover_lens(s, a), f(-s, a))
gap = max(abs(r - t) for r, t in rec.values())
rep("4 halving with the rung operation alone, 14 times, recovers the flow's lens E^-s(a)", gap < 1e-3,
    "; ".join(f"s={s},a={a}: {r:.5f} vs {t:.5f}" for (s, a), (r, t) in rec.items()))

# 5. the c = 1 slice alone does not pin the rung at a fixed height: same unit curve, different law
s = 0.5; W = f(-s, 1.0); EPS = 0.02
h = lambda u: u + EPS*np.sin(2*np.pi*u/W)**3   # h(W - u) = W - h(u), h(0)=0, h'(0)=1
def hinv(y):
    return brentq(lambda u: h(u) - y, -1, 5)
psi = lambda x: h(f(-s, x)); phi = lambda y: f(s, hinv(y))
Bt = lambda a, b: phi(psi(a) + psi(b))
aa = np.linspace(0.05, 0.95, 19)
unit_gap = max(abs(Bt(a, float(Es(s, W - f(-s, a)))) - 1) for a in aa)
law_gap = max(abs(Bt(a, b) - float(B(s, a, b))) for a in aa[::3] for b in aa[::3])
assoc = max(abs(Bt(Bt(p, q), r) - Bt(p, Bt(q, r))) for p, q, r in [(0.2, 0.3, 0.4), (0.5, 0.1, 0.25), (0.7, 0.6, 0.05)])
rep("5 a different associative rung-½ law has exactly the same unit curve", unit_gap < 1e-9 and law_gap > 1e-3 and assoc < 1e-9,
    f"unit-curve gap {unit_gap:.1e}, law differs by up to {law_gap:.4f}, associativity err {assoc:.1e}")
lvl_half = max(abs(Bt(a, float(Es(s, f(-s, 0.5) - f(-s, a)))) - 0.5) for a in np.linspace(0.05, 0.45, 9))
rep("5b ...but its level-½ curve differs: the other level sets carry the extra information", lvl_half > 1e-4, f"off by up to {lvl_half:.4f}")

json.dump({k: v for k, v in res.items()}, open("sheet_vs_flow_results.json", "w"), indent=1)
print(f"\n{sum(res.values())} of {len(res)} checks pass")
