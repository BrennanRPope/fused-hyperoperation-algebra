"""Checks for the Lee-Yang proof pieces (section N3, proof).  Light polish: exact series + zero counts in a disc.
G(s) = int_0^1 (1 - x^p)^a cos(s x) dx,  a = (n - 1)/p for the uniform n-dim p-ball.
"""
import json, time
import numpy as np
import mpmath as mp
RES = []; OUT = {}
def check(name, ok, detail=""):
    RES.append((name, bool(ok))); print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))
def ly_coeffs(p, a, K, dps):
    mp.mp.dps = dps; p = mp.mpf(p); a = mp.mpf(a)
    return [(-1)**k/mp.factorial(2*k)*mp.beta((2*k + 1)/p, a + 1)/p for k in range(K)]
def ly_eval(c, s):
    s2 = mp.mpmathify(s)**2; acc = mp.mpf(0)
    for ck in reversed(c): acc = acc*s2 + ck
    return acc
def zeros_in_disc(p, a, X=60.0, npts=2400):
    c = ly_coeffs(p, a, int(1.5*X) + 80, int(X/2) + 50)
    zs = [X*mp.expj(2*mp.pi*k/npts) for k in range(npts + 1)]
    ang = [float(mp.arg(ly_eval(c, z))) for z in zs]
    tot = float(np.sum(np.diff(np.unwrap(ang))))/(2*np.pi)
    rv = [ly_eval(c, -X + i*0.05) for i in range(int(2*X/0.05) + 1)]
    return round(tot), sum(1 for i in range(len(rv) - 1) if rv[i]*rv[i + 1] < 0)
t0 = time.perf_counter()
# P1  integration by parts: s G(s) = int_0^1 g(x) sin(s x) dx,  g = a p x^(p-1) (1 - x^p)^(a-1)
mp.mp.dps = 30; worst = 0
for p, a, s in [(1.5, 0.7, 3.1), (3.0, 0.4, 7.7), (1.2, 2.0, 5.0)]:
    G = mp.quad(lambda x: (1 - x**p)**a*mp.cos(s*x), [0, 1])
    I = mp.quad(lambda x: a*p*x**(p - 1)*(1 - x**p)**(a - 1)*mp.sin(s*x), [0, 1])
    worst = max(worst, abs(s*G - I))
check("P1 s*G(s) equals the sine transform of g = a p x^(p-1)(1-x^p)^(a-1)", worst < 1e-10, f"max gap {float(worst):.1e}")
# P2  g is nondecreasing on (0,1) exactly when a <= 1 (p >= 1)
x = np.linspace(1e-4, 1 - 1e-6, 20001)
def mono(p, a):
    g = x**(p - 1)*(1 - x**p)**(a - 1); return bool(np.all(np.diff(g) >= -1e-12*np.max(g)))
grid = [(p, a, mono(p, a)) for p in [1.0, 1.3, 2.0, 3.5, 6.0] for a in [0.3, 0.8, 1.0, 1.05, 1.5, 3.0]]
check("P2 g rises all the way to the rim when a <= 1 and turns down near the rim when a > 1",
      all(m == (a <= 1) for p, a, m in grid), f"{len(grid)} (p, a) pairs")
# P3  p = 4, integer a: G_4,a = (1 - D^2)^a F_a, where F_a is the p = 2 (Bessel) transform of (1 - x^2)^a
def series_fn(p, a, K=90): c = ly_coeffs(p, a, K, 40); return lambda s: ly_eval(c, s)
worst = 0
for a in [1, 2, 3]:
    G4 = series_fn(4, a); F = series_fn(2, a)
    for s in [0.7, 2.3, 5.9]:
        lhs = G4(s)
        # (1 - D^2)^a F = sum_j C(a,j) (-1)^j F^(2j)
        rhs = sum(mp.binomial(a, j)*(-1)**j*mp.diff(F, s, 2*j) for j in range(a + 1))
        worst = max(worst, abs(lhs - rhs))
check("P3 at p = 4 with whole a, G equals (1 - D^2)^a applied to the round-ball (Bessel) transform", worst < 1e-15,
      f"max gap {float(worst):.1e}")
# P4  zero counts in |s| < 60: proved cases and gap cases
cases = [(4.0, 2.0, "proved"), (4.0, 3.0, "proved"), (1.5, 0.95, "proved"), (3.0, 1.0, "proved"),
         (4.0, 1.5, "gap"), (4.0, 2.5, "gap"), (6.0, 2.0, "gap"), (6.0, 4.0, "gap"),
         (1.5, 1.2, "gap"), (3.0, 2.0, "gap"), (5.0, 4.5, "gap"), (6.0, 9.0, "gap (a > p, even p)")]
disc = {}
for p, a, kind in cases:
    disc[f"p={p:g}, a={a:g} ({kind})"] = zeros_in_disc(p, a)
print("    zeros in |s| < 60 (all, real):", disc)
check("P4 every zero inside |s| < 60 is real in all 12 cases, proved and gap", all(v[0] == v[1] for v in disc.values()))
OUT["disc"] = disc; OUT["seconds"] = round(time.perf_counter() - t0, 1)
json.dump(OUT, open("proof_checks.json", "w"), indent=1)
print(f"\n{sum(ok for _, ok in RES)} of {len(RES)} checks pass  ({OUT['seconds']} s)")
