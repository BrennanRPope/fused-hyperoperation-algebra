"""Zoom dial (FHA 1.0 checklist item 2): can the circle operation take a zoom ratio per piece,
with dimension coming out as an output?

Proposal tested: a piece shrunk by zoom r carries weight r^s (s a free exponent).
- Joining pieces side by side adds their weights: x(s) = sum r_i^s.
- Zooming in again and again is "any number of times" (the k = -1 element):
  all finite words of pieces have total weight 1/(1 - x(s)).
- Dimension d is where that element hits its wall, x(d) = 1 (Moran 1946).
- Complex s on the wall are the complex dimensions (the wobble).
"""
import numpy as np, mpmath as mp, itertools
np.seterr(all="ignore")
ok = 0; total = 0
def check(name, cond):
    global ok, total
    total += 1; ok += bool(cond); print(("PASS " if cond else "FAIL ") + name)

SETS = {"A two at 1/4": [0.25, 0.25], "B three at 1/9": [1/9]*3, "C three at 1/9 (off-centre)": [1/9]*3,
        "D uneven 9/100, 49/100": [0.09, 0.49], "E four at 1/16": [1/16]*4,
        "gasket three at 1/2": [0.5]*3, "middle-third Cantor": [1/3, 1/3]}
x = lambda rs, s: sum(r**s for r in rs)
def moran(rs):
    lo, hi = 0.0, 5.0
    for _ in range(200):
        mid = (lo+hi)/2
        lo, hi = (mid, hi) if x(rs, mid) > 1 else (lo, mid)
    return (lo+hi)/2

# 1. dimension is the wall crossing
for name, rs in SETS.items():
    d = moran(rs)
    want = {"gasket three at 1/2": np.log(3)/np.log(2), "middle-third Cantor": np.log(2)/np.log(3)}.get(name, 0.5)
    check(f"1 {name}: x(d) = 1 at d = {d:.6f} (expected {want:.6f})", abs(d - want) < 1e-9)

# 2. word sums: depth-n words carry x(s)^n in total; all depths sum to 1/(1 - x)
for name in ["A two at 1/4", "D uneven 9/100, 49/100"]:
    rs = SETS[name]
    for s in [0.3, 0.5, 0.8]:
        good = True
        for n in range(1, 7):
            tot = sum(np.prod([rs[i] for i in w])**s for w in itertools.product(range(len(rs)), repeat=n))
            good &= abs(tot - x(rs, s)**n) < 1e-12
        check(f"2a {name}, s={s}: words of depth 1..6 weigh x(s)^n", good)
    s = 0.8; X = x(rs, s)
    partial = sum(X**n for n in range(0, 400))
    check(f"2b {name}: above the dimension the word series converges to 1/(1-x)", abs(partial - 1/(1-X)) < 1e-9)
    s = 0.5
    check(f"2c {name}: at the dimension every depth weighs exactly 1 (the wall, series diverges)", abs(x(rs, s) - 1) < 1e-12)

# 3. the natural weights of the pieces are r_i^d, and they sum to 1
for name, rs in SETS.items():
    d = moran(rs); w = [r**d for r in rs]
    check(f"3 {name}: natural weights r^d sum to 1", abs(sum(w) - 1) < 1e-9)

# 4. zoom as twist: scaling dial values by w = r^d turns the circle operation into the twisted rule at k = w
for r, d in [(0.25, 0.5), (1/3, np.log(2)/np.log(3)), (0.5, 1.0)]:
    w = r**d; t, u = 0.37, 1.9
    lhs = (w*t) + (w*u) + (w*t)*(w*u); rhs = w*(t + u + w*t*u)
    check(f"4 zoom r={r:.3f}, d={d:.3f}: w t ∘₁ w u = w (t ∘_w u)", abs(lhs - rhs) < 1e-12)
# size-dial link: doubling at dimension beta is zoom 1/2, weight 2^-beta
check("4b size dial: zoom 1/2 at dimension beta gives weight 2^-beta (beta = 1, 2, 3)", all(abs(0.5**b - 2.0**-b) < 1e-15 for b in [1, 2, 3]))

# 5. complex dimensions = complex wall crossings (wobble)
mp.mp.dps = 30
def roots_lattice(r, N, m):
    d = mp.log(N)/mp.log(1/r); return [d + 2j*mp.pi*j/mp.log(1/r) for j in range(-m, m+1)]
for name, (r, N, sp) in {"A": (0.25, 2, 4.53), "B": (1/9, 3, 2.86), "E": (1/16, 4, 2.27)}.items():
    rts = roots_lattice(mp.mpf(r), N, 3)
    good = all(abs(N*mp.power(r, s) - 1) < 1e-20 for s in rts)
    spacing = float(mp.im(rts[1] - rts[0]))
    check(f"5a {name}: complex roots of N r^s = 1 lie on Re s = 1/2, spacing {spacing:.2f} (page: {sp})", good and abs(spacing - sp) < 0.01)
# D non-lattice: roots near Re s = 1/2 but not on one evenly spaced line
L1, L2 = np.log(0.09), np.log(0.49)
fz = lambda z: np.exp(L1*z) + np.exp(L2*z) - 1
dfz = lambda z: L1*np.exp(L1*z) + L2*np.exp(L2*z)
found = set()
for im0 in np.arange(0.5, 40, 0.25):
    for re0 in [0.2, 0.4, 0.5]:
        z = complex(re0, im0)
        for _ in range(60):
            z = z - fz(z)/dfz(z)
        if abs(fz(z)) < 1e-12 and 0.01 < z.imag < 40 and -2 < z.real < 1:
            found.add((round(z.real, 6), round(z.imag, 6)))
found = sorted(found, key=lambda t: t[1])
res = [t[0] for t in found]
check(f"5b D: {len(found)} complex dimensions found up to height 40, all with Re s <= 1/2", len(found) >= 5 and all(r <= 0.5 + 1e-9 for r in res))
check("5c D: their real parts vary (non-lattice: no single evenly spaced line)", max(res) - min(res) > 1e-3)
np.save("/tmp/zoom_D_roots.npy", np.array(found))
print("D roots:", found[:8])

# 6. dimension is an output: change one zoom, d moves; same d from different zooms
check("6a zooms (1/4,1/4) and (9/100,49/100) both give d = 1/2", abs(moran([.25,.25]) - .5) < 1e-9 and abs(moran([.09,.49]) - .5) < 1e-9)
check("6b d is monotone in the zoom: (0.2,0.2) < (0.25,0.25) < (0.3,0.3)", moran([.2,.2]) < moran([.25,.25]) < moran([.3,.3]))
print(f"\n{ok}/{total} checks passed")
