"""Rigidity check (DRAFT). A genus with exponential f gives, by fixed-point localization on CP^n,
phi_t = sum_i prod_{j!=i} 1/f(t_j - t_i). A rigid genus gives the same number for every choice of
torus weights t. Plane-dial exponentials (the shapes) should be rigid; the rung-2 exponential
E^2(t) = e^(e^t - 1) - 1 (ruler log(1+log(1+x))) should not be."""
import cmath, random
def f_dial(t, k, d):
    disc = cmath.sqrt(k*k - 4*d); a, b = (k+disc)/2, (k-disc)/2
    if abs(a-b) < 1e-12: return t/(1-a*t)
    E = cmath.exp((a-b)*t); return (E-1)/(a-b*E)
f_rung2 = lambda t: cmath.exp(cmath.exp(t)-1)-1
def loc(f, ts):
    n = len(ts); s = 0
    for i in range(n):
        p = 1
        for j in range(n):
            if j != i: p /= f(ts[j]-ts[i])
        s += p
    return s
random.seed(7)
print("CP^2 by localization, five random weight choices:")
for name, f in [("plane dial k=1,d=2 (angle)", lambda t: f_dial(t, 1, 2)),
                ("plane dial k=3,d=1 (radius)", lambda t: f_dial(t, 3, 1)),
                ("rung 2: e^(e^t-1)-1", f_rung2)]:
    vals = []
    for _ in range(5):
        ts = [random.uniform(-0.6, 0.6) for _ in range(3)]
        vals.append(loc(f, ts).real)
    print(f"  {name:30s}: " + "  ".join(f"{v:+.6f}" for v in vals))
