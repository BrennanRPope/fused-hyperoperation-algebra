"""Round 5, pair 2 add-on: search for complex dimensions of the prime fractal (roots of P(s) = 1, Re s > 1),
and the tail of the direct Dirichlet sum at s = 2 behind the identity check."""
import json, math, numpy as np, mpmath as mp
mp.mp.dps = 15
f = lambda z: mp.primezeta(z) - 1
roots = []
seeds = [(a, b) for b in np.arange(1.0, 200.0, 1.0) for a in (1.08, 1.22, 1.36)]
for a, b in seeds:
    try:
        r = mp.findroot(f, mp.mpc(a, b), tol=1e-18, maxsteps=40)
    except Exception:
        continue
    if 1.0 < r.real < 1.45 and r.imag > 0.5 and abs(f(r)) < 1e-9:
        rr = (round(float(r.real), 7), round(float(r.imag), 7))
        if all(abs(rr[0]-x) + abs(rr[1]-y) > 1e-5 for x, y in roots): roots.append(rr)
roots.sort(key=lambda p: p[1])
# largest |P(sigma+it)| level check: on 0 < t < 200 how close does P get to 1 at sigma = 1.05 and 1.2?
near = {}
for sg in (1.02, 1.05, 1.2):
    ts = np.arange(0.5, 200, 0.05)
    v = np.array([abs(complex(mp.primezeta(mp.mpc(sg, t))) - 1) for t in ts[::10]])
    near[str(sg)] = dict(t_min=float(ts[::10][v.argmin()]), min_dist=float(v.min()))
out = dict(complex_roots_t_lt_200=roots, closest_approach=near)
print(out)
json.dump(out, open('pair2_complexdims.json', 'w'), indent=1)
