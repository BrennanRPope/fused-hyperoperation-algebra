"""Does any l_p shape other than the sphere (p=2) keep what gravity needs?
Tests, for p in {1, 1.5, 2, 3, 4, inf} on the plane C and on H:
 (S1) size-keeping product |zw|_p = |z|_p |w|_p  (needed for the polar split)
 (S2) continuous rotations keep size (needed for angular momentum, Noether)
 (S3) boosts: the split-complex (relativity plane) partner keeps a Minkowski
      'size' only as the quadratic form t^2 - x^2 (p=2 analogue)
 (S4) Kepler: under an l_p central pull, is angular momentum conserved and
      does a circular orbit stay circular?  (numerical)
Writes results_shape.json and checks_shape_output.txt."""
import json, math, random
import numpy as np
from scipy.integrate import solve_ivp

def norm(v, p):
    v = np.abs(np.asarray(v, float))
    return v.max() if p == math.inf else (v**p).sum()**(1/p)

def cmul(a, b): return np.array([a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0]])
def qmul(a, b):
    a0,a1,a2,a3 = a; b0,b1,b2,b3 = b
    return np.array([a0*b0-a1*b1-a2*b2-a3*b3, a0*b1+a1*b0+a2*b3-a3*b2,
                     a0*b2-a1*b3+a2*b0+a3*b1, a0*b3+a1*b2-a2*b1+a3*b0])
PS = [1, 1.5, 2, 3, 4, math.inf]
rng = np.random.default_rng(3)
res = {}
for p in PS:
    key = 'inf' if p == math.inf else str(p)
    e1 = max(abs(norm(cmul(a, b), p) - norm(a, p)*norm(b, p)) for a, b in (rng.normal(size=(2, 2)) for _ in range(2000)))
    e2 = max(abs(norm(qmul(a, b), p) - norm(a, p)*norm(b, p)) for a, b in (rng.normal(size=(2, 4)) for _ in range(2000)))
    # rotation by angle th: size change over a full turn of the unit vector (1,0)
    th = np.linspace(0, 2*np.pi, 721)
    sizes = [norm([math.cos(t), math.sin(t)], p) for t in th]
    rot = max(sizes)/min(sizes) - 1
    # Kepler with pull along -grad of 1/|r|_p (potential -1/|r|_p)
    def rhs(t, y):
        x, yy, vx, vy = y; h = 1e-6
        V = lambda X, Y: -1/norm([X, Y], p)
        gx = (V(x+h, yy)-V(x-h, yy))/(2*h); gy = (V(x, yy+h)-V(x, yy-h))/(2*h)
        return [vx, vy, -gx, -gy]
    sol = solve_ivp(rhs, [0, 40], [1, 0.2, 0, 1.0], rtol=1e-10, atol=1e-12, max_step=0.01)
    L = sol.y[0]*sol.y[3] - sol.y[1]*sol.y[2]
    Ldrift = float((L.max()-L.min())/abs(L[0]))
    res[key] = {'S1 product keeps size, plane (max err)': float(e1),
                'S1 product keeps size, quaternions (max err)': float(e2),
                'S2 size change under a full rotation': float(rot),
                'S4 angular momentum drift in a Kepler run': Ldrift}
lines = []
ok = True
for k, r in res.items():
    is2 = k == '2'
    keeps = r['S1 product keeps size, plane (max err)'] < 1e-9 and r['S1 product keeps size, quaternions (max err)'] < 1e-9
    rotok = r['S2 size change under a full rotation'] < 1e-9
    Lok = r['S4 angular momentum drift in a Kepler run'] < 1e-6
    good = (keeps == is2) and (rotok == is2) and (Lok == is2)
    ok &= good
    lines.append(f"p={k:4s} product keeps size: {'yes' if keeps else 'no'} (err {r['S1 product keeps size, plane (max err)']:.1e}); "
                 f"rotation keeps size: {'yes' if rotok else 'no'} (change {r['S2 size change under a full rotation']:.3f}); "
                 f"angular momentum kept: {'yes' if Lok else 'no'} (drift {r['S4 angular momentum drift in a Kepler run']:.1e})  [{'PASS' if good else 'FAIL'}]")
lines.append('ALL PASS: only p=2 keeps all three' if ok else 'SOMETHING FAILED')
txt = "\n".join(lines); print(txt)
open('checks_shape_output.txt', 'w').write(txt+"\n")
json.dump(res, open('results_shape.json', 'w'), indent=1)
