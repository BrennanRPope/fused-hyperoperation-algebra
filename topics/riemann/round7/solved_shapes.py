"""Round 7, item 5: the hypothesis on shapes where it is a theorem, read through the dials.

Curvature dial local factor 1 + k u + d u^2, Delta = k^2 - 4d.  Roots on |u| = d^(-1/2) iff Delta <= 0.
 (a) Ramanujan graphs (Ihara): factor 1 - lambda u + q u^2 per adjacency eigenvalue; RH iff |lambda| <= 2 sqrt q.
 (b) Elliptic curve over F_p (Hasse): factor 1 - a_p u + p u^2; RH iff a_p^2 <= 4p (always true: Hasse 1933).
 (c) Lee-Yang (odds/fugacity dial): ferromagnetic Ising zeros in the fugacity lie on |z| = 1.
 Riemann's own local factor 1 - u has d = 0: the local test is empty; RH is global.
"""
import json, itertools, math
import numpy as np

OUT = {}
def petersen():
    E = [(i, (i+1) % 5) for i in range(5)] + [(5+i, 5+(i+2) % 5) for i in range(5)] + [(i, i+5) for i in range(5)]
    return 10, E
def bottleneck(m=10):
    # necklace: m beads (K4 minus an edge) joined in a ring through their two degree-2 vertices; 3-regular
    E = []
    for i in range(m):
        off = 4*i
        for a, b in itertools.combinations(range(4), 2):
            if (a, b) != (0, 1): E.append((off+a, off+b))
        E.append((off+1, (4*(i+1)) % (4*m)))
    return 4*m, E
def adj(n, E):
    A = np.zeros((n, n))
    for a, b in E: A[a, b] = A[b, a] = 1
    return A

OUT['graphs'] = {}
for name, (n, E) in [('petersen', petersen()), ('necklace10', bottleneck(10))]:
    A = adj(n, E); deg = A.sum(1)
    assert np.all(deg == 3)
    q = 2
    lam = np.sort(np.linalg.eigvalsh(A))[::-1]
    roots = []
    for l in lam[1:]:
        roots += list(np.roots([q, -l, 1]))     # q u^2 - l u + 1
    mods = np.abs(roots)
    ram = bool(np.all(np.abs(lam[1:]) <= 2*math.sqrt(q)+1e-12))
    print(f"(a) {name}: eigenvalues {np.round(lam, 4).tolist()}; Ramanujan bound 2 sqrt2 = {2*math.sqrt(2):.4f}; "
          f"Ramanujan {ram}; Ihara root moduli range [{mods.min():.4f}, {mods.max():.4f}] vs 1/sqrt2 = {1/math.sqrt(2):.4f}")
    OUT['graphs'][name] = {'eig': lam.tolist(), 'ramanujan': ram, 'roots': [[r.real, r.imag] for r in roots]}

# (b) elliptic curve y^2 = x^3 + x + 1
def count(p):
    sq = {}
    for y in range(p): sq[y*y % p] = sq.get(y*y % p, 0) + 1
    return 1 + sum(sq.get((x**3 + x + 1) % p, 0) for x in range(p))
from sympy import primerange
ap = []
for p in primerange(3, 2000):
    if p == 31: continue
    a = p + 1 - count(p)
    ap.append((p, a, a*a - 4*p))
worst = max(ap, key=lambda r: r[1]**2/(4*r[0]))
print(f"(b) y^2=x^3+x+1, {len(ap)} primes < 2000: all a_p^2 <= 4p: {all(r[2] <= 0 for r in ap)}; "
      f"largest |a_p|/(2 sqrt p) = {abs(worst[1])/(2*math.sqrt(worst[0])):.4f} at p = {worst[0]}")
OUT['elliptic'] = {'ap': ap, 'all_hasse': all(r[2] <= 0 for r in ap)}

# (c) Lee-Yang on the same two graphs
def lee_yang(n, E, a):
    """Z(z) = sum_S a^{|boundary S|} z^{|S|}, a = exp(-2J)."""
    coeff = np.zeros(n+1)
    for mask in range(1 << n):
        b = sum(((mask >> i) & 1) != ((mask >> j) & 1) for i, j in E)
        coeff[bin(mask).count('1')] += a**b
    return np.roots(coeff[::-1])
OUT['lee_yang'] = {}
for name, (n, E) in [('petersen', petersen()), ('necklace3', bottleneck(3))]:
    for J, a in [('ferro J=0.5', math.exp(-1.0)), ('antiferro J=-0.5', math.exp(1.0))]:
        r = lee_yang(n, E, a)
        print(f"(c) {name} {J}: |z| of Lee-Yang zeros in [{np.abs(r).min():.6f}, {np.abs(r).max():.6f}]")
        OUT['lee_yang'][f"{name} {J}"] = [[z.real, z.imag] for z in r]

# (d) the dial reading of zeta's own local factor
print("(d) zeta local factor 1 - u: (k, d) = (-1, 0), Delta = 1 > 0, root at |u| = 1 (Re s = 0): the local test says nothing.")
print("    Ramanujan-Petersson for a weight-2 newform: factor 1 - a_p u + p u^2 with Delta <= 0 (Deligne); its global L-function RH is open.")
json.dump(OUT, open('solved_shapes.json', 'w'))
