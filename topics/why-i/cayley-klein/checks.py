"""Checks for topics/why-i/cayley-klein (Nov 25, 2026). Run: python3 checks.py
Question: is the dial k of the weighted circle operation a o_k b = a + b + k ab the Cayley-Klein dial kappa,
does the plane's second dial d supply the second (signature) dial, and do chains give every Cayley-Klein dial?"""
import numpy as np, itertools
from scipy.linalg import expm
rng = np.random.default_rng(1125)
results = []
def check(name, ok, detail=""):
    results.append((name, bool(ok))); print(("PASS " if ok else "FAIL ") + name + (("  | " + detail) if detail else ""))
circ = lambda a, b, k: a + b + k*a*b
def tanadd(x, y, kap): return (x + y)/(1 - kap*x*y)                 # Cayley-Klein kappa-tangent addition
def C_(t, kap):  # kappa-cosine
    s = np.sqrt(complex(kap)); return np.cos(s*t) if kap != 0 else 1.0 + 0*t
def S_(t, kap):  # kappa-sine
    s = np.sqrt(complex(kap)); return np.sin(s*t)/s if kap != 0 else t + 0j
def kofkappa(kap): return 2*np.sqrt(complex(-kap))                   # [ours] k = 2 sqrt(-kappa), i.e. kappa = -k^2/4
def mob(T, k): return T/(1 - (k/2)*T)                                # Cayley/Moebius map, a = T + O(T^2)
KAPS = [-0.7, -0.25, 0.0, 0.25, 0.5]

# ---------- K. one dial: k is the Cayley-Klein kappa of a one-parameter group ----------
ok = True
for kap in KAPS:
    for _ in range(20):
        t1, t2 = rng.uniform(-.8, .8, 2)
        R = lambda t: np.array([[C_(t, kap), -kap*S_(t, kap)], [S_(t, kap), C_(t, kap)]])
        ok &= np.allclose(R(t1) @ R(t2), R(t1 + t2))
        T = lambda t: S_(t, kap)/C_(t, kap)
        ok &= np.isclose(T(t1 + t2), tanadd(T(t1), T(t2), kap))
check("K1 a Cayley-Klein 'rotation' with label kappa composes its tangent T = S_kappa/C_kappa by (x + y)/(1 - kappa x y): "
      "circle for kappa > 0, line for 0, hyperbola for kappa < 0 (Yaglom; Herranz-Ortega-Santander 2000)", ok)
ok = True
for kap in KAPS:
    k = kofkappa(kap)
    for _ in range(50):
        x, y = rng.uniform(-.5, .5, 2)
        ok &= np.isclose(mob(tanadd(x, y, kap), k), circ(mob(x, k), mob(y, k), k))
check("K2 [ours, a short derivation] the Cayley map a = T/(1 - (k/2) T) turns kappa-tangent addition into a o_k b exactly, with "
      "kappa = -k^2/4: real k gives kappa < 0, k = 0 gives kappa = 0, imaginary k gives kappa > 0", ok)
ok = True
for kap in [-0.7, -0.25, 0.25, 0.5]:
    k = kofkappa(kap)
    for t in rng.uniform(-.6, .6, 20):
        a = mob(S_(t, kap)/C_(t, kap), k)
        ok &= np.isclose(np.log(1 + k*a)/k, t)
check("K3 the hidden coordinate log(1 + k a)/k is the Cayley-Klein arc length (rapidity, angle) theta: the additive coordinate is the same", ok)
ok = True
for kap in [0.25, 0.5]:
    k = kofkappa(kap); T = rng.uniform(-3, 3, 50)
    ok &= np.allclose(np.abs(1 + k*mob(T, k)), 1)
for kap in [-0.25, -0.7]:
    k = kofkappa(kap); T = rng.uniform(-.9, .9, 50)/np.sqrt(-kap)
    ok &= np.all(np.isreal(1 + k*mob(T, k))) and np.all((1 + k*mob(T, k)).real > 0)
check("K4 the lift 1 + k a is a unit phase e^(2 i sqrt(kappa) theta) when kappa > 0 (imaginary k) and a positive number e^(k theta) when "
      "kappa < 0 (real k): our k = -i is the elliptic (rotation) case, k = -1 the hyperbolic (boost) case", ok)
c = 3.0; ok = True
for _ in range(50):
    u, v = rng.uniform(-.95*c, .95*c, 2)
    w = (u + v)/(1 + u*v/c**2)
    ok &= np.isclose(w, tanadd(u, v, -1/c**2))
    k = 2/c; a = lambda s: mob(s, k)
    ok &= np.isclose(a(w), circ(a(u), a(v), k)) and np.isclose(1 + k*a(u), (1 + u/c)/(1 - u/c))
check("K5 velocity addition is the kappa = -1/c^2 law (Levy-Leblond 1976); with k = 2/c the lift 1 + k a is the squared Doppler factor "
      "(Bondi's k-calculus); c = 1 gives a = u/(1 - u), and the repo's a = 2u/(1 - u) is the same law rescaled (k = 1)", ok)

# ---------- D. the plane's second dial d: a coordinate, not a second geometry ----------
F = lambda x, y, k, d: (x + y + k*x*y)/(1 - d*x*y)                  # Coleman-McGuinness rational formal group laws
ok = True
for _ in range(40):
    k, d = rng.uniform(-2, 2, 2); kap = d - k**2/4
    m = lambda y: y/(1 - (k/2)*y)                                     # x = m(y)
    y1, y2 = rng.uniform(-.2, .2, 2)
    ok &= np.isclose(m(tanadd(y1, y2, kap)), F(m(y1), m(y2), k, d))
check("D1 [ours, from two known results] every law (a + b + k ab)/(1 - d ab) is the kappa-tangent law in disguise, with kappa = d - k^2/4 "
      "(the same Cayley map). So (k, d) is two coordinates on ONE dial; lines of equal kappa are the parabolas d = kappa + k^2/4", ok)
ok = True
for _ in range(40):
    k, d = rng.uniform(-2, 2, 2); kap = d - k**2/4
    roots = np.roots([d, k, 1]) if abs(d) > 1e-12 else np.array([-1/k])
    nreal = np.sum(np.abs(np.imag(roots)) < 1e-9)
    ok &= (nreal == 2) == (kap < 0) and (nreal == 0) == (kap > 0)
check("D2 the walls of the (k, d) law (real roots of 1 + k x + d x^2, Coleman-McGuinness's discriminant k^2 - 4d = -4 kappa) are "
      "two for kappa < 0, none for kappa > 0: only the sign of kappa is a real invariant, as in Cayley-Klein", ok)
ok = True
for k, d in [(0.0, -1/9), (2/3, 0.0), (1.0, 0.25 - 1/9), (-2/3, 0.0)]:
    ok &= np.isclose(d - k**2/4, -1/9)
check("D3 the repo's relativity point (k, d) = (0, -1/c^2), pure o_k at k = 2/c, and every point on the parabola d = k^2/4 - 1/c^2 are the "
      "same geometry, kappa = -1/c^2; the repo's line q = k^2/4 - 1 (c = 1) is exactly this parabola", ok)
ok = True
for kap1 in [-0.5, 0.0, 0.5]:
    for kap2 in [-0.5, 0.0, 0.5]:
        P0 = np.array([[0, -kap1, 0], [1, 0, 0], [0, 0, 0]]); K = np.array([[0, 0, 0], [0, 0, -kap2], [0, 1, 0]])
        A, B = expm(0.7*P0), expm(0.9*K)
        ok &= not np.allclose(A @ B, B @ A)
check("D4 all nine Cayley-Klein groups of the plane are non-commutative (even Galilean), while every scalar law (k, d) is commutative: "
      "no scalar dial, k or d, can hold the second Cayley-Klein dial", ok)

# ---------- N. higher dimensions: one label per pair, but only N free dials ----------
def gens(labels, N):
    G = {}
    for (a, b), kab in labels.items():
        M = np.zeros((N + 1, N + 1)); M[b, a] = 1; M[a, b] = -kab; G[(a, b)] = M
    return G
def closure_residual(G):
    basis = np.array([M.ravel() for M in G.values()]).T; worst = 0.0
    for (p, q) in itertools.combinations(G.values(), 2):
        Cm = (p @ q - q @ p).ravel()
        coef, *_ = np.linalg.lstsq(basis, Cm, rcond=None); worst = max(worst, np.abs(basis @ coef - Cm).max())
    return worst
def ck_labels(kaps):
    N = len(kaps); return {(a, b): float(np.prod(kaps[a:b])) for a in range(N + 1) for b in range(a + 1, N + 1)}
ok_ck, bad = True, []
for N in [2, 3, 4]:
    for _ in range(10):
        kaps = rng.choice([-1.3, -0.4, 0.0, 0.6, 1.1], N)
        ok_ck &= closure_residual(gens(ck_labels(kaps), N)) < 1e-10
        free = {p: rng.uniform(-1.5, 1.5) for p in ck_labels(kaps)}
        bad.append(closure_residual(gens(free, N)))
check("N1 in N dimensions the motion generators J_ab (one per pair of directions, label kappa_ab) close into a Lie algebra when "
      "kappa_ac = kappa_ab kappa_bc, i.e. every label is a product of N dials kappa_1..kappa_N (Cayley-Klein, Herranz-Santander); "
      "with an independent label per pair they do not close", ok_ck and min(bad) > 1e-3,
      f"product rule: closes for N = 2, 3, 4; free labels: smallest closure failure {min(bad):.2f}")
ok = True
for N in [2, 3]:
    kaps = np.array([-0.6, 0.8, -0.3][:N]); L = ck_labels(kaps); G = gens(L, N)
    for (a, b), M in G.items():
        k = kofkappa(L[(a, b)])
        for _ in range(5):
            t1, t2 = rng.uniform(-.5, .5, 2)
            R1, R2, R12 = expm(t1*M), expm(t2*M), expm((t1 + t2)*M)
            T = lambda R: R[b, a]/R[a, a]
            ok &= np.isclose(mob(T(R12), k), circ(mob(T(R1), k), mob(T(R2), k), k))
check("N2 [ours] each pair's one-parameter group is a scalar o_k chain with its own k_ab = 2 sqrt(-kappa_ab), so a chain across the planes "
      "carries all N(N+1)/2 labels, but only N of the k's are free: k_ac^2 = -k_ab^2 k_bc^2/4", ok)
ok = True
for N in [2, 3]:
    kaps = rng.choice([-1.0, 0.0, 0.7], N); G = list(gens(ck_labels(kaps), N).values())
    for _ in range(10):
        A = expm(sum(rng.normal()*M for M in G)) - np.eye(N + 1); B = expm(sum(rng.normal()*M for M in G)) - np.eye(N + 1)
        ok &= np.allclose(np.eye(N + 1) + (A + B + A @ B), (np.eye(N + 1) + A) @ (np.eye(N + 1) + B))
check("N3 with matrix elements, o_k (k = 1, A = G - 1) is the group law of every Cayley-Klein group, all dials included: the dials live "
      "in the matrices, not in k", ok)
k = -1.0; ok = np.isclose(-(1j*k)**2/4, -(-k**2/4)) and np.isclose(-(-1j)**2/4, 0.25) and np.isclose(-(-1.0)**2/4, -0.25)
check("N4 multiplying k by i flips the sign of that plane's kappa (hyperbolic <-> elliptic): our k = -1 -> -i step is the Cayley-Klein "
      "analytic continuation of one label (Gromov's j: 1 -> i), the 'Wick rotation' of that plane", ok)
print(f"\n{sum(r[1] for r in results)}/{len(results)} pass")
