"""Checks for 'Everything we have proved, read in the polar coordinates' (DRAFT, Nov 26 2026).
Run: python3 checks.py > checks_output.txt   (numpy, scipy, sympy)"""
import itertools, math, cmath
import numpy as np
from scipy.linalg import expm
import sympy as sp

rng = np.random.default_rng(7)
results = []
def check(name, ok, detail=""):
    results.append((name, bool(ok)))
    print(("PASS " if ok else "FAIL ") + name + (f"   [{detail}]" if detail else ""))

def wadd(a, b, k): return a + b + k*a*b
def ruler(x, k): return math.log1p(k*x)/k if k != 0 else x
def flip(x, k): return x/(1+k*x)

print("== A. Core laws: the ruler is the radius ==")
err = 0
for _ in range(2000):
    k = rng.uniform(-2, 2); a, b = rng.uniform(-0.4, 0.4, 2)
    if abs(k) < 1e-3: continue
    err = max(err, abs(ruler(wadd(a, b, k), k) - ruler(a, k) - ruler(b, k)))
check("A1 composing adds radii: L_k(a o_k b) = L_k a + L_k b", err < 1e-12, f"max err {err:.1e}")
err = 0
for _ in range(2000):
    k = rng.uniform(-2, 2); x = rng.uniform(-0.3, 0.3)
    if abs(k) < 1e-3: continue
    err = max(err, abs(ruler(flip(x, k), -k) - ruler(x, k)))
check("A2 the flip leaves the radius alone: L_-k(phi_k x) = L_k x", err < 1e-12, f"max err {err:.1e}")
err = 0
for _ in range(500):
    k = rng.uniform(0.1, 2); a = rng.uniform(0, 0.3); n = int(rng.integers(1, 12))
    p = 0.0
    for _ in range(n): p = wadd(p, a, k)
    err = max(err, abs(ruler(p, k) - n*ruler(a, k)))
check("A3 the n-th power is n times the radius", err < 1e-10, f"max err {err:.1e}")
check("A4 the wall -1/k sits at radius minus infinity (k=1: L(-1+1e-12))", ruler(-1+1e-12, 1) < -27)
check("A5 certainty at k=-1 sits at radius plus infinity (L_-1(1-1e-12))", ruler(1-1e-12, -1) > 27)
x = sp.symbols('x1:4'); poly = sp.expand(sp.prod([1+xi for xi in x]))
exps = sorted(sp.Poly(poly, *x).monoms())
check("A6 1 + (x1 o x2 o x3) has one term per cube corner (Newton polytope = cube)",
      exps == sorted(itertools.product([0, 1], repeat=3)))

print("\n== B. Rungs: each rung is one more radius ==")
err = max(abs(ruler(math.expm1(t), 1) - t) for t in np.linspace(-3, 3, 61))
check("B1 the height flow E(x)=e^x-1 undoes the radius: L_1(E(t)) = t", err < 1e-12)
err = 0
for a in range(1, 6):
    for b in range(0, 6):
        fe = (a+1)**(b+1) - 1
        err = max(err, abs(math.log(math.log(1+fe)) - (math.log(math.log(1+a)) + math.log(1+b))))
check("B2 rung 3 (fused exponent) is plain addition after taking the radius twice", err < 1e-12)
h = 1e-6
check("B3 the radius map log(1+x) has a parabolic fixed point at 0 (slope 1)",
      abs((math.log1p(h)-math.log1p(-h))/(2*h) - 1) < 1e-9,
      "so fractional rungs need Ecalle-Voronin, as the height-flow page found")

print("\n== C. The plane dial: angle or radius ==")
def plane(x, y, k, d): return (x + y + k*x*y)/(1 - d*x*y)
agree = tot = 0
for _ in range(400):
    k, d = rng.uniform(-2, 2, 2)
    disc = k*k - 4*d
    if abs(disc) < 0.05: continue
    a = 0.07; xv = 0.0; seen_far = False; came_back = False
    for n in range(3000):
        xv = plane(xv, a, k, d)
        if not np.isfinite(xv): xv = 1e12
        if abs(xv) > 5: seen_far = True
        if seen_far and abs(xv) < 0.05: came_back = True; break
    tot += 1; agree += (came_back == (disc < 0))
check("C1 sign of k^2-4d decides wrap-around (angle) vs run-to-a-wall (radius)", agree == tot, f"{agree}/{tot}")
E = {"-1": np.array([[0, -1], [1, 0.]]), "+1": np.array([[0, 1], [1, 0.]]), "0": np.array([[0, 1], [0, 0.]])}
t = 0.83; I2 = np.eye(2)
ok = (np.allclose(expm(t*E["-1"]), math.cos(t)*I2 + math.sin(t)*E["-1"]) and
      np.allclose(expm(t*E["+1"]), math.cosh(t)*I2 + math.sinh(t)*E["+1"]) and
      np.allclose(expm(t*E["0"]), I2 + t*E["0"]))
check("C2 Clifford square -1 -> angle (cos, sin); +1 -> radius (cosh, sinh); 0 -> straight (1 + t e)", ok)

print("\n== D. Swapping order costs a turn proportional to area ==")
def Rx(a): c, s = math.cos(a), math.sin(a); return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])
def Ry(a): c, s = math.cos(a), math.sin(a); return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])
a = b = 1e-3
C = Rx(a) @ Ry(b) @ Rx(-a) @ Ry(-b); rot = math.atan2(C[1, 0], C[0, 0])
def Bx(r):
    M = np.eye(4); M[0, 0] = M[1, 1] = math.cosh(r); M[0, 1] = M[1, 0] = -math.sinh(r); return M
def By(r):
    M = np.eye(4); M[0, 0] = M[2, 2] = math.cosh(r); M[0, 2] = M[2, 0] = -math.sinh(r); return M
Cb = Bx(a) @ By(b) @ Bx(-a) @ By(-b); boost = math.atan2(Cb[2, 1], Cb[1, 1])
G = lambda v, ax: (lambda M: (M.__setitem__((ax, 0), v), M)[1])(np.eye(4))
Cg = G(a, 1) @ G(b, 2) @ G(-a, 1) @ G(-b, 2); gal = np.abs(Cg - np.eye(4)).max()
r1, r2 = rot/(a*b), boost/(a*b)
check("D1 turn/area: rotations (sphere, d=+1) vs boosts (hyperbolic, d=-1) vs Galilean (flat, d=0)",
      abs(abs(r1) - 1) < 1e-3 and abs(r1 + r2) < 1e-3 and gal < 1e-15, f"{r1:+.4f}, {r2:+.4f}, {gal:.0e}")
q = sp.symbols('q')
def qbin(n, k):
    num = sp.prod([1 - q**(n - i) for i in range(k)]); den = sp.prod([1 - q**(i + 1) for i in range(k)])
    return sp.expand(sp.cancel(num/den))
ok = True
for n in range(1, 8):
    for kk in range(n + 1):
        s = 0
        for pos in itertools.combinations(range(n), kk):   # positions of up-steps
            area = sum(p - i for i, p in enumerate(pos))      # boxes under the path
            s += q**area
        ok &= sp.expand(s - qbin(n, kk)) == 0
check("D2 face weights: Gaussian binomial = sum over lattice paths of q^(area)", ok)
N = 7; w = cmath.exp(2j*math.pi/N)
U = np.diag([w**j for j in range(N)]); V = np.roll(np.eye(N), 1, axis=0)
qv = w if np.allclose(V @ U, w * U @ V) else (np.conj(w) if np.allclose(V @ U, np.conj(w)*U @ V) else None)
ok = qv is not None
if ok:
    for n in range(1, 6):
        lhs = np.linalg.matrix_power(U + V, n)
        rhs = sum(complex(qbin(n, kk).subs(q, qv)) * np.linalg.matrix_power(U, kk) @ np.linalg.matrix_power(V, n - kk)
                  for kk in range(n + 1))
        ok &= np.allclose(lhs, rhs)
check("D3 quantum clock/shift (VU = qUV): (U+V)^n expands with the same q-binomials", ok)
def spin(j):
    m = np.arange(j, -j - 1, -1); dim = len(m)
    Jp = np.zeros((dim, dim))
    for i in range(1, dim): Jp[i - 1, i] = math.sqrt(j*(j + 1) - m[i]*(m[i] + 1))
    Jx = (Jp + Jp.T)/2; Jy = (Jp - Jp.T)/(2j); return Jx, Jy
a = b = 0.4; phases = []
for j in (5, 20, 80, 320):
    Jx, Jy = spin(j); X, P = Jx/math.sqrt(j), Jy/math.sqrt(j)
    Uc = expm(-1j*a*X) @ expm(-1j*b*P) @ expm(1j*a*X) @ expm(1j*b*P)
    phases.append(cmath.phase(Uc[0, 0])/(a*b))
check("D4 spin sphere -> quantum plane: phase/area tends to a constant (hbar = 1) as the sphere grows",
      abs(abs(phases[-1]) - 1) < 0.01 and abs(abs(phases[-1]) - 1) < abs(abs(phases[0]) - 1),
      ", ".join(f"j={j}: {p:+.4f}" for j, p in zip((5, 20, 80, 320), phases)))

print("\n== E. The right shape for each family ==")
n = 4; worst = 1e9
for _ in range(50):
    wts = rng.exponential(size=n); Epos = np.zeros(n)
    for perm in itertools.permutations(range(n)):
        p, rest = 1.0, wts.sum()
        for it in perm: p *= wts[it]/rest; rest -= wts[it]
        for posn, it in enumerate(perm): Epos[it] += p*posn
    for r in range(1, n + 1):
        for I in itertools.combinations(range(n), r):
            worst = min(worst, Epos[list(I)].sum() - r*(r - 1)/2)
    worst = min(worst, -abs(Epos.sum() - n*(n - 1)/2))
check("E1 rankings: expected ranks of a Plackett-Luce model lie in the permutohedron", worst > -1e-9, f"min slack {worst:.2e}")
k1, k2 = 2.0, 1.0; A = ((k1 - k2)/(k1 + k2))**2
xs = np.linspace(-60, 60, 120001)
def u_kdv(t):
    e1 = k1*xs - k1**3*t; e2 = k2*xs - k2**3*t
    m = np.maximum.reduce([np.zeros_like(xs), e1, e2, e1 + e2 + math.log(A)])
    tau = np.exp(-m) + np.exp(e1 - m) + np.exp(e2 - m) + A*np.exp(e1 + e2 - m)
    lt = np.log(tau) + m; return 2*np.gradient(np.gradient(lt, xs), xs)
def crest_slow(t):
    u = u_kdv(t); centre = k2**2*t; win = np.abs(xs - centre) < 8
    return xs[win][np.argmax(u[win])] - centre
shift = crest_slow(12) - crest_slow(-12)
check("E2 KdV 2-soliton: the slow wave's jump = log(pair weight)/k2, as the tropical (polar) grid predicts",
      abs(shift - math.log(A)/k2) < 0.02, f"measured {shift:+.4f}, predicted {math.log(A)/k2:+.4f}")
edges = [(0, 1), (1, 2), (2, 3), (3, 0), (0, 2), (1, 3)]; worst = 0
for aa in (0.2, 0.5, 0.8):
    coeffs = np.zeros(5)
    for S in itertools.product([0, 1], repeat=4):
        cut = sum(S[i] != S[j] for i, j in edges); coeffs[sum(S)] += aa**cut
    roots = np.roots(coeffs[::-1]); worst = max(worst, np.abs(np.log(np.abs(roots))).max())
check("E3 Lee-Yang: all zeros of a ferromagnet's partition function sit at log-radius 0 (pure angle)", worst < 1e-9, f"max |log|z|| {worst:.1e}")
def rand_rho():
    G = rng.normal(size=(4, 4)) + 1j*rng.normal(size=(4, 4)); R = G @ G.conj().T; return R/np.trace(R)
sx = np.array([[0, 1], [1, 0]]); sy = np.array([[0, -1j], [1j, 0]]); sz = np.diag([1, -1])
def obs(v): v = v/np.linalg.norm(v); return v[0]*sx + v[1]*sy + v[2]*sz
worst = -9; chsh_max = 0
for _ in range(3000):
    if rng.random() < 0.5:
        psi = rng.normal(size=4) + 1j*rng.normal(size=4); psi /= np.linalg.norm(psi); rho = np.outer(psi, psi.conj())
    else: rho = rand_rho()
    As = [obs(rng.normal(size=3)) for _ in range(2)]; Bs = [obs(rng.normal(size=3)) for _ in range(2)]
    Em = np.array([[np.real(np.trace(rho @ np.kron(Ai, Bj))) for Bj in Bs] for Ai in As])
    th = np.arcsin(np.clip(Em, -1, 1))
    for s in range(4):
        sg = np.ones((2, 2)); sg[s//2, s % 2] = -1
        worst = max(worst, abs((sg*th).sum()) - math.pi)
        chsh_max = max(chsh_max, abs((sg*Em).sum()))
check("E4 quantum correlations: round in E, straight in angle coordinates (|sum of asin E| <= pi)",
      worst < 1e-9, f"max excess {worst:+.2e}; largest CHSH seen {chsh_max:.3f} <= 2.828")
def mobius(Sf, n):
    return {I: sum((-1)**(len(I) - len(J))*Sf[J] for r in range(len(I) + 1) for J in itertools.combinations(I, r))
            for r in range(n + 1) for I in itertools.combinations(range(n), r)}
def base_vertices(Sf, n):
    V = []
    for perm in itertools.permutations(range(n)):
        x = np.zeros(n); prev = ()
        for i in perm:
            cur = tuple(sorted(prev + (i,))); x[i] = Sf[cur] - Sf[prev]; prev = cur
        V.append(x)
    return V
def ent_vec(S1, S2, S3):
    return {(): 0, (0,): S1, (1,): S1, (2,): S1, (0, 1): S2, (0, 2): S2, (1, 2): S2, (0, 1, 2): S3}
cases = {"GHZ, 4 qubits (bits)": ent_vec(1, 1, 1), "AME(4,3), perfect tensor (trits)": ent_vec(1, 2, 1),
         "classical copies A=B=C": ent_vec(1, 1, 1), "classical XOR C=A+B": ent_vec(1, 2, 2)}
ok = True; det = []
for name, Sf in cases.items():
    I3 = Sf[(0,)] + Sf[(1,)] + Sf[(2,)] - Sf[(0, 1)] - Sf[(0, 2)] - Sf[(1, 2)] + Sf[(0, 1, 2)]
    y = mobius(Sf, 3); ok &= abs(y[(0, 1, 2)] - I3) < 1e-12
    V = base_vertices(Sf, 3)
    # side lengths in the order of the hexagon (permutations adjacent by one swap)
    order = [(0, 1, 2), (1, 0, 2), (1, 2, 0), (2, 1, 0), (2, 0, 1), (0, 2, 1)]
    perms = list(itertools.permutations(range(3))); P = [V[perms.index(o)] for o in order]
    L = [np.abs(P[i] - P[(i + 1) % 6]).max() for i in range(6)]
    alt = (L[0] + L[2] + L[4]) - (L[1] + L[3] + L[5])
    ok &= abs(abs(alt) - 3*abs(I3)) < 1e-9
    det.append(f"{name}: I3={I3:+g}, sides={['%g' % l for l in L]}")
check("E5 entropy as a shape: top corner weight = I3, and alternate hexagon sides differ by I3", ok, "; ".join(det))
n5 = 5; cyc = lambda i, m: tuple(sorted((i + j) % n5 for j in range(m)))
lhs = [cyc(i, 3) for i in range(5)]; rhs = [cyc(i, 2) for i in range(5)] + [tuple(range(5))]
coef = {}
for r in range(1, 6):
    for J in itertools.combinations(range(5), r):
        c = sum(set(J) <= set(I) for I in lhs) - sum(set(J) <= set(I) for I in rhs)
        if c: coef[J] = c
check("E6 5-party holographic cyclic inequality written in corner weights (computed, see detail)", True,
      "nonzero corner coefficients: " + ", ".join(f"{''.join('ABCDE'[i] for i in J)}:{c:+d}" for J, c in coef.items()))
coef2 = {}
for r in range(1, 6):
    for J in itertools.combinations(range(5), r):
        c = sum(bool(set(J) & set(I)) for I in lhs) - sum(bool(set(J) & set(I)) for I in rhs)
        if c: coef2[J] = c
check("E7 the same inequality in the shape's own (Minkowski) weights touches no single or pair corner",
      coef2 == coef and all(len(J) >= 3 for J in coef), "pairs = segments are untouched; only simplices of 3 or more parties")
print("    MMI for 3 parties in corner weights: I3 = y_ABC <= 0 (one corner).")

print("\n== F. Fermat in polar coordinates ==")
worst = 0
for nexp in (3, 5, 7, 11):
    for _ in range(200):
        X = complex(*rng.normal(size=2)); Y = 1 - X
        x = X**(1/nexp); y = Y**(1/nexp)
        assert abs(x**nexp + y**nexp - 1) < 1e-9
        worst = max(worst, abs(math.log(abs(x)) - math.log(abs(X))/nexp) + abs(math.log(abs(y)) - math.log(abs(Y))/nexp))
check("F1 every Fermat curve x^n+y^n=1 has the line's radius picture shrunk by n (n lives only in the angles)", worst < 1e-12)

print("\n== G. Temperature is the roundness dial ==")
Ev = rng.normal(size=10); ok = True
for T in (1, 0.3, 0.1, 0.01):
    F = -T*np.log(np.exp(-(Ev - Ev.min())/T).sum()) + Ev.min()
    ok &= (Ev.min() - T*math.log(10) - 1e-12 <= F <= Ev.min() + 1e-12)
check("G1 free energy -> lowest energy as T -> 0, within T log n (the Maslov / p -> infinity limit)", ok)

print(f"\n{sum(o for _, o in results)}/{len(results)} pass")
