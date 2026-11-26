"""DRAFT checks: are the k^3 and k^4 Trotter error terms rungs or corner terms? (Affine Composition Theory)

Run: python3 checks_trotter.py > checks_trotter_output.txt   (numpy only; exact fractions for the algebra)

Set-up. N gates in a row, each 1 + k x_i, compose by our core law:
    1 + k (x_1 o_k ... o_k x_N) = (1 + k x_1)(1 + k x_2) ... (1 + k x_N)
which has one term per corner S of the N-cube (the ordered product of the x_i in S, weight k^|S|).
The generator ('radius') of the result is log of that product. The question: when the x_i do not commute,
what is left over at each order, and is it a new rung (another log) or a sum over corners?
"""
import itertools, math
from fractions import Fraction as F
import numpy as np

results = []
def check(name, ok, detail):
    results.append((name, bool(ok))); print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}")

from trotter_lib import *

print("="*78); print("1. THE LEFTOVER IS A SUM OVER THE CUBE'S ORDERINGS (exact)"); print("="*78)
# multilinear part of log(prod (1 + x_i)) computed straight from the log series, letters 0..N-1 each used once
for N in range(2, 7):
    Y = {}
    for r in range(1, N+1):
        for S in itertools.combinations(range(N), r):
            Y[S] = F(1)                      # corner S of the N-cube: ordered product, increasing order
    logp = {}; Ym = {(): F(1)}
    for m in range(1, N+1):
        Ym = {w: c for w, c in pmul(Ym, Y, N).items() if len(set(w)) == len(w)}
        logp = padd(logp, pscale(Ym, F((-1)**(m+1), m)))
    top = {w: c for w, c in logp.items() if len(w) == N}
    ok = len(top) == math.factorial(N) and all(c == weight(N, descents(w)) for w, c in top.items())
    check(f"E1 N={N}: generator of x_1 o ... o x_{N}, top level = sum over all {math.factorial(N)} orderings, weight (-1)^d / (N C(N-1,d)), d = descents",
          ok, "weights by d: " + ", ".join(f"d={d}: {weight(N, d)}" for d in range(N)))
for N in range(2, 9):
    tot = sum(weight(N, descents(s)) for s in itertools.permutations(range(N)))
    if tot != 0: break
check("E2 the weights add to zero at every level N >= 2: so when the inputs commute, the radius adds exactly (no error)", tot == 0, "N = 2..8")
eul = {N: [sum(1 for s in itertools.permutations(range(N)) if descents(s) == d) for d in range(N)] for N in (3, 4, 5)}
check("E3 how many orderings carry each weight: Eulerian numbers (1,4,1), (1,11,11,1), (1,26,66,26,1) = n! x volumes of the cube's slabs d <= x_1+...+x_n <= d+1",
      eul == {3: [1, 4, 1], 4: [1, 11, 11, 1], 5: [1, 26, 66, 26, 1]}, f"{eul}")
# slab volumes by sampling
rng = np.random.default_rng(1)
pts = rng.random((2_000_000, 4)); s = np.floor(pts.sum(1)).astype(int)
fr = np.bincount(s, minlength=4)[:4]/len(pts)
check("E3b sampling the 4-cube: slab volumes 1/24, 11/24, 11/24, 1/24", np.allclose(fr, np.array([1, 11, 11, 1])/24, atol=1e-3), f"{np.round(fr*24, 3)} / 24")

print(); print("="*78); print("2. THE k^3 AND k^4 TROTTER TERMS FROM THOSE CORNER WEIGHTS"); print("="*78)
A, B = L("A"), L("B")
bch3 = pscale(padd(comm(A, comm(A, B)), comm(B, comm(B, A))), F(1, 12))
bch4 = pscale(comm(B, comm(A, comm(A, B))), F(-1, 24))
m2 = magnus_piecewise([("A", 1), ("B", 1)], 2)
m3 = magnus_piecewise([("A", 1), ("B", 1)], 3)
m4 = magnus_piecewise([("A", 1), ("B", 1)], 4)
check("T1 k^2 term from the corner weights = (1/2)[A,B]", m2 == pscale(comm(A, B), F(1, 2)), f"{m2}")
check("T2 k^3 term from the corner weights = (1/12)([A,[A,B]] + [B,[B,A]])  (exact, as words)", m3 == bch3, f"{len(m3)} words agree")
check("T3 k^4 term from the corner weights = -(1/24)[B,[A,[A,B]]]  (exact, as words)", m4 == bch4, f"{len(m4)} words agree")

# numerical cross-check on SU(2): fit log(e^{kA} e^{kB}) in powers of k
def qmul(a, b):
    w1, x1, y1, z1 = a; w2, x2, y2, z2 = b
    return np.array([w1*w2-x1*x2-y1*y2-z1*z2, w1*x2+x1*w2+y1*z2-z1*y2, w1*y2-x1*z2+y1*w2+z1*x2, w1*z2+x1*y2-y1*x2+z1*w2])
def qexp(h):
    t = np.linalg.norm(h[1:]); return np.r_[math.cos(t), h[1:]/t*math.sin(t)] if t else np.array([1., 0, 0, 0])
def qlog(q):
    v = q[1:]; s = np.linalg.norm(v); return np.r_[0, v/s*math.atan2(s, q[0])]
def evalp(p, env):
    tot = np.zeros(4)
    for w, c in p.items():
        prod = np.array([1., 0, 0, 0])
        for let in w: prod = qmul(prod, env[let])
        tot += float(c)*prod
    return tot
Aq = np.array([0, 0.9, -0.3, 0.4]); Bq = np.array([0, 0.2, 0.7, -0.5])
ks = np.cos(np.pi*(np.arange(40)+0.5)/40)*0.25
vals = np.array([qlog(qmul(qexp(k*Aq), qexp(k*Bq))) for k in ks])
V = np.vander(ks, 12, increasing=True)
coef = np.linalg.lstsq(V, vals, rcond=None)[0]
e3 = np.linalg.norm(coef[3] - evalp(m3, {"A": Aq, "B": Bq})); e4 = np.linalg.norm(coef[4] - evalp(m4, {"A": Aq, "B": Bq}))
check("T4 numerically on qubit gates: the fitted k^3 and k^4 parts of log(e^{kA} e^{kB}) match the corner-weight terms",
      e3 < 1e-8 and e4 < 1e-8, f"k^3 miss {e3:.1e}, k^4 miss {e4:.1e} (sizes {np.linalg.norm(coef[3]):.3f}, {np.linalg.norm(coef[4]):.3f})")

# Strang splitting: e^{kA/2} e^{kB} e^{kA/2}: three time blocks, same weights
s2 = magnus_piecewise([("A", F(1, 2)), ("B", 1), ("A", F(1, 2))], 2)
s3 = magnus_piecewise([("A", F(1, 2)), ("B", 1), ("A", F(1, 2))], 3)
strang3 = padd(pscale(comm(A, comm(A, B)), F(-1, 24)), pscale(comm(B, comm(B, A)), F(1, 12)))
vals = np.array([qlog(qmul(qmul(qexp(k*Aq/2), qexp(k*Bq)), qexp(k*Aq/2))) for k in ks])
coef = np.linalg.lstsq(V, vals, rcond=None)[0]
check("T5 symmetric (Strang) splitting: the same corner weights on three time blocks kill the k^2 term and give the k^3 term -(1/24)[A,[A,B]] + (1/12)[B,[B,A]]",
      s2 == {} and s3 == strang3 and np.linalg.norm(coef[3] - evalp(s3, {"A": Aq, "B": Bq})) < 1e-8,
      f"k^2 term {s2}; k^3 matches as words and numerically ({np.linalg.norm(coef[3] - evalp(s3, {'A': Aq, 'B': Bq})):.1e})")

print(); print("="*78); print("3. RUNGS: does any higher term need another log (a higher rung)?"); print("="*78)
# single-B part of the series: coefficients B_m/m! (Bernoulli), the reciprocal of the corner-simplex series (e^z-1)/z
def single_B(p): return {w: c for w, c in p.items() if w.count("B") == 1}
bern = [F(1), F(1, 2), F(1, 6), F(0), F(-1, 30), F(0), F(1, 42)]   # B_m with B_1 = +1/2
ok = True; shown = []
for n in range(2, 7):
    mn = single_B(magnus_piecewise([("A", 1), ("B", 1)], n))
    m = n - 1
    # ad_A^m (B) as words
    adm = B
    for _ in range(m): adm = comm(A, adm)
    target = pscale(adm, bern[m] / math.factorial(m))
    ok &= (mn == target); shown.append(f"n={n}: {bern[m]/math.factorial(m)}")
check("R1 the part of each term with one B is (B_m/m!) ad_A^m B: Bernoulli numbers, so the k^4 single-B part is 0",
      ok, "; ".join(shown))
# z/(e^z - 1) * (e^z - 1)/z = 1, and (e^z-1)/z = sum z^n/(n+1)! are the corner-simplex volumes 1/(n+1)!
c = [F(1, math.factorial(n+1)) for n in range(8)]                    # (e^z-1)/z
b = [((-1)**m if m == 1 else 1)*bern[m]/math.factorial(m) if m < len(bern) else None for m in range(7)]
b = [F(1), F(-1, 2), F(1, 12), F(0), F(-1, 720), F(0), F(1, 30240)]  # z/(e^z-1)
prod = [sum(b[i]*c[n-i] for i in range(n+1)) for n in range(7)]
check("R2 those Bernoulli weights are the reciprocal series of the corner-simplex volumes 1/(n+1)! (the series that gives e)",
      prod == [1, 0, 0, 0, 0, 0, 0], "(z/(e^z-1)) x ((e^z-1)/z) = 1 to order 6")
# polynomial: every term is a finite word polynomial (no logs inside): check degrees
ok = all(all(len(w) == n for w in magnus_piecewise([("A", 1), ("B", 1)], n)) for n in range(2, 6))
check("R3 every order is a plain polynomial in the gates (words of length n): one log at the start, none inside. Not a new rung.",
      ok, "orders 2..5")


print(); print("="*78); print("POLISH: replace each estimate above with an exact closed form (cheap)"); print("="*78)
import time
t0 = time.perf_counter()
# P1: slab volumes exactly (Irwin-Hall CDF), instead of 2 million random points
def cdf_sum_uniform(n, x):  # P(x_1+...+x_n <= x), exact for integer x
    return sum(F((-1)**j*math.comb(n, j)*(x - j)**n) for j in range(0, x+1))/math.factorial(n)
slabs = {n: [cdf_sum_uniform(n, d+1) - cdf_sum_uniform(n, d) for d in range(n)] for n in (3, 4, 5)}
check("P1 polish E3b: exact slab volumes of the n-cube = Eulerian number / n!, no sampling",
      all(slabs[n] == [F(e, math.factorial(n)) for e in eul[n]] for n in (3, 4, 5)), f"n=4: {[str(v) for v in slabs[4]]}")
# P2: exact k-series of log(e^{kA} e^{kB}) on quaternions, rational arithmetic, instead of a least-squares fit
NDEG = 5
def qmulF(a, b):
    w1, x1, y1, z1 = a; w2, x2, y2, z2 = b
    return (w1*w2-x1*x2-y1*y2-z1*z2, w1*x2+x1*w2+y1*z2-z1*y2, w1*y2-x1*z2+y1*w2+z1*x2, w1*z2+x1*y2-y1*x2+z1*w2)
Z4 = (F(0),)*4; O4 = (F(1), F(0), F(0), F(0))
def qadd(a, b): return tuple(p+q for p, q in zip(a, b))
def qsc(a, c): return tuple(p*c for p in a)
def smul(P, Q):  # truncated series product, coefficients are quaternions
    R = [Z4]*(NDEG+1)
    for i, p in enumerate(P):
        for j, q in enumerate(Q[:NDEG+1-i]): R[i+j] = qadd(R[i+j], qmulF(p, q))
    return R
def sexp(h, c=F(1)):  # e^{k c h} as a series
    out, term = [O4], O4
    for n in range(1, NDEG+1): term = qsc(qmulF(term, h), c/n); out.append(term)
    return out
def slog(P):  # log(1+X), X = P - 1 has no constant term
    X = [Z4] + P[1:]; out = [Z4]*(NDEG+1); Xn = [O4] + [Z4]*NDEG
    for n in range(1, NDEG+1):
        Xn = smul(Xn, X); out = [qadd(o, qsc(x, F((-1)**(n+1), n))) for o, x in zip(out, Xn)]
    return out
def evalF(p, env):
    tot = Z4
    for w, c in p.items():
        prod = O4
        for let in w: prod = qmulF(prod, env[let])
        tot = qadd(tot, qsc(prod, c))
    return tot
AF = tuple(F(str(v)) for v in Aq); BF = tuple(F(str(v)) for v in Bq); envF = {"A": AF, "B": BF}
lg = slog(smul(sexp(AF), sexp(BF)))
check("P2 polish T4: exact k-series of log(e^{kA} e^{kB}) (rational arithmetic, no fit): its k^3 and k^4 parts equal the corner-weight terms with zero error",
      lg[3] == evalF(m3, envF) and lg[4] == evalF(m4, envF) and lg[2] == evalF(m2, envF), "k^2, k^3, k^4 agree exactly")
lgS = slog(smul(smul(sexp(AF, F(1, 2)), sexp(BF)), sexp(AF, F(1, 2))))
check("P3 polish T5: exact k-series for Strang splitting: k^2 part is exactly 0 and the k^3 part equals the corner-weight term with zero error",
      lgS[2] == Z4 and lgS[3] == evalF(s3, envF), "exact")
print(f"POLISH cost: {1000*(time.perf_counter() - t0):.1f} ms in total (the sampling it replaces drew 2,000,000 points)")

print(); n_pass = sum(ok for _, ok in results); print(f"{n_pass} of {len(results)} checks pass")
