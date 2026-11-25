"""Checks for semiring-audit/README.md (DRAFT, Nov 25 2026).

Each check confirms one identification in the audit: that a structure in our
Lean files or write-ups is a known semiring (or provably is not one).
Run: python3 checks.py > checks_output.txt
"""
import itertools, math
import numpy as np
import sympy as sp

rng = np.random.default_rng(7)
results = []

def check(name, ok, detail=""):
    results.append((name, bool(ok)))
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))

a, b, c, k, s = sp.symbols("a b c k s")

# ---------------------------------------------------------------- 1. the shifted pair
sadd = lambda x, y: x + y + 1           # Rungs.sadd, AbstractAlgebra.fadd, PvsNP.sAdd, LinearAlgebra.sadd, FusedPeano.sadd
madd = lambda x, y: x + y + x*y         # multAdd = circle operation at k = 1
check("S1 shift x -> 1+x carries sadd to +", sp.expand(1 + sadd(a, b) - ((1+a) + (1+b))) == 0)
check("S2 shift x -> 1+x carries multAdd to x", sp.expand(1 + madd(a, b) - (1+a)*(1+b)) == 0)
check("S3 multAdd distributes over sadd (madd_sadd / multAdd_fadd_distrib / multAdd_sAdd_distrib are one law)",
      sp.expand(madd(a, sadd(b, c)) - sadd(madd(a, b), madd(a, c))) == 0)
check("S4 zero of the shifted pair is -1 (absorbs under multAdd), one is 0",
      sp.expand(sadd(a, -1) - a) == 0 and sp.expand(madd(a, -1) + 1) == 0 and madd(a, 0) == a)
check("S5 multAdd does NOT distribute over plain + (multAdd_not_distrib_add)",
      madd(1, 1 + 1) != madd(1, 1) + madd(1, 1))

# ---------------------------------------------------------------- 2. general k: a generated (Pap) semiring
g = lambda x: 1 + k*x
ginv = lambda y: (y - 1)/k
oplus_k = lambda x, y: x + y + 1/k
circ_k = lambda x, y: x + y + k*x*y
check("G1 circ_k = g^-1(g(a) g(b)) with generator g(x) = 1 + kx", sp.simplify(ginv(g(a)*g(b)) - circ_k(a, b)) == 0)
check("G2 partner a + b + 1/k = g^-1(g(a) + g(b))", sp.simplify(ginv(g(a) + g(b)) - oplus_k(a, b)) == 0)
check("G3 circ_k distributes over the partner at every k (symbolic)",
      sp.simplify(circ_k(a, oplus_k(b, c)) - oplus_k(circ_k(a, b), circ_k(a, c))) == 0)
check("G4 at k = 1 the partner is exactly rung 1, sadd", sp.simplify(oplus_k(a, b).subs(k, 1) - sadd(a, b)) == 0)

# ---------------------------------------------------------------- 3. height coordinate: log semiring
E = lambda x: sp.exp(x) - 1
Einv = lambda y: sp.log(1 + y)
lse = lambda x, y: sp.log(sp.exp(x) + sp.exp(y))
check("H1 pulled back by E, multAdd is + : E^-1(E a o1 E b) = a + b",
      sp.simplify(sp.exp(Einv(madd(E(a), E(b)))) - sp.exp(a + b)) == 0)
check("H2 pulled back by E, sadd is log-sum-exp : E^-1(E a + E b + 1) = log(e^a + e^b)",
      sp.simplify(sp.exp(Einv(sadd(E(a), E(b)))) - sp.exp(lse(a, b))) == 0)
# numeric version over random reals
X = rng.normal(size=(3, 2000))
nE = np.expm1; nEi = np.log1p
lhs = nEi(nE(X[0]) + nE(X[1]) + 1); rhs = np.logaddexp(X[0], X[1])
check("H3 same, 2000 random reals", np.max(np.abs(lhs - rhs)) < 1e-10, f"max err {np.max(np.abs(lhs-rhs)):.1e}")
d = np.max(np.abs((X[0] + np.logaddexp(X[1], X[2])) - np.logaddexp(X[0]+X[1], X[0]+X[2])))
check("H4 log semiring law a + lse(b,c) = lse(a+b, a+c)", d < 1e-12, f"max err {d:.1e}")
h_vals = [1, 0.1, 0.01, 0.001]
lim = [h*np.logaddexp(2/h, 1.5/h) for h in h_vals]
check("H5 Maslov dequantization: h*lse(a/h, b/h) -> max(a,b) as h -> 0", abs(lim[-1] - 2) < 1e-6, str([round(v, 6) for v in lim]))

# ---------------------------------------------------------------- 4. ladders: Bennett (exp) vs ours (E)
def num_E(x): return np.expm1(x)
def num_Ei(y): return np.log1p(y)
def lift(op, f, fi): return lambda x, y: f(op(fi(x), fi(y)))
OPS = {}
OPS["lse"] = lambda x, y: np.logaddexp(x, y)
OPS["+"] = lambda x, y: x + y
OPS["x"] = lambda x, y: x * y
OPS["a+b+1"] = lambda x, y: x + y + 1
OPS["o1"] = lambda x, y: x + y + x*y
OPS["+(2)"] = lift(OPS["o1"], num_E, num_Ei)          # addition at height 2 = E_+(o1)
OPS["P2"] = lift(OPS["a+b+1"], num_E, num_Ei)         # partner at height 2 = E_+(a+b+1)
OPS["B3"] = lambda x, y: np.exp(np.log(x) * np.log(y))  # Bennett grade 3 on x > 1

def dist_err(mul, add, lo, hi, n=3000):
    A, B, C = rng.uniform(lo, hi, size=(3, n))
    l = mul(A, add(B, C)); r = add(mul(A, B), mul(A, C))
    m = np.isfinite(l) & np.isfinite(r)
    return float(np.max(np.abs(l[m] - r[m]) / (1 + np.abs(r[m]))))

e = dist_err(OPS["+"], OPS["lse"], -2, 2); check("L1 Bennett: + distributes over lse", e < 1e-10, f"{e:.1e}")
e = dist_err(OPS["x"], OPS["+"], -2, 2); check("L2 Bennett: x distributes over +", e < 1e-10, f"{e:.1e}")
e = dist_err(OPS["B3"], OPS["x"], 1.1, 3); check("L3 Bennett: a^(log b) distributes over x (x > 1)", e < 1e-10, f"{e:.1e}")
e = dist_err(OPS["o1"], OPS["+"], 0.1, 1); check("L4 ours: +(1) = o1 does NOT distribute over +(0) = +", e > 1e-2, f"{e:.2f}")
e = dist_err(OPS["+(2)"], OPS["o1"], 0.05, 0.6); check("L5 ours: +(2) does NOT distribute over +(1)", e > 1e-3, f"{e:.3f}")
e = dist_err(OPS["+(2)"], OPS["P2"], 0.02, 0.4); check("L6 ours: +(2) DOES distribute over its partner E_+(a+b+1)", e < 1e-9, f"{e:.1e}")

# ---------------------------------------------------------------- 5. rung 3 is not a semiring product
fexp = lambda x, y: (x + 1)**(y + 1) - 1
A_, B_, C_ = 2, 3, 1
check("R1 fexp right-distributes over multAdd: (a o b)^c = a^c o b^c (fexp_multAdd_right)",
      fexp(madd(A_, B_), C_) == madd(fexp(A_, C_), fexp(B_, C_)))
check("R2 fexp does not left-distribute over multAdd and is not commutative",
      fexp(A_, madd(B_, C_)) != madd(fexp(A_, B_), fexp(A_, C_)) and fexp(2, 3) != fexp(3, 2))

# ---------------------------------------------------------------- 6. witness fold and the k-sum
def foldW(kk, xs):
    r = 0
    for x in reversed(xs): r = x + r + kk*x*r
    return r
ok = True
for n in range(0, 7):
    for bits in itertools.product([0, 1], repeat=n):
        ok &= foldW(-1, bits) == (1 if any(bits) else 0)
        ok &= foldW(0, bits) == sum(bits)
        ok &= foldW(-2, bits) == (sum(bits) % 2)
check("W1 fold at k=-1 is OR (Boolean semiring), k=0 is count (N), k=-2 is XOR (GF(2)), all lists up to 6", ok)
gap = lambda kk, A, B, C: (A*C + B*C + kk*(C*A)*(C*B)) - C*(A + B + kk*A*B)
ok = all(gap(kk, A, B, C) == kk*A*B*C*(C - 1) for kk in range(-3, 3) for A in range(-2, 3) for B in range(-2, 3) for C in range(-2, 3))
check("W2 distributivity gap of AND over the k-sum is k*a*b*c*(c-1) (ksum_distrib_gap)", ok)
ok = all(gap(kk, A, B, C) == 0 for kk in (-1, -2) for A in (0, 1) for B in (0, 1) for C in (0, 1))
check("W3 so on 0/1 values AND distributes over OR and over XOR", ok)

# ---------------------------------------------------------------- 7. power set: sum-product and union convolution
n = 4
t = rng.uniform(-0.5, 1.5, n); u = rng.uniform(-0.5, 1.5, n)
subsets = list(itertools.product([0, 1], repeat=n))
W = lambda dial, A: np.prod([dial[i] for i in range(n) if A[i]])
check("P1 sum over subsets of prod of dials = prod(1 + t) (distributive law, factorization)",
      abs(sum(W(t, A) for A in subsets) - np.prod(1 + t)) < 1e-12)
conv = {A: 0.0 for A in subsets}
for A in subsets:
    for B in subsets:
        conv[tuple(x | y for x, y in zip(A, B))] += W(t, A) * W(u, B)
merged = t + u + t*u
err = max(abs(conv[A] - W(merged, A)) for A in subsets)
check("P2 union (OR) convolution of two dialled power sets = dials multAdded (conv_W)", err < 1e-12, f"{err:.1e}")
zeta = lambda f: {S: sum(f[A] for A in subsets if all(x <= y for x, y in zip(A, S))) for S in subsets}
zt = zeta({A: W(t, A) for A in subsets}); zu = zeta({A: W(u, A) for A in subsets}); zc = zeta(conv)
err = max(abs(zc[S] - zt[S]*zu[S]) for S in subsets)
check("P3 the zeta transform turns union convolution into pointwise product (Fourier meets Mobius)", err < 1e-12, f"{err:.1e}")

# ---------------------------------------------------------------- 8. Sugeno measure is a pseudo-additive measure
lam = 0.7; dens = rng.uniform(0.05, 0.3, 5)
mu = lambda idx: (np.prod([1 + lam*dens[i] for i in idx]) - 1)/lam
A_i, B_i = [0, 2], [1, 3, 4]
check("M1 Sugeno lambda-measure: mu(A u B) = mu(A) o_lambda mu(B) for disjoint A, B (lam_append)",
      abs(mu(A_i + B_i) - (mu(A_i) + mu(B_i) + lam*mu(A_i)*mu(B_i))) < 1e-12)
check("M2 and g(mu) = 1 + lambda*mu is multiplicative: mu is a g-decomposable measure",
      abs((1 + lam*mu(A_i + B_i)) - (1 + lam*mu(A_i))*(1 + lam*mu(B_i))) < 1e-12)

# ---------------------------------------------------------------- 9. min-plus and the soliton tau function
ok = all(A + min(B, C) == min(A + B, A + C) for A in range(6) for B in range(6) for C in range(6))
check("T1 min-plus distributivity on Nat (minplus_distrib)", ok)
heights = rng.normal(size=(4, 1000)) * 5
gap4 = np.log(np.exp(heights).sum(0)) - heights.max(0)
check("T2 log tau = lse of corner heights sits within log 4 of their max (solitons-rungs check 3a)",
      gap4.min() >= 0 and gap4.max() <= math.log(4) + 1e-12)

# ---------------------------------------------------------------- 10. min cut, softened
cuts = rng.uniform(0, 3, size=(6, 500))
for h in (1.0, 0.1, 0.01):
    soft = -h*np.log(np.exp(-cuts/h).sum(0))
    err = np.max(np.abs(soft - cuts.min(0)))
check("C1 soft min cut -h*lse(-c/h) -> min cut as h -> 0 (exploration part 2b)", err < 0.02, f"err at h=0.01: {err:.4f}")

print()
print(f"{sum(ok for _, ok in results)}/{len(results)} pass")
