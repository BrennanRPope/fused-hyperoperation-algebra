"""Checks for topics/why-i/reflection-positivity.  Run: python3 checks.py
Each check prints PASS or FAIL with the numbers behind it."""
import itertools
import os
import sys

import numpy as np


sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rp import slice_multiplier, gram, rp_violation, complexify, J_of, expm  # noqa: E402

rng = np.random.default_rng(20261127)
results = []


def check(name, ok, detail):
    results.append(ok)
    print(f"{'PASS' if ok else 'FAIL'}  {name}: {detail}")


def rand_sym(n, eigs):
    Q, _ = np.linalg.qr(rng.normal(size=(n, n)))
    return Q @ np.diag(eigs) @ Q.T


TOL = 1e-10
EIGS = [-1.0, 0.5, 2.0]
H = rand_sym(3, EIGS)
eps = 0.25
k_lo, k_hi = -1 / (eps * max(EIGS)), 1 / (eps * -min(EIGS))   # explicit-step window [-2, 4]

# R1  the chain's multiplier is the product of slice multipliers (matrices, order kept)
k = -0.7
A = [rand_sym(3, rng.normal(size=3)) * 0.3 for _ in range(4)]
chain = A[0]
for a in A[1:]:
    chain = a + chain + k * a @ chain               # a o_k chain, newest step on the left
lhs = np.eye(3) + k * chain
rhs = np.eye(3)
for a in A:
    rhs = (np.eye(3) + k * a) @ rhs
check("R1 chain multiplier = product of slice multipliers", np.abs(lhs - rhs).max() < TOL,
      f"max diff {np.abs(lhs - rhs).max():.1e}")

# R2  the sum over all paths of the product of slice entries is the matrix power
T = slice_multiplier(-1.0, eps, H)
N = 4
tot = np.zeros((3, 3))
for path in itertools.product(range(3), repeat=N + 1):
    w = 1.0
    for a, b in zip(path[:-1], path[1:]):
        w *= T[b, a]
    tot[path[-1], path[0]] += w
check("R2 sum over paths = T^N", np.abs(tot - np.linalg.matrix_power(T, N)).max() < TOL,
      f"max diff {np.abs(tot - np.linalg.matrix_power(T, N)).max():.1e}")

# R3  site reflection: RP holds at every real dial, including where T has negative eigenvalues
ks = [-6, -3, -2, -1, -0.3, 0.4, 2, 4, 7]
v = [rp_violation(gram(slice_multiplier(k, eps, H), 4)) for k in ks]
check("R3 site RP holds for every real k", max(v) < TOL,
      f"worst violation {max(v):.1e} over k in {ks}")

# R4  link reflection: RP holds exactly when 1 + k eps E >= 0 for every energy E
inside = [-1.99, -1, 0, 1, 3.99]
outside = [-3, -2.05, 4.05, 6]
vin = [rp_violation(gram(slice_multiplier(k, eps, H), 5)) for k in inside]
vout = [rp_violation(gram(slice_multiplier(k, eps, H), 5)) for k in outside]
check("R4 link RP holds iff every slice multiplier is >= 0",
      max(vin) < TOL and min(vout) > 1e-4,
      f"window [{k_lo:g}, {k_hi:g}]; inside worst {max(vin):.1e}; outside smallest {min(vout):.1e}")

# R5  link RP <=> T positive semidefinite: the Gram matrix's lowest eigenvalue has the sign
#     of T's lowest eigenvalue, across a sweep of real k
sweep = np.linspace(-3, 5, 33)
agree = True
for k in sweep:
    Tk = slice_multiplier(k, eps, H)
    G = gram(Tk, 3)
    gmin = np.linalg.eigvalsh((G + G.conj().T) / 2).min() / np.abs(G).max()
    tmin = np.linalg.eigvalsh(Tk).min()
    if abs(tmin) > 1e-9 and (gmin < -1e-12) != (tmin < 0):
        agree = False
check("R5 sign of lowest Gram eigenvalue = sign of lowest eigenvalue of T", agree,
      f"{len(sweep)} dial values from -3 to 5")

# R6  off the real axis RP fails for every reflection and every nonzero H
offs = [1j, -1j, -1 + 0.05j, -1 + 0.3j, 0.5 - 0.2j, -2j]
v6 = []
for z in offs:
    for N_ in (4, 5):
        v6.append(rp_violation(gram(slice_multiplier(z, eps, H), N_)))
check("R6 off the real axis RP fails (site and link)", min(v6) > 1e-3,
      f"smallest violation {min(v6):.2e} over k in {offs}")

# R7  why: off the axis the slice multiplier is not Hermitian, so even a single pair of
#     configurations gives a pairing <theta F, F> that is not a real number
z = -1 + 0.3j
Tz = slice_multiplier(z, eps, H)
pair = Tz[0, 1] * Tz[1, 0]                         # <theta F, F> for F = [x1 = 0, x2 = 1], N = 2
check("R7 off-axis pairing is not even real", abs(pair.imag) > 1e-3 and
      np.abs(Tz - Tz.conj().T).max() > 1e-3,
      f"<theta F,F> = {pair:.4f};  |T - T*| = {np.abs(Tz - Tz.conj().T).max():.3f}")

# R8  rewind slices (1 - k eps H)^-1: link RP for every step size iff k E <= 0 for every E
Hpos = rand_sym(3, [0.2, 1.0, 3.0])
Hneg = rand_sym(3, [-1.0, 1.0, 3.0])
eps_list = [0.01, 0.1, 0.5, 1, 2, 5, 20]
vpos = [rp_violation(gram(slice_multiplier(-1.0, e, Hpos, "rewind"), 5)) for e in eps_list]
vneg = [rp_violation(gram(slice_multiplier(-1.0, e, Hneg, "rewind"), 5)) for e in eps_list]
check("R8 rewind slices: RP at every step size iff H >= 0 (at k = -1)",
      max(vpos) < TOL and vneg[0] < TOL and max(vneg) > 1e-3,
      f"H>=0 worst {max(vpos):.1e}; H with E=-1: fine at small steps, "
      f"fails from eps={eps_list[int(np.argmax(np.array(vneg) > 1e-3))]} (wall at eps=1)")

# R9  continuum: e^{k t H} is PSD at every real k, so the lattice window is a slicing artefact
v9 = []
for k in (-10, -1, 1, 10):
    Tc = expm(k * eps * H)
    v9.append(rp_violation(gram(Tc, 5)))
check("R9 continuum slices e^(k eps H): link RP at every real k", max(v9) < TOL,
      f"worst {max(v9):.1e}")

# R10 reversible Markov chains (k = -1, the probability dial) are RP; a one-way cycle is not
pi = np.array([0.5, 0.3, 0.2])
S = rng.uniform(0.1, 1, size=(3, 3)); S = (S + S.T) / 2
Q = S * pi[:, None]                                # Q[b,a] = rate a->b = S[b,a] pi_b (detailed balance)
np.fill_diagonal(Q, 0); Q -= np.diag(Q.sum(axis=0))
P = np.eye(3) + 0.2 * Q                            # P[b,a] = Prob(a -> b), columns sum to 1
cyc = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]]) * 1.0
Qc = cyc - np.eye(3)
Pc = np.eye(3) + 0.2 * Qc
pic = np.ones(3) / 3
# stationary weights at the ends: start with pi(x_0); in the symmetric (sqrt pi) frame
D = np.diag(np.sqrt(pi)); Dc = np.diag(np.sqrt(pic))
vrev = rp_violation(gram(np.linalg.inv(D) @ P @ D, 4, u=np.sqrt(pi)))
vcyc = rp_violation(gram(np.linalg.inv(Dc) @ Pc @ Dc, 4, u=np.sqrt(pic)))
check("R10 reversible Markov chain RP; one-way cycle not", vrev < TOL and vcyc > 1e-3,
      f"reversible {vrev:.1e}; cycle {vcyc:.2e}")

# R11 the turn: continuing the Euclidean semigroup e^{-tH} to t = i s gives generators
#     A(H) = J H (hbar = 1) over the reals, with ONE J for every H, so (A1)-(A4) hold
n = 3
J = J_of(n)
Hs = [rand_sym(n, rng.normal(size=n)) for _ in range(5)] + [np.eye(n)]
ok = True
worst = 0.0
for Hh in Hs:
    s = 0.37
    U = complexify(expm(-1j * s * Hh))             # continued evolution, written over R^{2n}
    A_ = complexify(-1j * Hh)                      # its generator, over the reals
    HR = complexify(Hh + 0j)
    errs = [np.abs(U - expm(s * A_)).max(),
            np.abs(A_ - (-J) @ HR).max(),          # A = -J H  (i <-> -J as in the paper)
            np.abs(A_ + A_.T).max(),               # (A1)
            np.abs(A_ @ HR - HR @ A_).max(),       # (A2)
            np.abs(A_ @ A_ + HR @ HR).max()]       # (A3)
    for c in (-1.3, 0.4, 2.0):                     # (A4)
        Ac = complexify(-1j * (Hh + c * np.eye(n)))
        errs.append(np.abs(Ac - (A_ + c * complexify(-1j * np.eye(n)))).max())
    worst = max(worst, max(errs))
check("R11 Wick-continued RP family satisfies (A1)-(A4) with one shared J", worst < 1e-9,
      f"worst error {worst:.1e} over {len(Hs)} Hamiltonians")

# R12 the i is forced at the turn: e^{-isH} is a real matrix for all s only when H = 0
s_vals = np.linspace(0.1, 2, 7)
imag_parts = [max(np.abs(expm(-1j * s * Hh).imag).max() for s in s_vals) for Hh in Hs]
check("R12 continued evolution is never real for H != 0", min(imag_parts) > 1e-2,
      f"smallest max|Im| = {min(imag_parts):.3f}")

# R13 the energy-zero shift on the Euclidean side is one positive number for every H
c = 0.8
ratios = [expm(-0.5 * (Hh + c * np.eye(n))) @ np.linalg.inv(expm(-0.5 * Hh)) for Hh in Hs]
check("R13 Euclidean energy shift = common factor e^(-ct) for every H",
      max(np.abs(R - np.exp(-0.5 * c) * np.eye(n)).max() for R in ratios) < 1e-12,
      f"factor {np.exp(-0.5 * c):.6f}")

# R14 real QM's per-dynamics complex structures (polar J_H) differ between Hamiltonians,
#     so real QM is not the Wick image of a single RP family (paper's Remark example)
J0 = np.kron(np.eye(2), np.array([[0, -1], [1, 0]]))
Pm = np.eye(4)[[0, 2, 1, 3]]
J2 = Pm @ J0 @ Pm.T
check("R14 two real-QM dynamics carry different J's", np.abs(J0 - J2).max() > 0.5,
      f"|J1 - J2| = {np.abs(J0 - J2).max():.0f}")

# R15 the dial picture: real numbers a with 1 + k a > 0 form the identity's component of
#     the weighted-circle group: closed under o_k and under its inverse; the wall -1/k absorbs
k = -2.0
pts = rng.uniform(-3, 0.5, size=400)
pos = pts[1 + k * pts > 0]
comp = lambda a, b: a + b + k * a * b              # noqa: E731
inv = lambda a: -a / (1 + k * a)                    # noqa: E731
closed = all(1 + k * comp(a, b) > 0 for a in pos for b in pos[:40])
invok = all(1 + k * inv(a) > 0 and abs(comp(a, inv(a))) < 1e-9 for a in pos)
wall = -1 / k
absorbs = all(abs(comp(wall, b) - wall) < 1e-12 for b in pts)
check("R15 positive-multiplier side of the wall is a subgroup; the wall absorbs",
      closed and invok and absorbs, f"k = {k}, wall at a = {wall}")

print(f"\n{sum(results)}/{len(results)} checks pass")
