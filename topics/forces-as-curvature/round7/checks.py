"""DRAFT round 7 checks: can twisting the hidden space make left-handed matter?

Round 5 argued from textbook index theory that a twist (torsion) cannot make the net handedness of massless
matter nonzero, and that magnetic flux through the hidden space can. Here both are shown on explicit spectra.
The Dirac operator on each hidden space is built as a matrix from spin-j angular momentum matrices
(Peter-Weyl: functions on SU(2) = sum over j of (2j+1) copies of V_j) and diagonalised numerically.
These are illustrations of textbook theorems (Atiyah-Singer, Lichnerowicz), not independent tests of them.
Run: python3 checks.py  (numpy).  Writes results.json and checks_output.txt.
"""
import json, io, math
import numpy as np

results, passed, failed = {}, 0, 0
_log = io.StringIO()
def say(s):
    print(s); _log.write(s + "\n")
def check(name, ok, detail=""):
    global passed, failed
    passed, failed = (passed + 1, failed) if ok else (passed, failed + 1)
    say(("PASS " if ok else "FAIL ") + name + (("  | " + detail) if detail else ""))

def Jmats(j):
    m = np.arange(j, -j - 1, -1); d = len(m)
    Jp = np.zeros((d, d))
    for i in range(1, d):
        Jp[i - 1, i] = math.sqrt(j * (j + 1) - m[i] * (m[i] + 1))
    Jz = np.diag(m).astype(complex); Jx = (Jp + Jp.T) / 2; Jy = (Jp - Jp.T) / 2j
    return [Jx.astype(complex), Jy, Jz], Jp, m
sig = [np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]]), np.array([[1, 0], [0, -1]], complex)]
r = 1.0; JMAX = 4.0

# ------------------------------------------------------------------ N. the 3-sphere with a twist of strength t
say("== N. Dirac operator on the hidden 3-sphere, untwisted (t = 0) and twisted toward the flat, parallelizing twist (t = 1)")
def D3_block(j, t):
    """Dirac operator on S^3 restricted to one copy of V_j (x C^2 spinor): left-invariant frame e_a = -(2i/r) J_a,
    spin connection of the connection with twist t: (1 - t) times the Levi-Civita one (t = 1: Cartan-Schouten, flat)."""
    J, _, _ = Jmats(j)
    return (2 / r) * sum(np.kron(J[a], sig[a]) for a in range(3)) + 1.5 * (1 - t) / r * np.eye(2 * len(J[0]))
def spectrum_S3(t):
    ev = []
    for j2 in range(int(2 * JMAX) + 1):
        j = j2 / 2; w = np.linalg.eigvalsh(D3_block(j, t))
        ev += [(x, int(2 * j + 1)) for x in w]          # each eigenvalue comes with 2j+1 copies
    return ev
ev0 = spectrum_S3(0.0)
vals = {}
for x, mlt in ev0:
    key = round(x, 9); vals[key] = vals.get(key, 0) + mlt
ok, worst = True, 0.0
for n in range(0, 7):
    for sgn in (+1, -1):
        lam = sgn * (n + 1.5) / r; got = vals.get(round(lam, 9), 0)
        ok &= got == (n + 1) * (n + 2)
check("N1 untwisted 3-sphere: eigenvalues +-(n + 3/2)/r with multiplicity (n+1)(n+2) (Baer 1996), built from spin matrices",
      ok and min(abs(x) for x, _ in ev0) > 1.5 / r - 1e-9, f"smallest |lambda| = {min(abs(x) for x, _ in ev0):.6f} (Lichnerowicz bound 3/(2r))")
ev1 = spectrum_S3(1.0)
z1 = sum(m for x, m in ev1 if abs(x) < 1e-9)
check("N2 full parallelizing twist (t = 1, the sphere's curvature twisted away): exactly 2 massless modes appear, the constant spinors",
      z1 == 2, f"zero modes {z1}; untwisted: {sum(m for x, m in ev0 if abs(x) < 1e-9)}")
tgrid = np.linspace(-1.5, 3.5, 1001)
flow = []
for tt in tgrid:
    flow.append(sorted(x for x, _ in spectrum_S3(tt) if abs(x) < 6))
crossings = []
for j2 in range(0, 5):
    j = j2 / 2
    crossings += [1 + 4 * j / 3] + ([1 - 4 * (j + 1) / 3] if j > 0 else [])
zero_at = sorted(c for c in crossings if -1.5 <= c <= 3.5)
okc = all(sum(m for x, m in spectrum_S3(c) if abs(x) < 1e-9) > 0 for c in zero_at)
gap_away = min(min(abs(x) for x in f) for tt, f in zip(tgrid, flow) if min(abs(tt - c) for c in zero_at) > 0.01)
okc = okc and gap_away > 1e-3
check("N3 zero modes appear only at special twist strengths, t = 1 + 4j/3 and 1 - 4(j+1)/3", okc,
      "t = " + ", ".join(f"{c:.3f}" for c in zero_at) + f"; elsewhere smallest |lambda| >= {gap_away:.4f}")

# ------------------------------------------------------------------ O. handedness: hidden S^3 x S^1 (even-dimensional)
say("\n== O. Handedness needs an even number of hidden directions: hidden S^3 x S^1, chirality = 1 x tau_3")
tau = sig
def chiral_counts_S3S1(t, nmax=3, Rc=0.7):
    """D = D3 x tau_1 + (n/R) 1 x tau_2 on each (j, n) sector; it flips handedness (anticommutes with 1 x tau_3).
    Count zero modes of each handedness."""
    L = Rr = 0
    for j2 in range(int(2 * JMAX) + 1):
        j = j2 / 2
        for n in range(-nmax, nmax + 1):
            D = np.kron(D3_block(j, t), tau[0]) + (n / Rc) * np.kron(np.eye(2 * int(2 * j + 1)), tau[1])
            G5 = np.kron(np.eye(2 * int(2 * j + 1)), tau[2])
            assert np.abs(D @ G5 + G5 @ D).max() < 1e-12
            w, v = np.linalg.eigh(D)
            ker = v[:, np.abs(w) < 1e-9]
            if ker.shape[1]:
                ch = np.linalg.eigvalsh(ker.conj().T @ G5 @ ker)
                L += int(2 * j + 1) * int(np.sum(ch > 0.5)); Rr += int(2 * j + 1) * int(np.sum(ch < -0.5))
    return L, Rr
rowsO = []
for t in (0.0, 0.5, 1.0, 1 + 4 * 0.5 / 3, 1 - 4 * 1.5 / 3, 2.0):
    rowsO.append((t,) + chiral_counts_S3S1(t))
check("O1 hidden S^3 x S^1 at every twist strength: massless modes, when the twist makes any, come as many left- as right-handed",
      all(a == b for _, a, b in rowsO) and any(a > 0 for _, a, b in rowsO),
      "; ".join(f"t={t:.2f}: L{a} R{b}" for t, a, b in rowsO))

# ------------------------------------------------------------------ P. the way out: flux through a hidden 2-sphere
say("\n== P. The textbook way out: magnetic flux through a hidden 2-sphere (Randjbar-Daemi, Salam, Strathdee 1983)")
def monopole_dirac(nflux, jmax=8.0):
    """Dirac operator on S^2 = SU(2)/U(1) coupled to flux nflux. Spinor halves are sections of weight
    m- = (n-1)/2 and m+ = (n+1)/2; D+ acts as J_+/r inside each V_j (2j+1 copies). Returns D+ as a matrix."""
    mm, mp = (nflux - 1) / 2, (nflux + 1) / 2
    rows, cols, blocks = [], [], []
    js = [j / 2 for j in range(0, int(2 * jmax) + 1)]
    src = [(j, c) for j in js if abs(mm) <= j and (j - mm) % 1 == 0 for c in range(int(2 * j + 1))]
    dst = [(j, c) for j in js if abs(mp) <= j and (j - mp) % 1 == 0 for c in range(int(2 * j + 1))]
    Dp = np.zeros((len(dst), len(src)))
    for i, (jd, cd) in enumerate(dst):
        for k, (js_, cs) in enumerate(src):
            if jd == js_ and cd == cs:
                Dp[i, k] = math.sqrt(js_ * (js_ + 1) - mm * (mm + 1)) / r
    # top shell: keep only j < jmax on both sides so truncation does not fake modes
    keep_s = [k for k, (j, c) in enumerate(src) if j < jmax]; keep_d = [i for i, (j, c) in enumerate(dst) if j < jmax]
    return Dp[np.ix_(keep_d, keep_s)]
def kernels(Dp, tol=1e-9):
    s = np.linalg.svd(Dp, compute_uv=False) if Dp.size else np.array([])
    rank = int(np.sum(s > tol))
    return Dp.shape[1] - rank, Dp.shape[0] - rank       # ker D+ (one handedness), ker D- = ker D+^dagger (other)
rowsP = []
for nf in range(0, 5):
    kp, km = kernels(monopole_dirac(nf)); rowsP.append((nf, kp, km))
Dp0 = monopole_dirac(0); s0 = np.sort(np.linalg.svd(Dp0, compute_uv=False))
okS2 = abs(s0[0] - 1 / r) < 1e-12 and abs(s0[2] - 2 / r) < 1e-12
check("P1 no flux: the 2-sphere's Dirac eigenvalues are +-(l+1)/r, no massless modes", rowsP[0][1:] == (0, 0) and okS2,
      f"lowest {s0[0]:.6f}, {s0[2]:.6f}")
check("P2 flux n through the hidden 2-sphere: n massless modes, ALL of one handedness (net handedness = n = the flux)",
      all(kp == nf and km == 0 for nf, kp, km in rowsP[1:]), "; ".join(f"n={nf}: L{kp} R{km}" for nf, kp, km in rowsP))

say(f"\n{passed} passed, {failed} failed")
results.update({"S3_untwisted_smallest": min(abs(x) for x, _ in ev0), "S3_twisted_zero_modes": z1,
                "zero_mode_twists": zero_at, "flow_t": list(tgrid[::5]), "flow": [f for f in flow[::5]],
                "S3xS1_counts": rowsO, "S2_flux_counts": rowsP, "summary": {"passed": passed, "failed": failed}})
json.dump(results, open("results.json", "w"), indent=1)
open("checks_output.txt", "w").write(_log.getvalue())
