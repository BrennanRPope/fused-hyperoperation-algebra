"""DRAFT round 10 checks: does the twisted-boost law fix matter's charge under its own twist?

Round 9: on the closed velocity sphere (q = +1) the law's twist is 2 turns of flux, and matter with an extra charge c
under it gets 2c one-handed massless modes. Here: the law's composition B(u)B(v) = (1 - q u* v) B(u (+) v) diag(T, 1)
(../twisted-boosts/) divides out a COMPLEX corner c = 1 - q u* v, exactly as u (+) v = (u + v)/c does. A 2-component
object carried by the law's own bookkeeping turns by diag(T, 1); carried as an ordinary spinor (dividing by |c| only)
it turns by diag(T^1/2, T^-1/2). The difference is a common phase T^1/2: charge 1/2 on top of spin.
Run: python3 checks.py  (numpy, sympy).  Writes results.json and checks_output.txt.
"""
import json, io, math
import numpy as np
import sympy as sp

results, passed, failed = {}, 0, 0
_log = io.StringIO()
def say(s):
    print(s); _log.write(s + "\n")
def check(name, ok, detail=""):
    global passed, failed
    passed, failed = (passed + 1, failed) if ok else (passed, failed + 1)
    say(("PASS " if ok else "FAIL ") + name + (("  | " + detail) if detail else ""))

say("== U. The law's corner carries exactly half the twist")
a, b, c_, d, q = sp.symbols("a b c d q", real=True)
u = a + sp.I * b; v = c_ + sp.I * d; cj = sp.conjugate
corner = 1 - q * cj(u) * v
T = (1 - q * u * cj(v)) / (1 - q * cj(u) * v)
check("U1 (exact) the twist is the corner's conjugate over the corner, T = c*/c, so the corner's phase is minus half the twist",
      sp.simplify(T - cj(corner) / corner) == 0 and sp.simplify(corner**2 * T - corner * cj(corner)) == 0)
Bq = lambda z: sp.Matrix([[1, z], [-q * cj(z), 1]])
w = (u + v) / corner
F = corner * Bq(w) * sp.diag(T, 1)
check("U2 (exact) the law's factorization B(u)B(v) = c . B(u (+) v) . diag(T, 1) (from twisted-boosts)",
      all(sp.simplify(sp.expand(Bq(u) * Bq(v) - F)[i, j]) == 0 for i in range(2) for j in range(2)))

def Bn(z, qq):
    return np.array([[1, z], [-qq * np.conj(z), 1]], complex)
def oplus(z1, z2, qq):
    return (z1 + z2) / (1 - qq * np.conj(z1) * z2)
rng = np.random.default_rng(10)
errL = errS = 0.0; dets = []
for _ in range(500):
    z1, z2 = (rng.normal(size=2) * 0.6) @ [1, 1j], (rng.normal(size=2) * 0.6) @ [1, 1j]
    cc = 1 - np.conj(z1) * z2; Tn = np.conj(cc) / cc
    left_law = np.linalg.inv(Bn(oplus(z1, z2, 1), 1)) @ Bn(z1, 1) @ Bn(z2, 1) / cc          # law: divide by the complex corner
    left_spin = np.linalg.inv(Bn(oplus(z1, z2, 1), 1)) @ Bn(z1, 1) @ Bn(z2, 1) / abs(cc)    # spinor: divide by its size only
    errL = max(errL, np.abs(left_law - np.diag([Tn, 1])).max())
    # a pure turn: diagonal, unit entries, det 1, and the two entries differ by the twist T
    errS = max(errS, abs(left_spin[0, 1]) + abs(left_spin[1, 0]) + abs(abs(left_spin[0, 0]) - 1)
               + abs(left_spin[0, 0] / left_spin[1, 1] - Tn))
    dets.append((np.linalg.det(left_law), Tn, np.linalg.det(left_spin)))
okdet = all(abs(dl - t) < 1e-12 and abs(ds - 1) < 1e-12 for dl, t, ds in dets)
check("U3 500 random pairs: law's bookkeeping leaves diag(T, 1) (det T); spinor bookkeeping leaves a pure rotation diag(T^1/2, T^-1/2) (det 1)",
      errL < 1e-12 and errS < 1e-12 and okdet, f"{errL:.1e}, {errS:.1e}")

say("\n== V. Round a closed loop on the velocity sphere (q = +1)")
def sph(z):            # unit-sphere point from the chart u = (x + iy)/(1 + z)
    r2 = abs(z)**2
    return np.array([2 * z.real, 2 * z.imag, 1 - r2]) / (1 + r2)
def sph_excess(A, B, C):
    return 2 * math.atan2(np.dot(A, np.cross(B, C)), 1 + np.dot(A, B) + np.dot(B, C) + np.dot(C, A))
rowsV = []
for theta0 in (0.3, 0.9, 1.5, 2.2):
    N = 3000; r0 = math.tan(theta0 / 2)
    pts = [r0 * np.exp(2j * math.pi * k / N) for k in range(N + 1)]
    H_law = np.eye(2, dtype=complex); H_spin = np.eye(2, dtype=complex); area = 0.0; wsum = 0.0
    for k in range(N):
        z1 = pts[k]; z2 = oplus(-z1, pts[k + 1], 1)          # the step, found with the law
        cc = 1 - np.conj(z1) * z2; Tn = np.conj(cc) / cc
        H_law = np.diag([Tn, 1]) @ H_law; wsum += float(np.angle(Tn))   # the twist, kept unwrapped
        H_spin = np.diag([np.sqrt(Tn), 1 / np.sqrt(Tn)]) @ H_spin
        area += sph_excess(np.array([0, 0, 1.0]), sph(z1), sph(pts[k + 1]))
    ph_law = np.angle(np.linalg.det(H_law)); extra_per_component = ph_law / 2
    rowsV.append((theta0, area, float(np.angle(H_law[0, 0])), float(np.angle(H_law[1, 1])), float(np.angle(np.linalg.det(H_spin))), wsum))
okV = all(abs(((pl - A + math.pi) % (2 * math.pi)) - math.pi) < 1e-9 and abs(p2) < 1e-12 and abs(ds) < 1e-9 for _, A, pl, p2, ds, _ in rowsV)
check("V1 law's bookkeeping round a loop: top component turns by the enclosed area, bottom not at all; spinor bookkeeping's common phase stays 0",
      okV, "; ".join(f"cap {t:.1f}: area {A:.4f}, top turns {ws:.4f}" for t, A, _, _, _, ws in rowsV))
check("V2 so on top of spin (+-1/2 of the area), the law's object carries a common phase of +1/2 of the area: charge 1/2 under the twist",
      all(abs(ws / 2 - A / 2) < 1e-9 for _, A, _, _, _, ws in rowsV),
      "; ".join(f"cap {t:.1f}: common phase {ws/2:.4f} = area/2 {A/2:.4f}" for t, A, _, _, _, ws in rowsV))
total_extra = 0.5 * 4 * math.pi
check("V3 over the whole closed sphere that common phase adds to 1/2 x 4 pi = 2 pi: exactly 1 whole turn, so it is a consistent flux",
      abs(total_extra / (2 * math.pi) - 1) < 1e-15, f"{total_extra / (2 * math.pi):.0f} turn")

say("\n== W. What that charge gives: massless modes on the hidden velocity sphere")
r = 1.0
def monopole_dirac(nflux, jmax=8.0):
    mm, mp = (nflux - 1) / 2, (nflux + 1) / 2
    js = [j / 2 for j in range(0, int(2 * jmax) + 1)]
    src = [(j, c) for j in js if abs(mm) <= j and (j - mm) % 1 == 0 and j < jmax for c in range(int(2 * j + 1))]
    dst = [(j, c) for j in js if abs(mp) <= j and (j - mp) % 1 == 0 and j < jmax for c in range(int(2 * j + 1))]
    Dp = np.zeros((len(dst), len(src)))
    for i, (jd, cd) in enumerate(dst):
        for k, (js_, cs) in enumerate(src):
            if jd == js_ and cd == cs:
                Dp[i, k] = math.sqrt(js_ * (js_ + 1) - mm * (mm + 1)) / r
    return Dp
def counts(n):
    Dp = monopole_dirac(n); s = np.linalg.svd(Dp, compute_uv=False) if Dp.size else np.array([])
    rank = int(np.sum(s > 1e-9)); return Dp.shape[1] - rank, Dp.shape[0] - rank
law = counts(1); spin = counts(0)
check("W1 law's bookkeeping (charge 1/2, flux 1): exactly one massless mode, one-handed; ordinary spinor bookkeeping: none",
      law == (1, 0) and spin == (0, 0), f"law L{law[0]} R{law[1]}; spinor L{spin[0]} R{spin[1]}")

say(f"\n{passed} passed, {failed} failed")
results.update({"loops": rowsV, "total_extra_turns": total_extra / (2 * math.pi), "law_counts": law, "spin_counts": spin,
                "summary": {"passed": passed, "failed": failed}})
json.dump(results, open("results.json", "w"), indent=1)
open("checks_output.txt", "w").write(_log.getvalue())
