"""DRAFT round 36: where could hypercharge and colour come from? Round 33 found that the law's hidden 3-sphere gives the weak
SU(2) and a circle, but that hypercharge is not the circle weight and colour is missing. This round asks what symmetry the law
itself has, and in how many dimensions it can live. (C1) the law needs no algebra: on unit vectors of any space it is the
rotation in the plane of 1 and x that carries 1 to x, applied to y; this equals the sandwich x^(1/2) y x^(1/2) for complex
numbers, quaternions and octonions, and still lives on the sphere where the sedenion sandwich falls off it. (C2) so its
symmetries, the rotations that keep 1 fixed and respect the law, are all rotations of the other directions, SO(n-1), mirrors
included; octonions add nothing special. (C3) which hidden sphere holds colour, the weak force and hypercharge together with
one-handed matter: the first is the law on the 10-sphere, whose symmetry SO(10) is the textbook unified group; its spinor is
exactly one generation with the measured hypercharges, every anomaly cancelled. (C4) the weak mixing angle it fixes, 3/8, run
down to the Z mass with the Standard Model's couplings. C1 and C2 are numbers; C3 exact; C4 exact plus textbook one-loop running."""
import json, io, time, math, itertools
from fractions import Fraction as Fr
import numpy as np
from scipy.linalg import expm, null_space

T0 = time.time(); _log = io.StringIO(); passed = failed = 0; results = {}
def say(s): print(s); _log.write(s + "\n")
def check(name, ok, detail):
    global passed, failed
    passed += bool(ok); failed += (not ok)
    say(f"{'PASS' if ok else 'FAIL'} {name}  | {detail}")
rng = np.random.default_rng(0)
def runit(n): v = rng.normal(size=n); return v / np.linalg.norm(v)

# ---------------------------------------------------------------- C1
say("== C1. The law needs no algebra: it is a rotation in the plane of 1 and x")
def cd_mul(a, b):              # Cayley-Dickson product: reals -> complex -> quaternions -> octonions -> sedenions
    n = len(a)
    if n == 1: return a * b
    h = n // 2; a1, a2, b1, b2 = a[:h], a[h:], b[:h], b[h:]
    cj = lambda x: np.concatenate([x[:1], -x[1:]])
    return np.concatenate([cd_mul(a1, b1) - cd_mul(cj(b2), a2), cd_mul(b2, a1) + cd_mul(a2, cj(b1))])
def sqrt_unit(x):
    th = np.arccos(np.clip(x[0], -1, 1)); v = x.copy(); v[0] = 0; nv = np.linalg.norm(v)
    out = np.zeros_like(x); out[0] = np.cos(th / 2)
    if nv > 1e-15: out += np.sin(th / 2) * v / nv
    return out
def law_sandwich(x, y):
    s = sqrt_unit(x); return cd_mul(cd_mul(s, y), s)
def law(x, y):                 # the rotation in the plane of e0 = 1 and x that carries 1 to x, applied to y
    n = len(x); e0 = np.zeros(n); e0[0] = 1
    th = np.arccos(np.clip(x[0], -1, 1)); v = x.copy(); v[0] = 0; nv = np.linalg.norm(v)
    u = v / nv if nv > 1e-15 else np.zeros(n)
    a, b = y @ e0, y @ u
    return y + (a * (np.cos(th) - 1) - b * np.sin(th)) * e0 + (a * np.sin(th) + b * (np.cos(th) - 1)) * u
C1 = {}
for n, name in ((2, "complex"), (4, "quaternions"), (8, "octonions"), (16, "sedenions")):
    pairs = [(runit(n), runit(n)) for _ in range(400)]
    C1[name] = dict(n=n, rotation_vs_sandwich=float(max(np.linalg.norm(law(x, y) - law_sandwich(x, y)) for x, y in pairs)),
                    sandwich_off_sphere=float(max(abs(np.linalg.norm(law_sandwich(x, y)) - 1) for x, y in pairs)),
                    rotation_off_sphere=float(max(abs(np.linalg.norm(law(x, y)) - 1) for x, y in pairs)),
                    carries_1_to_x=float(max(np.linalg.norm(law(x, np.eye(n)[0]) - x) for x, y in pairs)))
ok1 = all(C1[k]["rotation_vs_sandwich"] < 1e-12 for k in ("complex", "quaternions", "octonions")) and C1["sedenions"]["sandwich_off_sphere"] > 0.1 \
      and all(C1[k]["rotation_off_sphere"] < 1e-12 and C1[k]["carries_1_to_x"] < 1e-12 for k in C1)
check("C1 (numbers; textbook: the sphere as a symmetric space, its transvections) the law x (+) y = x^(1/2) y x^(1/2) is the rotation in the plane of 1 and x that carries 1 to x, applied to y, leaving every direction across that plane alone. That holds for complex numbers, quaternions and octonions. With sedenions the sandwich falls off the sphere (they don't keep lengths), but the rotation is defined on every sphere. So the law needs no multiplication at all: only the round sphere and a chosen point 1",
      ok1, "; ".join(f"{k} (n={v['n']}): rotation vs sandwich {v['rotation_vs_sandwich']:.1e}, sandwich off sphere {v['sandwich_off_sphere']:.2f}, rotation off sphere {v['rotation_off_sphere']:.1e}, carries 1 to x {v['carries_1_to_x']:.1e}" for k, v in C1.items()))
results["C1"] = C1

# ---------------------------------------------------------------- C2
say("== C2. The law's symmetries: every rotation of the other directions, mirrors included")
def sym_dim(n, h=1e-5):
    basis = []
    for i in range(n):
        for j in range(i + 1, n):
            E = np.zeros((n, n)); E[i, j] = 1; E[j, i] = -1; basis.append(E)
    rows = []
    for _ in range(n + 4):
        x, y = runit(n), runit(n); z = law(x, y)
        rows.append(np.array([(law(expm(h * E) @ x, expm(h * E) @ y) - law(expm(-h * E) @ x, expm(-h * E) @ y)) / (2 * h) - E @ z for E in basis]).T)
    N = null_space(np.vstack(rows), rcond=1e-6)
    fixes1 = all(np.linalg.norm(sum(c * E for c, E in zip(N[:, k], basis))[:, 0]) < 1e-6 for k in range(N.shape[1]))
    return N.shape[1], fixes1
C2 = []
for n in (3, 4, 5, 6, 8, 11):
    d, f1 = sym_dim(n)
    Rf = np.eye(n); Rf[n - 1, n - 1] = -1                           # a mirror that keeps 1
    mirror = max(np.linalg.norm(law(Rf @ x, Rf @ y) - Rf @ law(x, y)) for x, y in [(runit(n), runit(n)) for _ in range(50)])
    C2.append(dict(sphere=n - 1, dim=d, so_n_minus_1=(n - 1) * (n - 2) // 2, fixes_1=f1, mirror=float(mirror)))
# octonion product automorphisms (G2) are only part of it
def g2_dim(h=1e-5):
    n = 8; basis = []
    for i in range(n):
        for j in range(i + 1, n):
            E = np.zeros((n, n)); E[i, j] = 1; E[j, i] = -1; basis.append(E)
    rows = []
    for _ in range(12):
        x, y = runit(n), runit(n); z = cd_mul(x, y)
        rows.append(np.array([E @ z - cd_mul(E @ x, y) - cd_mul(x, E @ y) for E in basis]).T)
    return null_space(np.vstack(rows), rcond=1e-8).shape[1]
g2 = g2_dim()
ok2 = all(r["dim"] == r["so_n_minus_1"] and r["fixes_1"] and r["mirror"] < 1e-12 for r in C2) and g2 == 14
check("C2 (numbers) the rotations that respect the law are exactly those that keep 1 fixed: all of SO(n-1) on the sphere S^(n-1), with mirrors too. On round 1's 3-sphere that is SO(3). On the octonions' 7-sphere it is SO(7), 21 dimensions, not just the 14 of G2 that respect octonion multiplication: the law doesn't see the octonions. And because mirrors respect it, the law on any sphere still picks no hand (as round 32 found on the 3-sphere)",
      ok2, "; ".join(f"S^{r['sphere']}: {r['dim']} (SO({r['sphere']}) has {r['so_n_minus_1']}), all keep 1: {r['fixes_1']}, mirror {r['mirror']:.0e}" for r in C2) + f"; octonion product's own symmetries: {g2} (G2)")
results["C2"] = dict(rows=C2, g2=g2)

# ---------------------------------------------------------------- C3
say("== C3. Which sphere holds colour, the weak force and hypercharge, with one-handed matter (exact)")
# (a) the size needed: SO(m) must hold SU(3) x SU(2) x U(1) and have complex (one-handed) spinors.
#     SU(3)'s smallest real representation is 6 (C^3); a commuting SU(2) needs at least 3 more directions; SO(m) spinors are
#     complex, with the two hands each other's conjugates, only for m = 2 mod 4 (textbook). So m >= 9 and m = 2 mod 4: m = 10.
cands = [m for m in range(2, 15) if m >= 9 and m % 4 == 2]
# (b) SO(10) = rotations of 10 directions, split 6 + 4: SO(6) = SU(4) (colour plus lepton number), SO(4) = SU(2)_L x SU(2)_R.
#     Its spinor of one hand: 16 weights (+-1/2)^5 with an even number of minus signs.
W = [w for w in itertools.product((Fr(1, 2), Fr(-1, 2)), repeat=5) if sum(1 for a in w if a < 0) % 2 == 0]
def charges(w):
    T3L = (w[3] - w[4]) / 2; T3R = (w[3] + w[4]) / 2; BL = Fr(2, 3) * (w[0] + w[1] + w[2])
    Y = T3R + BL / 2; return T3L, T3R, BL, Y, T3L + Y
rows = []
for w in W:
    T3L, T3R, BL, Y, Q = charges(w); nminus = sum(1 for a in w[:3] if a < 0)
    colour = "singlet" if nminus in (0, 3) else ("triplet" if nminus == 1 else "antitriplet")
    weak = "doublet" if T3L != 0 else "singlet"
    rows.append(dict(w=[str(a) for a in w], colour=colour, weak=weak, Y=Y, Q=Q, BL=BL))
def name(r):
    key = (r["colour"], r["weak"], r["Y"])
    return {("triplet", "doublet", Fr(1, 6)): "Q", ("antitriplet", "singlet", Fr(-2, 3)): "u^c", ("antitriplet", "singlet", Fr(1, 3)): "d^c",
            ("singlet", "doublet", Fr(-1, 2)): "L", ("singlet", "singlet", Fr(1)): "e^c", ("singlet", "singlet", Fr(0)): "nu^c"}.get(key, "?")
counts = {}
for r in rows: counts[name(r)] = counts.get(name(r), 0) + 1
sm_counts = {"Q": 6, "u^c": 3, "d^c": 3, "L": 2, "e^c": 1, "nu^c": 1}
TrY = sum(r["Y"] for r in rows); TrY3 = sum(r["Y"]**3 for r in rows)
TrT2Y = sum((charges([Fr(a) for a in r["w"]])[0])**2 * r["Y"] for r in rows)
# colour: T8-like generator of SU(3) inside SO(6), on triplets; Tr T_c^2 Y over coloured states, with T_c = (w0 - w1)/2
TrTc2Y = sum(((Fr(r["w"][0]) - Fr(r["w"][1])) / 2)**2 * r["Y"] for r in rows if r["colour"] != "singlet")
TrGrav = TrY
doublets = sum(1 for r in rows if r["weak"] == "doublet") // 2
ok3 = cands[0] == 10 and counts == sm_counts and TrY == 0 and TrY3 == 0 and TrT2Y == 0 and TrTc2Y == 0 and doublets % 2 == 0
check("C3 (exact; textbook SO(10) unification, Georgi 1974 and Fritzsch-Minkowski 1975) so the law on a 10-sphere has the symmetry SO(10). For colour (SU(3) needs 6 directions), a commuting weak SU(2) (3 more) and one-handed matter (spinors of SO(m) come in two distinct hands only for m = 2 mod 4), 10 is the smallest. Split its 10 directions 6 + 4: the 6 carry colour and lepton number, the 4 carry the weak force, left and right. One hand of its spinor is exactly one generation: Q, u^c, d^c, L, e^c and a right-handed neutrino, with the measured hypercharges, Y = T3R + (B - L)/2, fixed by the group with nothing to choose. Every anomaly cancels and the doublets come in an even number. This fills both of round 33's gaps at once: hypercharge is not the circle weight but a rotation of the 10 directions, and colour is the rotation of 6 of them",
      ok3, f"sizes that work up to 14: {cands}; the 16 states: {counts} (Standard Model generation plus nu^c: {sm_counts}); Tr Y = {TrY}, Tr Y^3 = {TrY3}, Tr T3L^2 Y = {TrT2Y}, Tr T_colour^2 Y = {TrTc2Y}; weak doublets {doublets} (even: no Witten anomaly); electric charges {sorted(set(str(r['Q']) for r in rows))}")
results["C3"] = dict(candidates=cands, states=[dict(w=r["w"], colour=r["colour"], weak=r["weak"], Y=str(r["Y"]), Q=str(r["Q"]), BL=str(r["BL"]), name=name(r)) for r in rows],
                     counts=counts, traces=dict(Y=str(TrY), Y3=str(TrY3), T2Y=str(TrT2Y), Tc2Y=str(TrTc2Y)), doublets=doublets)

# ---------------------------------------------------------------- C4
say("== C4. The weak mixing angle SO(10) fixes, run down to the Z mass")
TrT3 = sum(charges([Fr(a) for a in r["w"]])[0]**2 for r in rows); TrQ2 = sum(r["Q"]**2 for r in rows)
sin2_unif = TrT3 / TrQ2
MZ, ainv, s2, a3 = 91.1876, 127.95, 0.2312, 0.1179           # PDG review values at the Z mass (approximate)
a1inv = ainv * (1 - s2) * 3 / 5; a2inv = ainv * s2; a3inv = 1 / a3   # alpha_1 in the unified normalisation (Y scaled by sqrt(3/5))
b = (41 / 10, -19 / 6, -7.0)
def ainv_at(mu): L = math.log(mu / MZ); return [a1inv - b[0] / (2 * math.pi) * L, a2inv - b[1] / (2 * math.pi) * L, a3inv - b[2] / (2 * math.pi) * L]
def meet(i, j): return MZ * math.exp(2 * math.pi * ([a1inv, a2inv, a3inv][i] - [a1inv, a2inv, a3inv][j]) / (b[i] - b[j]))
m12, m23, m13 = meet(0, 1), meet(1, 2), meet(0, 2)
# prediction: impose one meeting point for all three (from alpha_em and alpha_3), solve for sin^2 at the Z mass
from scipy.optimize import brentq
def mismatch(s):
    A1, A2 = ainv * (1 - s) * 3 / 5, ainv * s
    M = MZ * math.exp(2 * math.pi * (A2 - a3inv) / (b[1] - b[2]))
    return A1 - b[0] / (2 * math.pi) * math.log(M / MZ) - (A2 - b[1] / (2 * math.pi) * math.log(M / MZ))
s2_pred = brentq(mismatch, 0.15, 0.3)
grid = np.logspace(math.log10(MZ), 19, 200); run = np.array([ainv_at(m) for m in grid])
ok4 = sin2_unif == Fr(3, 8) and m12 < m13 < m23 and m23 / m12 > 100 and abs(s2_pred - s2) > 0.01
check("C4 (exact 3/8; textbook one-loop running) SO(10) fixes the weak mixing angle where it holds: sin^2 theta_W = Tr T3^2 / Tr Q^2 = 3/8 over one generation. Run down with the Standard Model's own particles, the three couplings never meet at one point: hypercharge meets the weak force near 10^13 GeV, hypercharge meets colour near 2 x 10^14 and the weak force meets colour near 10^17. Forcing a single meeting predicts sin^2 theta_W = 0.21 at the Z mass, where 0.231 is measured. So SO(10) broken straight to the Standard Model misses; it needs something in between (in textbook models an intermediate stage such as SU(4) x SU(2) x SU(2), or new particles), and the law doesn't supply that yet. Compare round 33's reading of the 3-sphere, which gave 1/2",
      ok4, f"sin^2 at unification = {sin2_unif}; one-loop meeting scales (GeV): hypercharge-weak {m12:.2e}, hypercharge-colour {m13:.2e}, weak-colour {m23:.2e}; single-meeting prediction sin^2(MZ) = {s2_pred:.4f} against measured {s2}")
results["C4"] = dict(sin2_unif=str(sin2_unif), meet12=m12, meet13=m13, meet23=m23, s2_pred=s2_pred, s2_meas=s2, grid=grid.tolist(), run=run.tolist())

say(f"\n{passed} / {passed + failed} pass  ({time.time() - T0:.1f} s)")
results["summary"] = dict(passed=passed, failed=failed)
json.dump(results, open("results.json", "w"), indent=1, default=str)
open("checks_output.txt", "w").write(_log.getvalue())
