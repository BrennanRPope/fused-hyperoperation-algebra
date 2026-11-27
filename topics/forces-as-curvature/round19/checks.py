"""DRAFT round 19 checks: is the law at q = -1 the boosts of special relativity, as at q = +1 it is the 3-sphere's group?

Round 18 found that at q = +1 the law's 3D sum is the unit-quaternion group's symmetric product, x (+) y = h.y.h, with h the
law's half kick of x. Here the same is tried at q = -1, with points placed on the velocity hyperboloid and written as 2x2
positive Hermitian matrices (the 4-velocities of special relativity), and with the Lorentz group SL(2,C) in place of SU(2).
Run: python3 checks.py  (numpy).  Writes results.json and checks_output.txt.
"""
import json, io, math, time
import numpy as np

results, passed, failed = {}, 0, 0
_log = io.StringIO()
def say(s):
    print(s, flush=True); _log.write(s + "\n")
def check(name, ok, detail=""):
    global passed, failed
    passed, failed = (passed + 1, failed) if ok else (passed, failed + 1)
    say(("PASS " if ok else "FAIL ") + name + (("  | " + detail) if detail else ""))
T_START = time.time()

# ------------------------------------------------------------------ the law (rounds 12-13, unchanged) and the two placements
def oplus(a, b, qq):
    ab = (a * b).sum(-1, keepdims=True); aa = (a * a).sum(-1, keepdims=True); bb = (b * b).sum(-1, keepdims=True)
    return ((1 - 2 * qq * ab - qq * bb) * a + (1 + qq * aa) * b) / (1 - 2 * qq * ab + qq**2 * aa * bb)
def half_kick(g, qq):
    return g / (1 + np.sqrt(1 + qq * (g * g).sum(-1, keepdims=True)))
SIG = np.array([[[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]]])
def herm(X):                      # 4-vector (X0, X1, X2, X3) -> 2x2 Hermitian X0 + X.sigma ; det = X0^2 - |X|^2
    return X[..., 0, None, None] * np.eye(2) + np.einsum('...i,ijk->...jk', X[..., 1:], SIG)
def unherm(H):
    return np.stack([np.real(H[..., 0, 0] + H[..., 1, 1]) / 2, np.real(H[..., 0, 1]), -np.imag(H[..., 0, 1]),
                     np.real(H[..., 0, 0] - H[..., 1, 1]) / 2], -1)
def to_h3(x):                     # the law's chart at q = -1 (the unit ball) -> the velocity hyperboloid X0^2 - |X|^2 = 1
    r2 = (x * x).sum(-1, keepdims=True); return np.concatenate([(1 + r2), 2 * x], -1) / (1 - r2)
Pm = lambda x: herm(to_h3(x))     # a point as a positive Hermitian matrix of determinant 1 (a 4-velocity)
hk = lambda x: half_kick(x, -1.0)
Bm = lambda x: Pm(hk(x))          # the law's half kick, written the same way
def ball(rg, n, rmax=0.9):
    v = rg.normal(size=(n, 3)); return v / np.linalg.norm(v, axis=1)[:, None] * rmax * rg.uniform(0, 1, (n, 1))**(1 / 3)
rg = np.random.default_rng(19)
X1, Y1, Z1 = ball(rg, 2000), ball(rg, 2000), ball(rg, 2000)

# ================================================================== H. The law at q = -1
say("== H. Is the law at q = -1 the boosts of special relativity?")
S = oplus(X1, Y1, -1.0); H_ = Bm(X1); PX, PY, PS = Pm(X1), Pm(Y1), Pm(S)
scale = np.abs(PS).max()
sand = np.abs(H_ @ PY @ H_ - PS).max() / scale
grp = np.abs(PX @ PY - H_ @ PS @ np.linalg.inv(H_)).max() / scale
noneq = np.abs(PX @ PY - PS).max() / scale
check("H1 (exact) the law's sum at q = -1 is the symmetric product of positive Hermitian matrices: with points on the velocity hyperboloid written as 2x2 matrices, x (+) y lands at B.Y.B, where B is the law's half kick of x. The plain product differs only by a turn (conjugation by B). This is round 18's statement with SL(2,C) in place of SU(2)",
      sand < 1e-12 and grp < 1e-12 and noneq > 0.05,
      f"2000 random pairs (relative to the largest entry): B.Y.B vs the law {sand:.0e}; X.Y vs B (x (+) y) B^-1 {grp:.0e}; X.Y vs the law, for contrast, {noneq:.2f}")
results["sandwich"] = dict(sandwich=float(sand), group_vs_law=float(grp), plain_vs_law=float(noneq))

sq = np.abs(H_ @ H_ - PX).max() / np.abs(PX).max()
herm_ok = np.abs(H_ - np.conj(np.swapaxes(H_, -1, -2))).max()
detB = np.abs(np.linalg.det(H_) - 1).max()
rest = herm(np.array([1.0, 0, 0, 0]))
boosted_rest = np.abs(H_ @ rest @ np.conj(np.swapaxes(H_, -1, -2)) - PX).max() / np.abs(PX).max()
v_law = np.linalg.norm(to_h3(X1)[:, 1:], axis=1) / to_h3(X1)[:, 0]; r = np.linalg.norm(X1, axis=1)
e_vel = np.abs(v_law - 2 * r / (1 + r * r)).max()
check("H2 (exact) the law's half kick is the pure boost to x: a positive Hermitian matrix of determinant 1 whose square is x's 4-velocity, so it carries rest to x. Then x (+) y = B.Y.B-dagger is y seen after boosting by x, i.e. relativistic velocity addition, with the chart radius r standing for the speed 2r/(1 + r^2)",
      sq < 1e-12 and herm_ok < 1e-12 and detB < 1e-12 and boosted_rest < 1e-12 and e_vel < 1e-12,
      f"B^2 vs x's 4-velocity {sq:.0e}; B Hermitian {herm_ok:.0e}, det B - 1 {detB:.0e}; B rest B-dagger vs X {boosted_rest:.0e}; speed vs 2r/(1+r^2) {e_vel:.0e}")
results["boost"] = dict(square=float(sq), hermitian=float(herm_ok), det=float(detB), rest_to_x=float(boosted_rest), speed=float(e_vel))

U = np.linalg.inv(Bm(S)) @ Bm(X1) @ Bm(Y1)
unit = np.abs(U @ np.conj(np.swapaxes(U, -1, -2)) - np.eye(2)).max(); detU = np.abs(np.linalg.det(U) - 1).max()
lhs = oplus(X1, oplus(Y1, Z1, -1.0), -1.0)
Rz = unherm(U @ herm(np.concatenate([np.zeros((len(Z1), 1)), Z1], 1)) @ np.conj(np.swapaxes(U, -1, -2)))[:, 1:]
rhs = oplus(S, Rz, -1.0); e_gyr = np.abs(lhs - rhs).max()
wig = []
for a, b in ((0.2, 0.3), (0.3, 0.5), (0.5, 0.5), (0.6, 0.7), (0.8, 0.8), (0.9, 0.2), (0.95, 0.95)):
    x = np.array([[a, 0, 0.0]]); y = np.array([[0, b, 0.0]])
    Uw = np.linalg.inv(Bm(oplus(x, y, -1.0))) @ Bm(x) @ Bm(y)
    th = 2 * math.acos(min(1.0, float(np.real(np.trace(Uw[0]))) / 2))
    g1, g2 = (1 + a * a) / (1 - a * a), (1 + b * b) / (1 - b * b)
    wig.append(dict(r1=a, r2=b, gamma1=g1, gamma2=g2, law_turn=th, wigner=math.acos((g1 + g2) / (1 + g1 * g2))))
e_wig = max(abs(w["law_turn"] - w["wigner"]) for w in wig)
check("H3 (exact) the law's twist at q = -1 is the Thomas-Wigner rotation: B(x)B(y) = B(x (+) y).U with U a rotation (unitary, det 1), and that same U is the turn by which the law fails to associate, x (+) (y (+) z) = (x (+) y) (+) (U z). For two boosts at right angles its angle is the textbook Wigner angle, cos = (g1 + g2)/(1 + g1 g2)",
      unit < 1e-12 and detU < 1e-12 and e_gyr < 1e-10 and e_wig < 1e-10,
      f"U unitary {unit:.0e}, det U - 1 {detU:.0e}; associativity up to U, 2000 triples {e_gyr:.0e}; seven right-angle pairs vs Wigner's formula {e_wig:.0e}")
results["wigner"] = dict(unitary=float(unit), det=float(detU), gyration=float(e_gyr), right_angle_pairs=wig)

# H4: the dial's sign picks the group.  The map y -> h y h (q = +1, quaternions) or Y -> B Y B (q = -1, matrices) is linear on R^4;
# find the quadratic form each preserves.
def qmul(a, b):
    a0, av = a[..., :1], a[..., 1:]; b0, bv = b[..., :1], b[..., 1:]
    return np.concatenate([a0 * b0 - (av * bv).sum(-1, keepdims=True), a0 * bv + b0 * av + np.cross(av, bv)], -1)
def to_s3(x):
    r2 = (x * x).sum(-1, keepdims=True); return np.concatenate([(1 - r2), 2 * x], -1) / (1 + r2)
E4 = np.eye(4); eu, mi = np.eye(4), np.diag([1.0, -1, -1, -1])
worst_eu, worst_mi, cross_eu, cross_mi = 0.0, 0.0, 1e9, 1e9
for x in ball(rg, 50, 0.85):
    h = to_s3(half_kick(x[None], 1.0))[0]
    Lp = np.stack([qmul(qmul(h[None], e[None]), h[None])[0] for e in E4], 1)
    B = Bm(x[None])[0]
    Lm = np.stack([unherm(B @ herm(e) @ B) for e in E4], 1)
    worst_eu = max(worst_eu, np.abs(Lp.T @ eu @ Lp - eu).max()); worst_mi = max(worst_mi, np.abs(Lm.T @ mi @ Lm - mi).max())
    cross_eu = min(cross_eu, np.abs(Lp.T @ mi @ Lp - mi).max()); cross_mi = min(cross_mi, np.abs(Lm.T @ eu @ Lm - eu).max())
qs = [-2.0, -0.5, 0.3, 1.7]
e_scale = max(np.abs(oplus(X1 * 0.5, Y1 * 0.5, qq) - oplus(math.sqrt(abs(qq)) * X1 * 0.5, math.sqrt(abs(qq)) * Y1 * 0.5,
                                                       math.copysign(1.0, qq)) / math.sqrt(abs(qq))).max() for qq in qs)
small = [1e-1, 1e-2, 1e-3, 1e-4]
dev = [float(np.abs(oplus(X1, Y1, qq) - (X1 + Y1)).max()) for qq in small]
slope = float(np.polyfit(np.log(small), np.log(dev), 1)[0])
check("H4 the dial's sign picks the group. At q = +1 the law's moves keep w^2 + x^2 + y^2 + z^2 (rotations of 4D Euclidean space, Spin(4)); at q = -1 they keep t^2 - x^2 - y^2 - z^2 (the Lorentz group, Spin(1,3)), and neither keeps the other's form. Any other q is one of these rescaled (round 16), and as q -> 0 the law becomes plain vector addition, with the difference shrinking in proportion to q",
      worst_eu < 1e-12 and worst_mi < 1e-10 and cross_eu > 0.01 and cross_mi > 0.01 and e_scale < 1e-12 and abs(slope - 1) < 0.02,
      f"50 moves each: q = +1 keeps the Euclidean form to {worst_eu:.0e}, q = -1 the Minkowski form to {worst_mi:.0e}; "
      f"each breaks the other's by at least {min(cross_eu, cross_mi):.2f}; rescaling for q = {qs} {e_scale:.0e}; "
      f"law minus plain sum at q = 1e-1 ... 1e-4: " + ", ".join(f"{d:.1e}" for d in dev) + f" (slope {slope:.3f})")
results["signature"] = dict(euclid_kept=float(worst_eu), minkowski_kept=float(worst_mi), cross=float(min(cross_eu, cross_mi)),
                            rescale=float(e_scale), small_q=small, deviation=dev, slope=slope)

def law_dist_h(a, b):             # the law's own distance at q = -1 (round 13's formula with q -> -1)
    return 2 * np.arctanh(np.linalg.norm(oplus(-a, b, -1.0), axis=-1))
t1, t2 = rg.uniform(0.2, 3, 2000), rg.uniform(0.2, 3, 2000)
D = law_dist_h(X1, Y1)
interval = np.real(np.linalg.det(t2[:, None, None] * Pm(Y1) - t1[:, None, None] * Pm(X1)))
milne = t1**2 + t2**2 - 2 * t1 * t2 * np.cosh(D)
e_milne = np.abs(interval - milne).max() / np.abs(milne).max()
inside = np.all(np.real(np.linalg.det(t1[:, None, None] * Pm(X1))) > 0) and np.all(to_h3(X1)[:, 0] > 0)
same_x = np.abs(np.real(np.linalg.det((t2 - t1)[:, None, None] * Pm(X1))) - (t2 - t1)**2).max()
check("H5 (exact) special relativity's spacetime is the cone over the law's q = -1 space. Scaling each 4-velocity by a proper time tau fills the future light cone, and the spacetime interval between events (tau1, x) and (tau2, y) is tau1^2 + tau2^2 - 2 tau1 tau2 cosh d, with d the law's own distance from x to y (the Milne form of flat spacetime). Along a fixed x, the interval is the proper time squared",
      e_milne < 1e-12 and inside and same_x < 1e-10,
      f"2000 pairs of events: interval vs Milne form with the law's distance {e_milne:.0e}; all events inside the future cone: {bool(inside)}; along fixed x {same_x:.0e}")
results["milne"] = dict(interval=float(e_milne), inside=bool(inside), along_x=float(same_x))

# for the picture: a row of the law's kicks at q = -1 along one direction, placed on three slices tau = 1, 2, 3
kick = np.array([[0.25, 0.0, 0.0]]); row = [np.zeros((1, 3))]
for _ in range(6):
    row.append(oplus(kick, row[-1], -1.0)); row.insert(0, oplus(-kick, row[0], -1.0))
row = np.vstack(row); H3 = to_h3(row)
gaps = law_dist_h(row[:-1], row[1:])
check("H6 (exact) kicks repeated by the law at q = -1 step along a line of the hyperboloid by equal law-distances (equal rapidity steps), so each slice of the cone carries the law's own lattice, spaced in proportion to tau",
      np.abs(gaps - gaps[0]).max() < 1e-12,
      f"13 points; law-distance between neighbours {gaps[0]:.6f} each, spread {np.abs(gaps - gaps[0]).max():.0e}; rapidity of one kick 2 artanh(0.25) = {2 * math.atanh(0.25):.6f}")
results["row"] = dict(t=H3[:, 0].tolist(), x=H3[:, 1].tolist(), step=float(gaps[0]))

results["summary"] = {"passed": passed, "failed": failed, "seconds": round(time.time() - T_START)}
say(f"\n{passed} passed, {failed} failed  ({time.time() - T_START:.0f} s)")
json.dump(results, open("results.json", "w"), indent=1, default=float)
open("checks_output.txt", "w").write(_log.getvalue())
