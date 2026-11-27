"""DRAFT round 9 checks: can the law's own twist act as flux through a hidden 2-sphere?

Round 7: twist makes massless matter only in left-right pairs; FLUX through a hidden 2-sphere makes one hand.
Here the flux is taken from the twisted-boost law itself (../twisted-boosts/): dial boosts B_q(u) = [[1, u], [-q u*, 1]]
on half-velocities compose as u (+) v = (u + v)/(1 - q u* v) with a twist e^{i w} = (1 - q u v*)/(1 - q u* v),
and w = q x (signed area) with curvature K = q.  At the rotation end of the dial (q = +1) velocity space closes into a
sphere; at relativity's end (q = -1) it is an open disc.
Run: python3 checks.py  (numpy).  Writes results.json and checks_output.txt.
"""
import json, io, math, itertools
import numpy as np

results, passed, failed = {}, 0, 0
_log = io.StringIO()
def say(s):
    print(s); _log.write(s + "\n")
def check(name, ok, detail=""):
    global passed, failed
    passed, failed = (passed + 1, failed) if ok else (passed, failed + 1)
    say(("PASS " if ok else "FAIL ") + name + (("  | " + detail) if detail else ""))

def oplus(u, v, q):
    return (u + v) / (1 - q * np.conj(u) * v)
def twist_angle(u, v, q):
    return np.angle((1 - q * u * np.conj(v)) / (1 - q * np.conj(u) * v))
def tri_twist(a, b, c, q):
    """twist of the velocity triangle (a, b, c): move a to rest with the law (u = (-a)(+)b, w = (-a)(+)c),
    then the triangle (rest, u, u(+)v) with v = (-u)(+)w has twist angle(u, v)"""
    u = oplus(-a, b, q); w = oplus(-a, c, q); v = oplus(-u, w, q)
    return twist_angle(u, v, q)

# ------------------------------------------------------------------ geometry helpers (independent of the law)
def sph_excess(A, B, C):
    """signed spherical excess (area) of the geodesic triangle A, B, C on the unit sphere (Van Oosterom-Strackee)"""
    num = np.dot(A, np.cross(B, C)); den = 1 + np.dot(A, B) + np.dot(B, C) + np.dot(C, A)
    return 2 * math.atan2(num, den)
def icosphere(level):
    t = (1 + 5**0.5) / 2
    V = [(-1, t, 0), (1, t, 0), (-1, -t, 0), (1, -t, 0), (0, -1, t), (0, 1, t), (0, -1, -t), (0, 1, -t),
         (t, 0, -1), (t, 0, 1), (-t, 0, -1), (-t, 0, 1)]
    V = [np.array(v, float) / np.linalg.norm(v) for v in V]
    F = [(0, 11, 5), (0, 5, 1), (0, 1, 7), (0, 7, 10), (0, 10, 11), (1, 5, 9), (5, 11, 4), (11, 10, 2), (10, 7, 6),
         (7, 1, 8), (3, 9, 4), (3, 4, 2), (3, 2, 6), (3, 6, 8), (3, 8, 9), (4, 9, 5), (2, 4, 11), (6, 2, 10), (8, 6, 7), (9, 8, 1)]
    for _ in range(level):
        cache = {}; F2 = []
        def mid(i, j):
            key = (min(i, j), max(i, j))
            if key not in cache:
                m = V[i] + V[j]; V.append(m / np.linalg.norm(m)); cache[key] = len(V) - 1
            return cache[key]
        for a, b, c in F:
            ab, bc, ca = mid(a, b), mid(b, c), mid(c, a)
            F2 += [(a, ab, ca), (b, bc, ab), (c, ca, bc), (ab, bc, ca)]
        F = F2
    return np.array(V), F
def chart_north(p):    # stereographic from the south pole: |u| = tan(theta/2), metric 2|du|/(1+|u|^2), K = +1
    return complex(p[0], p[1]) / (1 + p[2])
def chart_south(p):    # from the north pole, orientation kept: u' = 1/u on the overlap
    return complex(p[0], -p[1]) / (1 - p[2])

say("== R. The law's twist at the rotation end (q = +1) is a closed sphere's worth of flux")
rng = np.random.default_rng(9)
Rm = np.linalg.qr(rng.normal(size=(3, 3)))[0]       # random orientation, so no vertex sits on a pole
V, F = icosphere(4)
V = V @ Rm.T
errs, total, per_face = [], 0.0, []
for a, b, c in F:
    A, B, C = V[a], V[b], V[c]
    zc = (A + B + C)[2]
    ch = chart_north if zc > 0 else chart_south      # use the chart away from its projection point
    w = tri_twist(ch(A), ch(B), ch(C), 1.0)
    e = sph_excess(A, B, C)
    errs.append(abs(w - e)); total += w; per_face.append((float(zc), float(w), float(e)))
check("R1 every one of 5120 small velocity triangles: the law's twist = its area on the unit sphere (K = q = +1)",
      max(errs) < 1e-12, f"worst {max(errs):.1e}")
check("R2 total twist over the closed velocity sphere = 4 pi: a flux of 2 whole turns (2 x 2 pi)",
      abs(total - 4 * math.pi) < 1e-9, f"total / 2pi = {total / (2 * math.pi):.12f}")
# the two charts disagree by a turn that winds twice round the equator: that winding IS the flux number
phis = np.linspace(0, 2 * math.pi, 2001)
frame = []
for ph in phis:
    p = np.array([math.cos(ph), math.sin(ph), 0.0]); u = chart_north(p)
    frame.append(np.angle(-1 / u**2))                    # d u' / d u for u' = 1/u: how a velocity arrow turns between charts
wind = np.sum(np.angle(np.exp(1j * np.diff(np.unwrap(frame))))) / (2 * math.pi)
check("R3 the twist can't be smoothed away: going round the equator, the two charts' frames differ by a turn that winds 2 times",
      abs(abs(wind) - 2) < 1e-9, f"winding {wind:+.6f}")
# relativity's end: open disc, the twist piles up without bound
rows_disc = []
nt = 120
for R in (0.5, 0.8, 0.9, 0.95, 0.99):
    nr = 40
    tot = 0.0
    for i in range(nr):
        r0, r1 = R * i / nr, R * (i + 1) / nr
        for j in range(nt):
            t0, t1 = 2 * math.pi * j / nt, 2 * math.pi * (j + 1) / nt
            p00, p01 = r0 * np.exp(1j * t0), r0 * np.exp(1j * t1); p10, p11 = r1 * np.exp(1j * t0), r1 * np.exp(1j * t1)
            tot += tri_twist(p00, p10, p11, -1.0) + (tri_twist(p00, p11, p01, -1.0) if i > 0 else 0.0)
    # the triangles tile the geodesic 120-gon with corners at |u| = R; its area (textbook): (n-2) pi - n alpha,
    # with cosh(rho) = cot(pi/n) cot(alpha/2) and rho = 2 artanh(R) the corners' distance from rest
    rho = 2 * math.atanh(R); alpha = 2 * math.atan(1 / (math.cosh(rho) * math.tan(math.pi / nt)))
    poly = (nt - 2) * math.pi - nt * alpha
    disc = 4 * math.pi * R**2 / (1 - R**2)              # area of the whole disc |u| < R
    rows_disc.append((R, tot, -poly, -disc))
check("R4 relativity's end (q = -1): the law's twist over a patch of the open velocity disc = minus its area, exactly; the disc's area grows without bound toward light speed, so there is no whole-number flux",
      all(abs(t - e) < 1e-9 for _, t, e, _ in rows_disc),
      ", ".join(f"|u|<{R}: disc {d / (2*math.pi):+.1f} turns" for R, _, _, d in rows_disc))

say("\n== T. Matter that feels the twist as a charge, on the hidden velocity sphere")
r = 1.0
def monopole_dirac(nflux, jmax=8.0):
    """Dirac operator on S^2 with extra flux nflux (as in round 7): D+ = J_+/r between weights (n-1)/2 and (n+1)/2"""
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
    Dp = monopole_dirac(n)
    s = np.linalg.svd(Dp, compute_uv=False) if Dp.size else np.array([])
    rank = int(np.sum(s > 1e-9))
    return Dp.shape[1] - rank, Dp.shape[0] - rank
flux_number = round(total / (2 * math.pi))           # 2, from R2
rowsT = []
for c in (0.0, 0.5, 1.0, -0.5):
    n = int(round(c * flux_number)); L, Rr = counts(n)
    rowsT.append((c, n, L, Rr))
check("T1 matter feeling the twist only through geometry (its spin, charge 0): no massless modes, as round 7",
      rowsT[0][2:] == (0, 0))
check("T2 matter with a charge c under the law's twist (flux 2c): 2c massless modes, all one hand (c = 1/2: one; c = 1: two)",
      rowsT[1][2:] == (1, 0) and rowsT[2][2:] == (2, 0), "; ".join(f"c={c:+.1f}: L{L} R{Rr}" for c, n, L, Rr in rowsT))
check("T3 the opposite charge gives the opposite hand", rowsT[3][2:] == (0, 1))

say(f"\n{passed} passed, {failed} failed")
results.update({"faces": len(F), "worst_face_err": max(errs), "total_twist": total, "flux_number": flux_number,
                "winding": wind, "disc": rows_disc, "charges": rowsT,
                "per_face_sample": per_face[::40], "summary": {"passed": passed, "failed": failed}})
json.dump(results, open("results.json", "w"), indent=1)
open("checks_output.txt", "w").write(_log.getvalue())
