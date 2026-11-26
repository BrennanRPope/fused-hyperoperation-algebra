"""Second pass (Nov 26 2026): does each associahedron equation have its own random-geodesic probability meaning?

Model (known, from the probability-entanglement page): a Poisson soup of hyperbolic geodesics in the disk with the
Crofton (kinematic) intensity dtheta dphi / (4 sin^2((theta-phi)/2)) on endpoint pairs, at rate lam (lam = 1 is the
page's choice). Two boundary arcs are "linked" if some geodesic has one end in each.
Run: python3 checks2.py > checks2_output.txt   (writes sim2.npz for figures2.py)
"""
import itertools, math, json
import numpy as np

rng = np.random.default_rng(1126)
passed = total = 0
results = {}


def check(name, ok, detail=""):
    global passed, total
    total += 1; passed += bool(ok)
    print(("PASS " if ok else "FAIL ") + name + ("  | " + detail if detail else ""))
    results[name] = {"ok": bool(ok), "detail": detail}


def soup(U, eps, lam=1.0):
    ce = 1 / math.tan(eps / 2)
    M = lam * math.pi * ce
    N = rng.poisson(M, U)
    ids = np.repeat(np.arange(U), N)
    tot = int(N.sum())
    # ordered pairs (theta, theta + d), d in (eps, 2pi - eps), density (1/2) * 1/(4 sin^2(d/2)) per unit theta;
    # CDF in d: (cot(eps/2) - cot(d/2)) / (2 cot(eps/2))
    w = rng.uniform(0, 1, tot)
    cot_half = ce - 2 * ce * w
    d = 2 * np.arctan2(1.0, cot_half)                # arccot, in (0, 2pi)
    th = rng.uniform(0, 2 * np.pi, tot)
    ph = (th + d) % (2 * np.pi)
    return ids, th, ph


def side_of(ang, alphas):
    """index k of the arc [alpha_k, alpha_{k+1}) containing ang (alphas sorted, wraps)."""
    k = np.searchsorted(alphas, ang, side="right") - 1
    return k % len(alphas)


def linked_matrix(U, n, ids, s1, s2):
    L = np.zeros((U, n, n), bool)
    L[ids, s1, s2] = True
    L[ids, s2, s1] = True
    return L


def zc(a):
    return np.exp(1j * a)


def u_arcs(x1, x2, y1, y2):
    return ((x1 - y2) * (x2 - y1) / ((x1 - y1) * (x2 - y2))).real


def chords(n):
    return [(i, j) for i in range(n) for j in range(i + 2, n) if not (i == 0 and j == n - 1)]


def between(i, j, n):
    """sides strictly between i and j going forward from i."""
    return [(i + t) % n for t in range(1, (j - i) % n)]


def crosses(X, Y, n):
    i, j = X; a, b = Y
    bi = set(between(i, j, n)); bo = set(between(j, i, n))
    return (a in bi and b in bo) or (a in bo and b in bi)


# ---------------------------------------------------------------- the polygon
n = 6
while True:
    alphas = np.sort(rng.uniform(0, 2 * np.pi, n))
    gaps = np.diff(np.r_[alphas, alphas[0] + 2 * np.pi])
    if gaps.min() > 0.25:
        break
z = zc(alphas)
U, eps = 20000, 0.02
ids, th, ph = soup(U, eps)
L = linked_matrix(U, n, ids, side_of(th, alphas), side_of(ph, alphas))
se = lambda p: math.sqrt(max(p * (1 - p), 1e-9) / U)
print(f"hexagon arc lengths {np.round(gaps, 3)}; {U} universes, {len(th)/U:.1f} geodesics each on average")

# P1 single channels
print("\n== P. Each coordinate as a probability")
CH = chords(n)
exact_u = {}; sim_u = {}
worst = 0
for (i, j) in CH:
    exact_u[(i, j)] = u_arcs(z[i], z[(i + 1) % n], z[j], z[(j + 1) % n])
    sim_u[(i, j)] = 1 - L[:, i, j].mean()
    worst = max(worst, abs(sim_u[(i, j)] - exact_u[(i, j)]) / se(exact_u[(i, j)]))
check("P1 u_ij = P(no geodesic links side i to side j), all 9 chords of the hexagon (known: the probability page's clock)",
      worst < 4.5, f"worst deviation {worst:.2f} standard errors")

# P2 gluing = independence, P3 u-equation from frequencies, P4 both = neither
w_ind = w_eq = w_bn = 0
rows = []
for X in CH:
    i, j = X
    Ys = [Y for Y in CH if crosses(X, Y, n)]
    b = between(i, j, n); d = between(j, i, n)
    block_link = L[:, b][:, :, d].any(axis=(1, 2))
    p_noblock = 1 - block_link.mean()
    prod = np.prod([sim_u[Y] for Y in Ys])
    w_ind = max(w_ind, abs(p_noblock - prod) / se(p_noblock))
    no_ac = ~L[:, i, j]
    w_eq = max(w_eq, abs(no_ac.mean() + p_noblock - 1) / (se(no_ac.mean()) + se(p_noblock)))
    both = (L[:, i, j] & block_link).mean(); neither = (no_ac & ~block_link).mean()
    w_bn = max(w_bn, abs(both - neither) / (se(both) + se(neither)))
    rows.append({"chord": f"{i}{j}", "both": both, "neither": neither, "exactly_one": 1 - both - neither,
                 "u": exact_u[X], "u_block": 1 - exact_u[X]})
check("P2 P(the two blocks are unlinked) = product of the crossing u's: the channels are independent (Poisson), "
      "so gluing is our union rule as an event", w_ind < 4.5, f"worst {w_ind:.2f} s.e.")
check("P3 each u-equation, read off simulated frequencies: P(no link i-j) + P(no link between the blocks) = 1",
      w_eq < 4.5, f"worst {w_eq:.2f} s.e.")
check("P4 the probability meaning (ours): for every chord, 'both crossing channels linked' and 'neither linked' are "
      "equally likely", w_bn < 4.5, f"worst {w_bn:.2f} s.e.; e.g. chord {rows[0]['chord']}: both {rows[0]['both']:.4f}, "
      f"neither {rows[0]['neither']:.4f}")
# exact equivalence (algebra): given independence, P(both) = P(neither) <=> u + v = 1
u_, v_ = np.meshgrid(np.linspace(0.01, .99, 99), np.linspace(.01, .99, 99))
eq = np.isclose((1 - u_) * (1 - v_), u_ * v_); line = np.isclose(u_ + v_, 1)
check("P5 with independent channels, 'P(both) = P(neither)' holds exactly when u + v = 1: the reading is equivalent to "
      "the u-equation, not weaker", np.array_equal(eq, line))

# R rate test
print("\n== R. The equations pin the clock's rate")
X = CH[0]; i, j = X
b = between(i, j, n); d = between(j, i, n)
rate_rows = []
okR = True
for lam in (0.5, 0.75, 1.0, 1.5, 2.0):
    Ur = 12000
    ids_r, th_r, ph_r = soup(Ur, eps, lam)
    Lr = linked_matrix(Ur, n, ids_r, side_of(th_r, alphas), side_of(ph_r, alphas))
    s_sim = (1 - Lr[:, i, j].mean()) + (1 - Lr[:, b][:, :, d].any(axis=(1, 2)).mean())
    s_ex = exact_u[X] ** lam + (1 - exact_u[X]) ** lam
    rate_rows.append((lam, s_sim, s_ex))
    okR &= abs(s_sim - s_ex) < 0.02
print("   rate, simulated sum, exact u^lam + (1-u)^lam:", [(l, round(a, 4), round(e, 4)) for l, a, e in rate_rows])
check("R1 at rate lam the sum is u^lam + (1-u)^lam, matching the simulation", okR)
check("R2 the sum is 1 only at rate 1: above 1 for slower clocks, below 1 for faster (the associahedron fixes the "
      "clock rate that the probability page chose by hand; ours)",
      all((e > 1) == (l < 1) for l, _, e in rate_rows if l != 1) and abs(rate_rows[2][2] - 1) < 1e-12)

# C certainty: facets of the associahedron
print("\n== C. Facets as certainty: crossing channels exclude each other's certainty")
cert = []
for s1len in (0.4, 0.15, 0.05, 0.015):
    base = np.array([0.0, 1.0, 1.0 + s1len, 2.6, 3.9, 5.1])
    zz = zc(base)
    u02 = u_arcs(zz[0], zz[1], zz[2], zz[3])                 # sides 0 and 2, side 1 between
    u_cross = [u_arcs(zz[1], zz[2], zz[k], zz[(k + 1) % 6]) for k in (3, 4, 5)]
    e2 = min(0.5 * s1len, 0.02)
    Uc = 3000
    ids_c, th_c, ph_c = soup(Uc, e2)
    Lc = linked_matrix(Uc, 6, ids_c, side_of(th_c, base), side_of(ph_c, base))
    p02 = Lc[:, 0, 2].mean(); p1x = Lc[:, 1][:, [3, 4, 5]].any(axis=1).mean()
    cert.append((s1len, 1 - u02, 1 - np.prod(u_cross), p02, p1x))
print("   side-1 length, P(link 0-2) exact, P(side 1 linked to 3,4,5) exact, simulated both:",
      [tuple(round(x, 4) for x in c) for c in cert])
check("C1 as the side between them shrinks, sides 0 and 2 become surely linked (u -> 0, a facet of the associahedron) "
      "and every crossing channel becomes surely unlinked",
      cert[-1][1] > 0.97 and cert[-1][2] < 0.03 and abs(cert[-1][3] - cert[-1][1]) < 0.03)
check("C2 so two crossing channels can never both be certain, which is why the associahedron's corners are the "
      "non-crossing sets of chords (triangulations): the maximal sets of channels that can be certain together (reading ours)",
      all(abs((1 - a) + (1 - bb) - 1) < 1e-12 for _, a, bb, _, _ in cert))

# N larger polygons, exact only
print("\n== N. Same meaning for larger polygons (exact, from the cross-ratios)")
okN = True
for nn in (7, 8, 9):
    al = np.sort(rng.uniform(0, 2 * np.pi, nn)); zn = zc(al)
    for Xn in chords(nn):
        i2, j2 = Xn
        u = u_arcs(zn[i2], zn[(i2 + 1) % nn], zn[j2], zn[(j2 + 1) % nn])
        v = np.prod([u_arcs(zn[a], zn[(a + 1) % nn], zn[bq], zn[(bq + 1) % nn]) for (a, bq) in chords(nn) if crosses(Xn, (a, bq), nn)])
        okN &= abs((1 - u) * (1 - v) - u * v) < 1e-12
check("N1 'both equally likely as neither' holds for every chord of 7-, 8- and 9-gons", okN)

np.savez("sim2.npz", alphas=alphas, rate=np.array(rate_rows), cert=np.array(cert),
         su=np.array([sim_u[c] for c in CH]), eu=np.array([exact_u[c] for c in CH]),
         chords=np.array(CH), th=th[ids < 1], ph=ph[ids < 1],
         both=np.array([r["both"] for r in rows]), neither=np.array([r["neither"] for r in rows]),
         one=np.array([r["exactly_one"] for r in rows]), U=U)
print(f"\n{passed}/{total} checks passed")
json.dump({"passed": passed, "total": total, "results": results, "both_neither": rows}, open("results2.json", "w"), indent=1)
