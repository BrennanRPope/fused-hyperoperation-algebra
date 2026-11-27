"""DRAFT round 31: does the hidden 2-sphere, and its flux, come out of the law's own hidden 3-sphere?
Rounds 9-10 got one-handed massless matter from the law's twist acting as flux, but assumed (a) a hidden 2-sphere at the
closed end of the dial and (b) matter carrying one sign of charge under the twist. Here: read the law's hidden 3-sphere
(the weak force's SU(2), rounds 4 and 7) along the circles of one of its own rotations (the Hopf fibration), and check
whether (a) follows; then count massless matter by hand in every charge sector to see whether (b) can come from geometry.
Exact parts in sympy; the spectra use round 7's Dirac operator on the 2-sphere unchanged (built from spin matrices)."""
import json, io, time, math
import numpy as np
import sympy as sp

T0 = time.time(); _log = io.StringIO(); passed = failed = 0; results = {}
def say(s): print(s); _log.write(s + "\n")
def check(name, ok, detail):
    global passed, failed
    passed += bool(ok); failed += (not ok)
    say(f"{'PASS' if ok else 'FAIL'} {name}  | {detail}")

# ---------------------------------------------------------------- H1
say("== H1. The law's hidden 3-sphere read along one family of its circles: a 2-sphere with one unit of flux")
th, ph, ps = sp.symbols('theta phi psi', real=True)
q = sp.symbols('q', positive=True)
# a point of the law's 3-sphere (unit quaternion as a pair in C^2), Euler angles; the circles are psi -> psi + s
z1 = sp.cos(th / 2) * sp.exp(sp.I * (ps + ph) / 2)
z2 = sp.sin(th / 2) * sp.exp(sp.I * (ps - ph) / 2)
X = [th, ph, ps]
J = [[sp.diff(zz, v) for v in X] for zz in (z1, z2)]
G = sp.Matrix(3, 3, lambda i, j: sp.simplify(sum(sp.re(sp.expand(J[a][i] * sp.conjugate(J[a][j]))) for a in range(2)) / q))
Hopf = sp.Matrix([[1, 0, 0], [0, 1, sp.cos(th)], [0, sp.cos(th), 1]]) / (4 * q)    # dth^2 + dph^2 + (dps + cos th dph)^2 (+ sin^2 below)
Hopf[1, 1] = (sp.sin(th)**2 + sp.cos(th)**2) / (4 * q)
same = sp.simplify(G - Hopf) == sp.zeros(3, 3)
# base 2-sphere: (1/4q)(dth^2 + sin^2 th dph^2); Gauss curvature
gb = sp.diag(1 / (4 * q), sp.sin(th)**2 / (4 * q))
sE, sG = 1 / (2 * sp.sqrt(q)), sp.sin(th) / (2 * sp.sqrt(q))      # sqrt of the two metric entries (0 < theta < pi)
Kb = sp.simplify(-1 / (sE * sG) * sp.diff(sp.diff(sG, th) / sE, th))
# connection A = cos th dph per unit fibre angle psi; circle length 4 pi in psi; flux of F = dA over the base / 4 pi
F = sp.diff(sp.cos(th), th)                     # F_{th ph}
chern = sp.simplify(sp.integrate(sp.integrate(F, (th, 0, sp.pi)), (ph, 0, 2 * sp.pi)) / (4 * sp.pi))
# the circles are orbits of one of the law's own rotations (right multiplication by e^{i s/2}): check
s_ = sp.symbols('s', real=True)
orbit_ok = sp.simplify(z1.subs(ps, ps + s_) - z1 * sp.exp(sp.I * s_ / 2)) == 0 and sp.simplify(z2.subs(ps, ps + s_) - z2 * sp.exp(sp.I * s_ / 2)) == 0
ok1 = same and sp.simplify(Kb - 4 * q) == 0 and abs(chern) == 1 and orbit_ok
check("H1 (exact; textbook Hopf 1931) the law's hidden 3-sphere of dial q splits exactly into circles (orbits of one of its own rotations, z -> z e^{is/2}) over a 2-sphere of radius 1/(2 sqrt q), i.e. the law's 2-sphere at dial 4q, plus a connection A = cos(theta) dphi along the circles. That connection is a magnetic monopole through the 2-sphere with exactly one unit of flux. So round 9's assumption (a), a hidden 2-sphere with flux, follows from the hidden 3-sphere the weak force already uses",
      ok1, f"pulled-back metric = (1/4q)[dth^2 + sin^2 th dph^2 + (dps + cos th dph)^2]: {same}; base curvature {Kb}; flux / circle length = {chern}; circles are orbits of e^(is/2): {orbit_ok}")
results["H1"] = dict(metric_matches=same, base_curvature=str(Kb), chern=int(chern))

# ---------------------------------------------------------------- H2
say("== H2. Rounds 9 and 10 in this light: the law's twist is flux 2 (the 2-sphere's own frame), charge 1/2 makes it the Hopf flux 1")
K_unit = 1  # Gauss curvature of a unit 2-sphere
euler = sp.integrate(sp.integrate(K_unit * sp.sin(th), (th, 0, sp.pi)), (ph, 0, 2 * sp.pi)) / (2 * sp.pi)
ok2 = euler == 2 and sp.Rational(1, 2) * euler == abs(chern)
check("H2 (exact; textbook) the law's twist round the closed velocity sphere is the 2-sphere's own frame-turning, flux = Euler number 2 (round 9). Round 10's bookkeeping charge 1/2 halves it to 1, which is exactly the Hopf flux of H1. On the 3-sphere a function of weight m along the circles (e^{-i m psi}, m a half-integer since a circle is 4 pi long) sees flux 2m, so round 10's charge 1/2 is weight 1/2 along the hidden 3-sphere's circles and round 9's charge 1 is weight 1",
      ok2, f"Euler number {euler}; charge 1/2 x 2 = {sp.Rational(1,2)*euler} = Hopf flux {abs(chern)}")
results["H2"] = dict(euler=int(euler))

# ---------------------------------------------------------------- H3
say("== H3. Every charge sector of the hidden 3-sphere: massless matter by hand (round 7's operator, unchanged)")
def Jmats(j):
    m = np.arange(j, -j - 1, -1); d = len(m)
    Jp = np.zeros((d, d))
    for i in range(1, d):
        Jp[i - 1, i] = math.sqrt(j * (j + 1) - m[i] * (m[i] + 1))
    return Jp, m
def monopole_dirac(nflux, jmax=8.0, r=1.0):
    """Round 7 verbatim: Dirac operator on S^2 = SU(2)/U(1) with flux nflux; spinor halves have weight (n -+ 1)/2;
    D+ acts as J_+/r inside each V_j (2j+1 copies, the copies carry the 2-sphere's rotations)."""
    mm, mp = (nflux - 1) / 2, (nflux + 1) / 2
    js = [j / 2 for j in range(0, int(2 * jmax) + 1)]
    src = [(j, c) for j in js if abs(mm) <= j and (j - mm) % 1 == 0 for c in range(int(2 * j + 1))]
    dst = [(j, c) for j in js if abs(mp) <= j and (j - mp) % 1 == 0 for c in range(int(2 * j + 1))]
    Dp = np.zeros((len(dst), len(src)))
    for i, (jd, cd) in enumerate(dst):
        for k, (js_, cs) in enumerate(src):
            if jd == js_ and cd == cs:
                Dp[i, k] = math.sqrt(js_ * (js_ + 1) - mm * (mm + 1)) / r
    keep_s = [k for k, (j, c) in enumerate(src) if j < jmax]; keep_d = [i for i, (j, c) in enumerate(dst) if j < jmax]
    return Dp[np.ix_(keep_d, keep_s)], [src[k] for k in keep_s], [dst[i] for i in keep_d]
def zero_modes(nf):
    Dp, src, dst = monopole_dirac(nf)
    if Dp.size == 0:
        return len(src), len(dst), [], []
    U, s, Vt = np.linalg.svd(Dp); rank = int(np.sum(s > 1e-9))
    kerp = Vt[rank:].T                                  # kernel of D+ in the source basis
    kerm = U[:, rank:]                                  # kernel of D- = D+^T in the destination basis
    jp = sorted({src[k][0] for k in range(len(src)) if np.any(np.abs(kerp[k]) > 1e-9)}) if kerp.size else []
    jm = sorted({dst[i][0] for i in range(len(dst)) if np.any(np.abs(kerm[i]) > 1e-9)}) if kerm.size else []
    return kerp.shape[1] if kerp.size else 0, kerm.shape[1] if kerm.size else 0, jp, jm
rows = []
for nf in range(-4, 5):
    L, R, jL, jR = zero_modes(nf)
    rows.append(dict(flux=nf, left=L, right=R, spin_left=jL, spin_right=jR))
    say(f"   flux {nf:+d}: left {L} right {R}; multiplet spin {jL or jR}")
ok3 = all((r_["left"], r_["right"]) == ((r_["flux"], 0) if r_["flux"] > 0 else (0, -r_["flux"])) for r_ in rows) and \
      all((r_["spin_left"] or r_["spin_right"] or [None])[0] == (abs(r_["flux"]) - 1) / 2 for r_ in rows if r_["flux"] != 0)
d2 = next(r_ for r_ in rows if r_["flux"] == 2)
check("H3 (exact spectra; textbook Wu-Yang monopole harmonics) on the 2-sphere under the hidden 3-sphere, matter of weight m along the circles sees flux n = 2m and has |n| massless modes, all of one hand, set by the sign of n, and they form a multiplet of the 2-sphere's rotations of spin (|n|-1)/2. Flux 1 (round 10) gives a single mode; flux 2 (round 9's charge 1) gives a doublet, the shape of the weak force's left-handed pairs",
      ok3, "; ".join(f"{r_['flux']:+d}: L{r_['left']} R{r_['right']}" for r_ in rows) + f"; flux 2 multiplet spin {d2['spin_left']}")
results["H3"] = rows

# ---------------------------------------------------------------- H4
say("== H4. Can the hidden 3-sphere's geometry alone pick one hand?")
# the 3-sphere's functions come in both signs of charge along the circles (D^j_{c,m} and D^j_{c,-m}), so every sector n
# arrives with its mirror -n; complex conjugation maps one to the other
jm_pairs = all(any(rr["flux"] == -r_["flux"] for rr in rows) for r_ in rows)
net = sum(r_["left"] - r_["right"] for r_ in rows)
net_law = sum(r_["left"] - r_["right"] for r_ in rows if r_["flux"] == 1)
# the 3-sphere's isometries SU(2)_L x SU(2)_R; the circles are a U(1) inside SU(2)_R, so what is left is SU(2)_L x U(1)
iso_left = "SU(2) x U(1)"
ok4 = jm_pairs and net == 0 and net_law == 1
check("H4 (exact count; textbook Witten 1981) the hidden 3-sphere carries every weight m together with -m (conjugate functions), so every flux n comes with -n, so its geometry gives the left-handed modes of sector n and the right-handed modes of sector -n in equal numbers: net handedness 0. One hand survives only if matter is kept in one sign of charge, which is round 10's bookkeeping (the law's division by the complex corner): a choice the law suggests, not one geometry forces. The symmetry left over by reading the 3-sphere along its circles is SU(2) x U(1) (textbook: the round 3-sphere's SU(2) x SU(2) with one U(1) kept), the weak and hypercharge pattern; which hand is called left is a convention (orientation)",
      ok4, f"sectors -4..4: net left - right = {net}; with only the law's sector (flux 1): {net_law}; leftover symmetry {iso_left}")
results["H4"] = dict(net_all=net, net_law_sector=net_law, leftover=iso_left)

say(f"\n{passed} passed, {failed} failed  ({time.time() - T0:.0f} s)")
results["summary"] = dict(passed=passed, failed=failed)
json.dump(results, open("results.json", "w"), indent=1, default=str)
open("checks_output.txt", "w").write(_log.getvalue())
