"""DRAFT round 32: is there a reason inside the law to keep one hand of matter?
Round 31 left one assumption under rounds 9-10: matter must keep one sign of weight along the hidden 3-sphere's circles,
because the geometry carries every sector with its mirror. Here we look for anything in the law that breaks the mirror:
(K1) whether the law's own sum prefers left or right; (K2) the Dirac operator on the law's hidden 3-sphere for the whole
f-family of round 7 / torsion-dial (f = 1 the law, f = 0, 2 the flat twists); (K3) its spectral asymmetry (eta invariant);
(K4) whether any member can give net one-handed massless matter once the 3 hidden directions are added to our 4.
Exact parts in sympy; spectra by representation theory with numeric spin matrices; eta by Hurwitz zeta, checked by direct sums."""
import json, io, time
import numpy as np
import sympy as sp
import mpmath as mp

T0 = time.time(); _log = io.StringIO(); passed = failed = 0; results = {}
def say(s): print(s); _log.write(s + "\n")
def check(name, ok, detail):
    global passed, failed
    passed += bool(ok); failed += (not ok)
    say(f"{'PASS' if ok else 'FAIL'} {name}  | {detail}")

# ---------------------------------------------------------------- K1
say("== K1. Does the law's own sum prefer a hand? Its mirror (x -> conjugate) acting on the f-family")
a0, a1, a2, a3, b0, b1, b2, b3 = sp.symbols('a0:4 b0:4', real=True)
Q = sp.Quaternion
# exact check on the half-power form: use x = h^2 with h a general quaternion, so x^(1/2) = h (the f = 1 member is h y h)
h = Q(a0, a1, a2, a3); y = Q(b0, b1, b2, b3)
def conj(u): return Q(u.a, -u.b, -u.c, -u.d)
def eq(u, v): return all(sp.expand(s) == 0 for s in (u.a - v.a, u.b - v.b, u.c - v.c, u.d - v.d))
x = h * h
law = h * y * h                                      # f = 1: x^(1/2) y x^(1/2)
mirror_law = conj(h) * conj(y) * conj(h)             # same rule applied to the mirrored inputs (conj(x) = conj(h)^2)
law_even = eq(conj(law), mirror_law)
f0 = x * y; f2 = y * x                               # the flat twists: all of x first (f = 0), all last (f = 2)
swap_02 = eq(conj(f0), conj(y) * conj(x)) and eq(conj(f2), conj(x) * conj(y))
f0_not_even = not eq(conj(f0), conj(x) * conj(y))
# the mirror reverses orientation of the 3 hidden directions (det of x -> conj on the imaginary part = -1)
det_mirror = sp.Matrix(np.diag([-1, -1, -1])).det()
# numeric: general member x^(1-f/2) y x^(f/2) (torsion-dial T2); mirror sends f -> 2 - f
rng = np.random.default_rng(32)
def qmul(p, r):
    a, b, c, d = p; e, f_, g, k = r
    return np.array([a*e - b*f_ - c*g - d*k, a*f_ + b*e + c*k - d*g, a*g - b*k + c*e + d*f_, a*k + b*g - c*f_ + d*e])
def qpow(p, t):
    p = p / np.linalg.norm(p); ang = np.arccos(np.clip(p[0], -1, 1)); v = p[1:]; nv = np.linalg.norm(v)
    if nv < 1e-15: return np.array([1.0, 0, 0, 0])
    return np.concatenate([[np.cos(t * ang)], np.sin(t * ang) * v / nv])
def qc(p): return p * np.array([1, -1, -1, -1])
def member(xx, yy, f): return qmul(qmul(qpow(xx, 1 - f / 2), yy), qpow(xx, f / 2))
errs = {}
for f in (0.0, 0.3, 1.0, 1.7, 2.0):
    e_same, e_flip = 0, 0
    for _ in range(200):
        xx = rng.normal(size=4); xx /= np.linalg.norm(xx); yy = rng.normal(size=4); yy /= np.linalg.norm(yy)
        lhs = qc(member(xx, yy, f))
        e_same = max(e_same, np.abs(lhs - member(qc(xx), qc(yy), f)).max())
        e_flip = max(e_flip, np.abs(lhs - member(qc(xx), qc(yy), 2 - f)).max())
    errs[f] = (e_same, e_flip)
ok1 = law_even and swap_02 and f0_not_even and det_mirror == -1 and errs[1.0][0] < 1e-12 and all(v[1] < 1e-12 for v in errs.values()) \
      and all(errs[f][0] > 1e-2 for f in (0.0, 0.3, 1.7, 2.0))
check("K1 (exact + numeric; ours on the law, the algebra is textbook) the law's sum x^(1/2) y x^(1/2) reads the same backwards, so the mirror of the hidden 3-sphere (x -> conjugate, which flips all three hidden directions) maps the law to itself: the law has no preferred hand. The mirror sends every other member f to 2 - f, so the flat twists xy (f = 0) and yx (f = 2) are each other's mirror image, and only f = 1 is its own",
      ok1, f"law mirror-even (exact): {law_even}; xy <-> yx under mirror (exact): {swap_02}; mirror det {det_mirror}; "
           f"max |mirror(x+_f y) - (x̄ +_f ȳ)| / |... - (x̄ +_(2-f) ȳ)|: " + ", ".join(f"f={f}: {v[0]:.1e}/{v[1]:.1e}" for f, v in errs.items()))
results["K1"] = dict(law_mirror_even=law_even, flat_twists_swap=swap_02, errs={str(k): v for k, v in errs.items()})

# ---------------------------------------------------------------- K2
say("== K2. The Dirac operator on the law's hidden 3-sphere for every member f: the law's spectrum is symmetric, the others are shifted")
sig = [sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -sp.I], [sp.I, 0]]), sp.Matrix([[1, 0], [0, -1]])]
f = sp.symbols('f', real=True)
# frame on the unit 3-sphere with [e_i, e_j] = 2 eps_ijk e_k; member f: nabla_{e_i} e_j = f eps_ijk e_k (Gamma_ijk = f eps_ijk)
# spinor term in nabla_i psi: (s/4) Gamma_ijk sigma_j sigma_k psi; Dirac D = sigma_i nabla_i. The sign s of this term depends on
# the Clifford sign convention, so we fix it geometrically below (parallel spinors must exist at f = 0 AND f = 2).
Tterm = sp.zeros(2, 2)
for i in range(3):
    for j in range(3):
        for k in range(3):
            Tterm += sp.Rational(1, 4) * f * sp.LeviCivita(i, j, k) * sig[i] * sig[j] * sig[k]
Tterm = sp.simplify(Tterm)                                     # = (3 i f / 2) * identity
torsion_const = sp.simplify(-sp.I * Tterm[0, 0])               # contribution to -iD from s = +1
def Jm(j):
    m = np.arange(j, -j - 1, -1); d = len(m); Jp = np.zeros((d, d))
    for r in range(1, d): Jp[r - 1, r] = np.sqrt(j * (j + 1) - m[r] * (m[r] + 1))
    Jz = np.diag(m); Jx = (Jp + Jp.T) / 2; Jy = (Jp - Jp.T) / (2 * 1j)
    return [Jx, Jy, Jz]
S = [np.array(sg.tolist(), dtype=complex) for sg in sig]
def spectrum(fv, s, jmax=6):
    """-iD_f on functions of the 3-sphere in each sector of spin j (Peter-Weyl; left copy multiplicity 2j+1).
    The frame fields act as e_i = -2i J_i, so -i sigma_i e_i = -2 sigma.J ; spinor term adds -i*s*Tterm."""
    ev = {}; ZS[(fv, s)] = {}
    tc = complex(torsion_const.subs(f, fv))
    for n in range(0, 2 * jmax + 1):
        j = n / 2; J = Jm(j); d = int(round(2 * j + 1))
        M = -2 * sum(np.kron(S[i], J[i]) for i in range(3)) + s * tc * np.eye(2 * d)
        for lam in np.linalg.eigvalsh((M + M.conj().T) / 2):
            key = round(lam, 9); ev[key] = ev.get(key, 0) + d
            if abs(lam) < 1e-7: ZS.setdefault((fv, s), {}).setdefault(j, 0); ZS[(fv, s)][j] += d
    return ev
ZS = {}   # massless modes by sector j, per (f, sign)
def zero_modes(ev): return sum(m for l, m in ev.items() if abs(l) < 1e-7)
# geometric fix of the sign: at f = 0 the frame itself is parallel, so frame-constant spinors (sector j = 0) are parallel;
# at f = 2 the right-handed frame is parallel, and it differs from ours by the 3-sphere's own rotation, which turns spinors
# as spin 1/2, so its constant spinors sit in sector j = 1/2. Each flat twist must have exactly these 2 massless modes.
sign_ok = {}
for s in (+1, -1):
    spectrum(0.0, s); spectrum(2.0, s)
    sign_ok[s] = ZS.get((0.0, s)) == {0.0: 2} and ZS.get((2.0, s)) == {0.5: 2}
s_geo = [s for s in (+1, -1) if sign_ok[s]]
s_ = s_geo[0] if len(s_geo) == 1 else None
# textbook (Baer 1996): Levi-Civita Dirac on the unit 3-sphere has eigenvalues +-(m + 3/2), multiplicity (m+1)(m+2)
ev1 = spectrum(1.0, s_) if s_ else {}
text_ok = all(ev1.get(round(sg * (m + 1.5), 9), 0) == (m + 1) * (m + 2) for m in range(0, 8) for sg in (1, -1))
# general f: the whole spectrum is the law's spectrum moved by c(f) = 3(1 - f)/2 = -(3/4) x torsion 2(f - 1)
shift_ok = True
for fv in (0.0, 0.4, 1.3, 2.0):
    evf = spectrum(fv, s_); c = 1.5 * (1 - fv)
    for m in range(0, 8):
        for sg in (1, -1):
            if evf.get(round(sg * (m + 1.5) + c, 9), 0) != (m + 1) * (m + 2): shift_ok = False
zm0, zm2 = zero_modes(spectrum(0.0, s_)), zero_modes(spectrum(2.0, s_)); zm1 = zero_modes(ev1)
ok2 = sp.simplify(torsion_const - sp.Rational(3, 2) * f) == 0 and s_ is not None and text_ok and shift_ok and zm0 == 2 and zm2 == 2 and zm1 == 0
check("K2 (exact spectra; textbook ingredients: Cartan-Schouten 1926, Baer 1996) on the law's hidden 3-sphere the Dirac operator of member f is the law's own Dirac operator plus the constant 3(1 - f)/2, i.e. -(3/4) x torsion. Only one sign convention gives the flat twists f = 0 and f = 2 exactly their known parallel spinors (frame-constant, in the sectors j = 0 and j = 1/2) and nothing else, and with it the law (f = 1) gives the textbook spectrum +-(m + 3/2), exactly symmetric with no massless mode. f = 0 and f = 2 each have 2 massless modes; every other member's spectrum is lopsided",
      ok2, f"spinor term = {Tterm[0,0]} x 1 -> {torsion_const} in -iD; massless modes by sector for sign +1 / -1 at f = 0: {ZS.get((0.0, 1))} / {ZS.get((0.0, -1))}, at f = 2: {ZS.get((2.0, 1))} / {ZS.get((2.0, -1))}; geometric sign: {s_geo}; "
           f"f = 1 matches +-(m+3/2), mult (m+1)(m+2): {text_ok}; shift 3(1-f)/2 at f = 0, 0.4, 1.3, 2: {shift_ok}; massless modes f=0,1,2: {zm0}, {zm1}, {zm2}")
results["K2"] = dict(shift="3(1-f)/2", sign=s_, massless={"0": zm0, "1": zm1, "2": zm2})

# ---------------------------------------------------------------- K3
say("== K3. Spectral asymmetry (the eta invariant) of each member: zero for the law, odd under the mirror")
c, s = sp.symbols('c s')
a_p, a_m = sp.Rational(3, 2) + c, sp.Rational(3, 2) - c
# eigenvalues +(m + 3/2) + c (x = m + a_p) and -(m + 3/2) + c (|.| = m + a_m), multiplicity (m+1)(m+2) = (x -+ c)^2 - 1/4, |c| < 3/2
def hz0(k, a):   # Hurwitz zeta(-k, a) = -B_{k+1}(a)/(k+1)
    return -sp.bernoulli(k + 1, a) / (k + 1)
eta0 = sp.expand((hz0(2, a_p) - 2 * c * hz0(1, a_p) + (c**2 - sp.Rational(1, 4)) * hz0(0, a_p))
                 - (hz0(2, a_m) + 2 * c * hz0(1, a_m) + (c**2 - sp.Rational(1, 4)) * hz0(0, a_m)))
# check the Hurwitz bookkeeping against a direct sum where it converges (s = 5), at two shifts
def eta_direct(cv, sv, M=200000):
    m = np.arange(M, dtype=float); w = (m + 1) * (m + 2)
    return float(np.sum(w * ((m + 1.5 + cv)**(-sv) - (m + 1.5 - cv)**(-sv))))
def eta_hurwitz(cv, sv):
    ap, am = 1.5 + cv, 1.5 - cv
    return float(mp.zeta(sv - 2, ap) - 2 * cv * mp.zeta(sv - 1, ap) + (cv**2 - .25) * mp.zeta(sv, ap)
                 - (mp.zeta(sv - 2, am) + 2 * cv * mp.zeta(sv - 1, am) + (cv**2 - .25) * mp.zeta(sv, am)))
hz_err = max(abs(eta_direct(cv, 5) - eta_hurwitz(cv, 5)) / abs(eta_hurwitz(cv, 5)) for cv in (0.4, -1.1))
cont_err = max(abs(eta_hurwitz(cv, 0) - float(eta0.subs(c, cv))) for cv in (0.4, -1.1, 1.3))
eta_f = sp.factor(sp.expand(eta0.subs(c, sp.Rational(3, 2) * (1 - f))))
odd = sp.expand(eta0 + eta0.subs(c, -c)) == 0
acc = sp.solve(sp.Eq(eta_f, 0), f)
acc_ok = sorted(acc, key=float) == sorted([1 - 1 / sp.sqrt(3), 1, 1 + 1 / sp.sqrt(3)], key=float)
ok3 = acc_ok and odd and eta0.subs(c, 0) == 0 and hz_err < 1e-4 and cont_err < 1e-10 and all(eta0.subs(c, cv) != 0 for cv in (sp.Rational(1, 2), -1))
check("K3 (exact; APS 1975 eta invariant, our family) the spectral asymmetry of member f is eta(0) = " + str(eta_f) + " for 0 < f < 2. It is odd under the mirror f -> 2 - f and zero for the law (f = 1). It also vanishes at two accidental points f = 1 +- 1/sqrt(3), where the regularised sum cancels though the spectrum itself stays lopsided; everywhere else it is nonzero. So the lopsided members do break the mirror, but as a parity-odd number (a mass-sign / Chern-Simons-type effect), not as massless matter of one hand",
      ok3, f"eta(0) in c = 3(1-f)/2: {sp.factor(eta0)}; odd in c: {odd}; zeros in f: {acc}; Hurwitz vs direct sum at s=5 rel err {hz_err:.1e}; closed form vs continued zeta at s=0 err {cont_err:.1e}")
results["K3"] = dict(eta0_c=str(sp.factor(eta0)), eta0_f=str(eta_f), zeros=[str(z) for z in acc])

# ---------------------------------------------------------------- K4
say("== K4. Our 4 directions plus the 3 hidden ones: no member gives net one-handed massless matter")
I2 = np.eye(2); Z2 = np.zeros((2, 2))
# Euclidean 4D gammas, chiral form
g = [np.block([[Z2, -1j * S[k]], [1j * S[k], Z2]]) for k in range(3)] + [np.block([[Z2, I2], [I2, Z2]])]
g5 = g[0] @ g[1] @ g[2] @ g[3]
if np.allclose(g5 @ g5, -np.eye(4)): g5 = 1j * g5
G7 = [np.kron(gg, I2) for gg in g] + [np.kron(g5, S[i]) for i in range(3)]
cliff = max(np.abs(G7[A] @ G7[B] + G7[B] @ G7[A] - 2 * (A == B) * np.eye(8)).max() for A in range(7) for B in range(7))
P7 = G7[0]
for A in range(1, 7): P7 = P7 @ G7[A]
no_7d_hand = np.allclose(P7, P7[0, 0] * np.eye(8))       # product of all seven is a number: no one-hand condition in 7D
# contrast: our 4 plus a hidden 2-sphere (round 31's base) has 6 gammas whose product is not a number, so a one-hand condition exists
G6 = [np.kron(gg, I2) for gg in g] + [np.kron(g5, S[i]) for i in range(2)]
P6 = G6[0]
for A in range(1, 6): P6 = P6 @ G6[A]
six_has_hand = not np.allclose(P6, P6[0, 0] * np.eye(8)) and np.allclose(P6 @ P6, np.eye(8) * (P6 @ P6)[0, 0])
# a hidden mode with -iD3 chi = lam chi gives, in momentum k, the 4D operator i gamma.k + i lam gamma5: mass^2 = k^2 + lam^2
lam = 0.7; k4 = rng.normal(size=4)
Op = sum(1j * k4[m] * g[m] for m in range(4)) + 1j * lam * g5
mass_ok = np.allclose(np.linalg.eigvalsh(Op.conj().T @ Op), np.dot(k4, k4) + lam**2)
# the 4D massless states are (any 4D spinor) x (hidden massless mode); their net handedness = tr(gamma5) x (number of hidden
# massless modes) = 0, whatever that number is
fs = [0.0, 0.5, 1.0, 1.5, 2.0, 2 + 2 / 3, -2 / 3]
hands = {}
for fv in fs:
    zm = zero_modes(spectrum(fv, s_)); hands[fv] = (zm, int(round(np.trace(g5).real)) * zm)
ok4 = cliff < 1e-12 and no_7d_hand and six_has_hand and mass_ok and all(v[1] == 0 for v in hands.values()) and hands[2 + 2 / 3][0] > 0
check("K4 (exact algebra; textbook: odd dimensions have no Weyl condition, Witten 1981) with the 3 hidden directions added to our 4 the product of all seven gamma matrices is a plain number, so there is no one-hand condition at all. A hidden mode of eigenvalue lam becomes a 4D particle of mass |lam|; a hidden massless mode becomes a full 4D Dirac particle, one left and one right half. So for every member f, including the lopsided ones and the flat twists with their massless modes, net handedness is zero. (With a hidden 2-sphere instead, 6 gammas, a one-hand condition does exist; that is round 31's setting, where it must be chosen.)",
      ok4, f"7D Clifford err {cliff:.0e}; product of seven gammas = {P7[0,0]:.0f} x 1: {no_7d_hand}; product of six is not a number (a hand condition exists in 4+2): {six_has_hand}; mass^2 = k^2 + lam^2: {mass_ok}; "
           "(hidden massless modes, left - right) by f: " + ", ".join(f"{fv:.2f}: {v}" for fv, v in hands.items()))
results["K4"] = dict(product_of_seven=str(P7[0, 0]), net_hand={f"{k:.3f}": v for k, v in hands.items()})

say(f"\n{passed} / {passed + failed} pass  ({time.time() - T0:.1f} s)")
results["summary"] = dict(passed=passed, failed=failed)
json.dump(results, open("results.json", "w"), indent=1, default=str)
open("checks_output.txt", "w").write(_log.getvalue())
