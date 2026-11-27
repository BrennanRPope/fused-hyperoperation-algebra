"""DRAFT round 33: with round 31's 2-sphere reading taken as a premise, does the weak force's left-handed pattern follow?
Premise (stated, not derived; round 32 showed the law can't supply it): matter lives on our 4 directions plus the 2-sphere at
the base of the law's hidden 3-sphere, as 6-dimensional spinors of one chosen hand, and carries weight w along the 3-sphere's
circles (flux n = 2w through the 2-sphere). Tests: (W1) what each sector gives; (W2) whether the measured hypercharges fit;
(W3) whether the gauge anomalies can cancel; (W4) the one number the law's round 3-sphere fixes here, the ratio of the two
couplings, against the measured weak mixing angle. Spectra use round 7's Dirac operator on the 2-sphere, copied verbatim from
round 31; the rest is exact sympy; W4 uses textbook one-loop running."""
import json, io, time, math
import numpy as np
import sympy as sp

T0 = time.time(); _log = io.StringIO(); passed = failed = 0; results = {}
def say(s): print(s); _log.write(s + "\n")
def check(name, ok, detail):
    global passed, failed
    passed += bool(ok); failed += (not ok)
    say(f"{'PASS' if ok else 'FAIL'} {name}  | {detail}")
R_ = sp.Rational

# ---------------------------------------------------------------- W1
say("== W1. What each sector gives under the premise: one multiplet, one hand, and its size tied to its weight")
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
for nf in [n for n in range(-6, 7) if n != 0]:
    L, R, jL, jR = zero_modes(nf)
    j = (jL or jR)[0]
    rows.append(dict(flux=nf, weight=nf / 2, left=L, right=R, spin=j))
# right-handed modes rewritten as left-handed (conjugate): weight -> -w, SU(2) multiplets are their own conjugates
leftform = [(r["spin"], r["weight"] if r["left"] else -r["weight"]) for r in rows]
tie = all(abs(abs(w) - (j + 0.5)) < 1e-12 for j, w in leftform)
one_sign = all(w > 0 for j, w in leftform)
count_ok = all(r["left"] + r["right"] == abs(r["flux"]) and (r["left"] == 0 or r["right"] == 0) for r in rows)
ok1 = tie and one_sign and count_ok
check("W1 (exact spectra; textbook monopole harmonics, round 31's H3) under the premise, sector w gives one multiplet of spin |w| - 1/2, all of one hand (left for w > 0, right for w < 0, with the 6D hand fixed). Written as left-handed fields (a right-handed field is the conjugate of a left-handed one of weight -w), every massless field has w > 0, and its weight is tied to its size: |w| = spin + 1/2. Doublets always have |w| = 1, singlets always |w| = 1/2, and each weight gives exactly one multiplet",
      ok1, "; ".join(f"w={r['weight']:+.1f}: L{r['left']} R{r['right']} spin {r['spin']}" for r in rows) + f"; |w| = spin + 1/2 for all: {tie}; left-form weights all > 0: {one_sign}")
results["W1"] = dict(rows=rows, leftform=leftform)

# ---------------------------------------------------------------- W2
say("== W2. The measured hypercharges against the premise (Y = c x w)")
# textbook Standard Model, one generation, all as left-handed Weyl fields: (name, colours, SU(2) spin, Y)
SM = [("Q", 3, R_(1, 2), R_(1, 6)), ("u^c", 3, 0, R_(-2, 3)), ("d^c", 3, 0, R_(1, 3)), ("L", 1, R_(1, 2), R_(-1, 2)), ("e^c", 1, 0, R_(1))]
ratio = {n: Y / (j + R_(1, 2)) for n, _, j, Y in SM}            # premise needs Y / (spin + 1/2) = c for all (one 6D hand)
absratio = {n: abs(v) for n, v in ratio.items()}                # or |Y| / (spin + 1/2) = |c| if either 6D hand were allowed
c = sp.symbols('c')
sol_one = sp.solve([sp.Eq(Y, c * (j + R_(1, 2))) for _, _, j, Y in SM], c)
sol_any = sp.solve([sp.Eq(abs(Y), c * (j + R_(1, 2))) for _, _, j, Y in SM], c)
# right-handed neutrino (Y = 0) would need w = 0, which gives no massless mode at all
nu_R = [r for r in rows if r["flux"] == 0] == []
# best single c for leptons only (L, e^c): fails too
sol_lep = sp.solve([sp.Eq(Y, c * (j + R_(1, 2))) for n, _, j, Y in SM if n in ("L", "e^c")], c)
ok2 = sol_one == [] and sol_any == [] and sol_lep == [] and len(set(absratio.values())) == 5 and nu_R
check("W2 (exact; Standard Model hypercharges are textbook) the premise ties each field's hypercharge to its size, Y / (spin + 1/2) = c, one number for all. The measured pattern breaks that tie: the two doublets Q and L have different hypercharges (1/6 and -1/2) though both are doublets, and the singlets' hypercharges are 2 to 4 times their doublet's (e^c to L: 2; u^c and d^c to Q: 4 and 2), where the premise gives 1/2. No c fits, not even for the leptons alone or allowing either 6D hand per field. So hypercharge is not the circle weight: it needs a charge of its own",
      ok2, "Y/(spin+1/2): " + ", ".join(f"{n} {v}" for n, v in ratio.items()) + f"; solutions for c: one 6D hand {sol_one}, either hand {sol_any}, leptons only {sol_lep}; neutral singlet impossible (w = 0 has no massless mode): {nu_R}")
results["W2"] = dict(ratio={k: str(v) for k, v in ratio.items()})

# ---------------------------------------------------------------- W3
say("== W3. Gauge anomalies: can the premise's matter be consistent?")
def anomalies(fields):          # fields: (multiplicity, SU(2) spin, Y), all left-handed Weyl
    dimj = lambda j: int(2 * j + 1)
    T2 = lambda j: sp.Rational(int(2 * j + 1) * j * (j + 1), 3)     # Dynkin index sum of T3^2 over the multiplet
    trY = sum(m * dimj(j) * Y for m, j, Y in fields)
    trY3 = sum(m * dimj(j) * Y**3 for m, j, Y in fields)
    trT2Y = sum(m * T2(j) * Y for m, j, Y in fields)
    doublets = sum(m for m, j, Y in fields if j == R_(1, 2))
    return dict(grav=trY, Y3=trY3, T2Y=trT2Y, doublets=doublets, witten_ok=(doublets % 2 == 0))
sm = anomalies([(col, j, Y) for _, col, j, Y in SM])
cc = sp.symbols('c', positive=True)
# premise: any content from one 6D hand, all left-form weights w > 0, Y = c w. Generic: n_k copies of weight k/2 (spin (k-1)/2).
nks = sp.symbols('n1:5', nonnegative=True, integer=True)
prem = anomalies([(nks[k], R_(k, 2), cc * R_(k + 1, 2)) for k in range(4)])     # n1 singlets (w = 1/2), n2 doublets (w = 1), ...
positive = all(sp.simplify(prem[key]).is_nonnegative for key in ("grav", "Y3", "T2Y"))
round9 = anomalies([(1, R_(1, 2), cc * 1)])                       # round 9's flux-2 sector alone: one doublet
ok3 = sm["grav"] == 0 and sm["Y3"] == 0 and sm["T2Y"] == 0 and sm["witten_ok"] and positive and not round9["witten_ok"] \
      and all(sp.simplify(round9[k]) != 0 for k in ("grav", "Y3", "T2Y"))
check("W3 (exact; textbook anomaly conditions, Witten's global SU(2) anomaly 1982) the Standard Model's generation cancels every anomaly (Tr Y, Tr Y^3, Tr T^2 Y all 0) with an even number of doublets (4). Under the premise every left-handed field has hypercharge of one sign, so Tr Y, Tr Y^3 and Tr T^2 Y are sums of positive terms: they cancel only if there is no matter at all. Round 9's flux-2 sector alone is one doublet, an odd number, so it also has Witten's global anomaly",
      ok3, f"SM: Tr Y {sm['grav']}, Tr Y^3 {sm['Y3']}, Tr T^2 Y {sm['T2Y']}, doublets {sm['doublets']}; premise (n_k copies of weight k/2): Tr Y = {sp.factor(prem['grav'])}, Tr T^2 Y = {sp.factor(prem['T2Y'])}, all >= 0: {positive}; round 9 alone: {round9}")
results["W3"] = dict(sm={k: str(v) for k, v in sm.items()}, premise={k: str(v) for k, v in prem.items()}, round9={k: str(v) for k, v in round9.items()})

# ---------------------------------------------------------------- W4
say("== W4. The one number the law's round 3-sphere fixes here: the two couplings' ratio, against the weak mixing angle")
# Killing vectors of the round 3-sphere: left multiplication by exp(i s sigma_a / 2) (the weak SU(2)) and right multiplication by
# exp(i s sigma_3 / 2) (the circles, weight w = T3 of the right SU(2)). In a Kaluza-Klein reading each coupling is set by the
# length of its Killing vector, 1/g^2 proportional to its mean square (textbook: Weinberg 1983). Compare the lengths.
x0, x1, x2, x3 = sp.symbols('x0:4', real=True)
qv = sp.Matrix([[x0 + sp.I * x3, x2 + sp.I * x1], [-x2 + sp.I * x1, x0 - sp.I * x3]])   # unit quaternion as SU(2) matrix
sg = [sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -sp.I], [sp.I, 0]]), sp.Matrix([[1, 0], [0, -1]])]
def sqlen(M):   # |M|^2 in the round metric on the unit 3-sphere: (1/2) tr(M M^dagger)
    return sp.simplify(sp.expand((M * M.H).trace() / 2))
unit = {x0**2: 1 - x1**2 - x2**2 - x3**2}
lens_L = [sp.simplify(sqlen(sp.I * sg[a] / 2 * qv).subs(unit)) for a in range(3)]
len_R = sp.simplify(sqlen(qv * sp.I * sg[2] / 2).subs(unit))
equal = all(sp.simplify(l - len_R) == 0 for l in lens_L)
sin2_law = R_(1, 2) if equal else None                          # g' = g with Y = T3_R  ->  sin^2 = g'^2/(g^2 + g'^2) = 1/2
# measured at the Z mass (PDG review values, approximate): 1/alpha = 127.95, sin^2 theta_W = 0.2312 (MS-bar)
MZ, ainv, s2 = 91.1876, 127.95, 0.2312
a2inv, aYinv = ainv * s2, ainv * (1 - s2)
bY, b2 = 41 / 6, -19 / 6                                        # textbook one-loop Standard Model coefficients (Y unnormalised)
def scale(aYi, a2i): return MZ * math.exp(2 * math.pi * (aYi - a2i) / (bY - b2))
mu_star = scale(aYinv, a2inv)
spread = [scale(ainv * (1 + e1) * (1 - s2 * (1 + e2)), ainv * (1 + e1) * s2 * (1 + e2)) for e1 in (-.005, .005) for e2 in (-.005, .005)]
MPl = 1.22e19
ok4 = equal and mu_star > MPl and min(spread) > MPl
check("W4 (exact lengths; textbook Kaluza-Klein couplings and one-loop running) on the law's round 3-sphere the weak SU(2)'s Killing vectors and the circles' Killing vector have exactly the same length, so if hypercharge were the circle weight the two couplings would be equal and sin^2 theta_W = 1/2 where the hidden space opens up. Measured, sin^2 theta_W = 0.231 at the Z mass. Running the Standard Model's couplings up, they become equal only at about 5 x 10^20 GeV, some 40 times above the Planck mass, where the running can't be trusted. So even this reading fails, independently of W2",
      ok4, f"squared lengths: weak SU(2) {lens_L}, circles {len_R}; equal: {equal}; sin^2 (law) = {sin2_law}; equal couplings at {mu_star:.2e} GeV (inputs +-0.5%: {min(spread):.1e} to {max(spread):.1e}); Planck {MPl:.2e} GeV")
results["W4"] = dict(equal_lengths=equal, sin2_law=str(sin2_law), mu_star=mu_star, spread=[min(spread), max(spread)], a2inv=a2inv, aYinv=aYinv, MPl=MPl)

say(f"\n{passed} / {passed + failed} pass  ({time.time() - T0:.1f} s)")
results["summary"] = dict(passed=passed, failed=failed)
json.dump(results, open("results.json", "w"), indent=1, default=str)
open("checks_output.txt", "w").write(_log.getvalue())
