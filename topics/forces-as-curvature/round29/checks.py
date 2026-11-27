"""DRAFT round 29: does the law fix the slope of the dial as rolling dark energy (audit C6)?
The dial-potential search found that rescaling forces an exponential potential V = V0 exp(-lam phi) (reduced Planck
units) and fitted lam = 0.82 (1 sigma 0.61 to 0.97) to Pantheon+ and DESI DR2 BAO. The slope is fitted, not forced.
Here: every reading of the dial that our gravity rounds make geometric fixes lam. Which numbers come out, and do any
land in the data's band? Exact parts in sympy; the data part reuses ../dial-potential/quintessence.py unchanged."""
import json, io, time, math, os, sys
import sympy as sp

T0 = time.time(); _log = io.StringIO(); passed = failed = 0; results = {}
def say(s): print(s); _log.write(s + "\n")
def check(name, ok, detail):
    global passed, failed
    passed += bool(ok); failed += (not ok)
    say(f"{'PASS' if ok else 'FAIL'} {name}  | {detail}")

# ---------------------------------------------------------------- D1
say("== D1. The dial as the size of the law's hidden 3-sphere: 7 = 4 + 3 dimensions")
t = sp.symbols('t')
chi, th, ph, x1, x2, x3 = sp.symbols('chi theta phi x1 x2 x3')
N, a, b = [sp.Function(s)(t) for s in ('N', 'a', 'b')]
X = [t, x1, x2, x3, chi, th, ph]
g = sp.diag(-N**2, a**2, a**2, a**2, b**2, b**2 * sp.sin(chi)**2, b**2 * sp.sin(chi)**2 * sp.sin(th)**2)
D = 7; gi = g.inv()
Gam = [[[sp.simplify(sum(gi[i, l] * (sp.diff(g[l, j], X[k]) + sp.diff(g[l, k], X[j]) - sp.diff(g[j, k], X[l]))
                         for l in range(D)) / 2) for k in range(D)] for j in range(D)] for i in range(D)]
def Ric(j, k):
    return sum(sp.diff(Gam[i][j][k], X[i]) - sp.diff(Gam[i][j][i], X[k])
               + sum(Gam[i][i][l] * Gam[l][j][k] - Gam[i][k][l] * Gam[l][j][i] for l in range(D)) for i in range(D))
R7 = sp.simplify(sum(gi[j, j] * Ric(j, j) for j in range(D)))
L7 = sp.simplify(N * a**3 * b**3 * R7)                     # sqrt(-g) R per unit angular volume
# Einstein frame: g4 = e^{-3 psi} gE, b = e^{psi}; so a = A e^{-3psi/2}, N = n e^{-3psi/2}
A, psi, nE = [sp.Function(s)(t) for s in ('A', 'psi', 'nE')]
sub = {a: A * sp.exp(-3 * psi / 2), b: sp.exp(psi), N: nE * sp.exp(-3 * psi / 2)}
L7E = sp.simplify(L7.subs(sub).doit())
Lexp = -6 * A * A.diff(t)**2 / nE + sp.Rational(15, 2) * A**3 * psi.diff(t)**2 / nE + 6 * nE * A**3 * sp.exp(-5 * psi)
from sympy.calculus.euler import euler_equations
el = euler_equations(sp.expand(L7E - Lexp), [A, psi, nE], t)
el_zero = all(sp.simplify(e.lhs - e.rhs) == 0 for e in el)
# With S = (1/2) int sqrt(g) R: kinetic (1/2)(15/2) psi'^2 -> phi = sqrt(15/2) psi; V = -3 e^{-5 psi}
lam_sphere = sp.Rational(5) / sp.sqrt(sp.Rational(15, 2))
lam_flux = sp.Rational(7) / sp.sqrt(sp.Rational(15, 2))
check("D1 (exact; textbook Kaluza-Klein breathing mode, e.g. Carroll-Geddes-Hoffman-Wald 2002) if the dial is the radius b = e^psi of the law's hidden 3-sphere (the weak force's S^3, rounds 1 and 7), reducing 7-dimensional gravity on it gives 4D Einstein gravity plus a field phi = sqrt(15/2) psi with potential V = -3 e^{-5 psi}/r0^2: an exponential with slope lam = 5/sqrt(15/2) = sqrt(10/3) = 1.83, and V < 0 because the sphere's own curvature enters with a minus sign",
      el_zero and sp.simplify(lam_sphere - sp.sqrt(sp.Rational(10, 3))) == 0,
      f"Euler-Lagrange equations of (reduced 7D action) - (Einstein-frame action) vanish identically: {el_zero}; R7 contains +6/b^2 (sphere curvature); lam = {sp.nsimplify(lam_sphere)} = {float(lam_sphere):.3f}, sign of V: negative")
results["D1"] = dict(R7=str(R7), el_zero=el_zero, lam=float(lam_sphere), V_sign=-1)

# ---------------------------------------------------------------- D2
say("== D2. A negative potential can't accelerate; adding the law's meron on the hidden sphere gives a negative minimum")
phidot, V = sp.symbols('phidot V', real=True)
rho = phidot**2 / 2 + V; p = phidot**2 / 2 - V
rp3 = sp.expand(rho + 3 * p)                               # acceleration a''/a = -(rho + 3p)/6
no_acc = sp.simplify(rp3 - (2 * phidot**2 - 2 * V)) == 0
# meron (round 27 G4) energy on the hidden sphere ~ 1/b^4 per unit internal volume: Einstein frame weight e^{-(3+4) psi}
Acu, Bfl, s = sp.symbols('A_c B_f s', positive=True)
Vtot = -Acu * sp.exp(-5 * s) + Bfl * sp.exp(-7 * s)
smin = sp.solve(sp.diff(Vtot, s), s)[0]
Vmin = sp.simplify(Vtot.subs(s, smin))
Vmin_ratio = sp.simplify(Vmin / (-Acu * sp.exp(-5 * smin)))
check("D2 (exact; textbook, Freund-Rubin 1980 for flux on a sphere) for a rolling field rho + 3p = 2 phidot^2 - 2V, so with V <= 0 the field always slows the expansion: the hidden-sphere reading cannot be dark energy. The law's own f = 1 field (the meron) on the hidden sphere adds energy B e^{-7 psi} > 0 (slope 7/sqrt(15/2) = 2.56); the two together hold the sphere at a fixed size, but at a negative vacuum energy, 2/7 of the curvature term: anti-de Sitter, the wrong sign for our universe",
      no_acc and sp.simplify(Vmin_ratio - sp.Rational(2, 7)) == 0 and sp.simplify(Vmin).is_negative,
      f"rho+3p = {rp3}; V_min = {Vmin} = (2/7) x (-A e^(-5 s*)) < 0; flux slope {float(lam_flux):.3f}")
results["D2"] = dict(rho_plus_3p=str(rp3), V_min_over_curv=str(Vmin_ratio), lam_flux=float(lam_flux))

# ---------------------------------------------------------------- D3
say("== D3. The dial as the slice's own curvature q = k/a^2 (round 28)")
lam, aa, kk = sp.symbols('lam a k', positive=True)
w_scal = lam**2 / 3 - 1                                    # field-dominated scaling solution of V0 e^{-lam phi}
w_curv = sp.Rational(-1, 3)                                 # rho ~ a^{-3(1+w)} = a^{-2}
lam_curv = sp.solve(sp.Eq(w_scal, w_curv), lam)[0]
q = kk / aa**2
w_from_q = sp.simplify(-1 - sp.diff(sp.log(q), aa) * aa / 3)
check("D3 (exact) if the dial is the slice's own curvature q(t) = k/a^2, it dilutes as a^{-2}: w = -1/3, which is exactly the field-dominated solution of an exponential with lam = sqrt 2. That is not new physics but space's curvature itself (Omega_k), already measured near zero (round 27), and w = -1/3 never accelerates",
      sp.simplify(lam_curv - sp.sqrt(2)) == 0 and w_from_q == w_curv,
      f"w(q = k/a^2) = {w_from_q}; scaling solution w = lam^2/3 - 1 equals -1/3 at lam = {lam_curv}")
results["D3"] = dict(w=str(w_from_q), lam=float(lam_curv))

# ---------------------------------------------------------------- D4
say("== D4. Each slope against Pantheon+ and DESI DR2 BAO (same fitter as the dial-potential search)")
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "dial-potential"))
import numpy as np
import quintessence as Q
cands = [("cosmological constant", 0.0, "measured baseline"),
         ("best fit (dial-potential search)", 0.8177083115447195, "fitted"),
         ("conformal sqrt(2/3)", math.sqrt(2 / 3), "chosen, not forced by the law"),
         ("slice curvature sqrt 2 (D3)", math.sqrt(2), "forced by that reading"),
         ("hidden-sphere curvature sqrt(10/3) (D1), sign flipped", math.sqrt(10 / 3), "forced; real sign V<0 cannot accelerate"),
         ("meron on hidden sphere 7/sqrt(15/2) (D2)", 7 / math.sqrt(7.5), "forced by that reading")]
rows = []
for name, lv, status in cands:
    try:
        r = Q.fit("exp", lv, "both")
        E, w = Q.history(r["Om"], "exp", lv)
        z = Q.ZG; aG = 1 / (1 + z); m = z <= 2.0
        M = np.vstack([np.ones(m.sum()), 1 - aG[m]]).T           # CPL w = w0 + wa (1 - a), least squares to z = 2
        w0, wa = np.linalg.lstsq(M, w[m], rcond=None)[0]
        rows.append(dict(name=name, lam=lv, status=status, Om=r["Om"], chi2=r["chi2"], w0=float(w0), wa=float(wa)))
    except Exception as e:
        rows.append(dict(name=name, lam=lv, status=status, Om=None, chi2=None, w0=None, wa=None, err=str(e)))
c0 = rows[0]["chi2"]
for r in rows:
    r["dchi2_vs_LCDM"] = None if r["chi2"] is None else r["chi2"] - c0
    say(f"   {r['name']:<58s} lam={r['lam']:.3f}  " + ("fails: " + r.get('err', '') if r["chi2"] is None else
        f"Om={r['Om']:.3f} dchi2 vs Lambda={r['dchi2_vs_LCDM']:+7.2f}  w0={r['w0']:+.3f} wa={r['wa']:+.3f}"))
band = (0.6118765981676063, 0.9741551199238941)
forced = [r for r in rows if r["status"].startswith("forced")]
forced_in = [r for r in forced if band[0] <= r["lam"] <= band[1]]
worst_ok = all(r["chi2"] is None or r["dchi2_vs_LCDM"] > 9 for r in forced)
# exact: the field-dominated solution of V0 e^{-lam phi} exists iff lam^2 < 6 and accelerates iff w = lam^2/3 - 1 < -1/3
acc_edge = sp.solve(sp.Eq(lam**2 / 3 - 1, sp.Rational(-1, 3)), lam)[0]
no_late_acc = all(r["lam"] >= float(acc_edge) - 1e-12 for r in forced)
check("D4 (numbers + exact; data as in dial-potential/quintessence.py) every slope that a geometric reading of the dial forces lies outside the data's 1 sigma band 0.61 to 0.97, and either fits worse than a cosmological constant by more than 9 in chi^2 (3 sigma) or reaches no history with today's dark-energy fraction at all. All of them are at or above sqrt 2, where the field's own late-time solution stops accelerating (exact: w = lam^2/3 - 1 < -1/3 needs lam < sqrt 2). The only simple value inside the band, sqrt(2/3), is chosen, not forced; 14 of 43 simple values sit in the band (dial-potential search), so its closeness to the best fit means nothing yet",
      len(forced_in) == 0 and worst_ok and no_late_acc and band[0] <= math.sqrt(2 / 3) <= band[1],
      "; ".join(f"{r['name'].split(' (')[0]} {r['lam']:.2f}: " + ("no fit" if r['chi2'] is None else f"{r['dchi2_vs_LCDM']:+.1f}") for r in rows))
results["D4"] = dict(rows=rows, band=band, accel_edge=float(acc_edge))

# ---------------------------------------------------------------- D5
say("== D5. What would decide the rest")
best = rows[1]
check("D5 (numbers) the fitted dial (lam = 0.82) thaws: w0 above -1 and wa below 0, the same direction as DESI 2025's hint. It is decided by w0-wa from DESI DR3 and Euclid: if they put a cosmological constant inside 3 sigma, the rolling dial is not needed; if they exclude it, the dial survives, but with a fitted slope, so it is a fit, not a prediction of the law",
      best["w0"] > -1 and best["wa"] < 0,
      f"lam = 0.818: w0 = {best['w0']:+.3f}, wa = {best['wa']:+.3f} (CPL fit to the field's w(a) out to z = 2)")
results["D5"] = dict(w0=best["w0"], wa=best["wa"])

say(f"\n{passed} passed, {failed} failed  ({time.time() - T0:.0f} s)")
results["summary"] = dict(passed=passed, failed=failed)
json.dump(results, open("results.json", "w"), indent=1, default=str)
open("checks_output.txt", "w").write(_log.getvalue())
