"""DRAFT round 28: does anything in the law force matter's rest slices to be its 3-spheres (closed space)?
Exact, symbolic. Checks that one de Sitter spacetime (one Lambda = 3q0) carries closed, flat and open slicings,
that the law's dial on the slices is q(t) = k/a(t)^2 with any sign, and that with matter k is free initial data. Also tests the audit C4 row: the law torsion vs Einstein-Cartan."""
import json, io, time
import sympy as sp

T0 = time.time(); _log = io.StringIO(); passed = failed = 0; results = {}
def say(s): print(s); _log.write(s + "\n")
def check(name, ok, detail):
    global passed, failed
    passed += bool(ok); failed += (not ok)
    say(f"{'PASS' if ok else 'FAIL'} {name}  | {detail}")

t, H, k, rho0, Lam, a0 = sp.symbols('t H k rho0 Lambda a0', positive=True)
kk = sp.symbols('kk')

# ================= S1. one de Sitter, three slicings
say("== S1. One de Sitter spacetime, three ways to slice it")
slicings = {+1: sp.cosh(H * t) / H, 0: sp.exp(H * t), -1: sp.sinh(H * t) / H}
fr = {}
for kv, a in slicings.items():
    res = sp.simplify((sp.diff(a, t) / a)**2 + kv / a**2 - H**2)        # Friedmann with Lambda = 3H^2, no matter
    acc = sp.simplify(sp.diff(a, t, 2) / a - H**2)
    fr[kv] = (res, acc)
ok1 = all(r == 0 and a_ == 0 for r, a_ in fr.values())
# the three are the same spacetime: each is a slicing of the hyperboloid -X0^2 + X1^2 + ... + X4^2 = 1/H^2
X = sp.symbols('X0:5'); chi, r = sp.symbols('chi r', positive=True)
emb = {
 +1: [sp.sinh(H * t) / H, sp.cosh(H * t) / H * sp.cos(chi), sp.cosh(H * t) / H * sp.sin(chi)],
 0: [sp.sinh(H * t) / H + H * r**2 * sp.exp(H * t) / 2, sp.cosh(H * t) / H - H * r**2 * sp.exp(H * t) / 2, sp.exp(H * t) * r],
 -1: [sp.sinh(H * t) / H * sp.cosh(chi), sp.sinh(H * t) / H * sp.sinh(chi), sp.cosh(H * t) / H],
}
# 2D reduction (dS_2 in R^{2,1}): -X0^2 + X1^2 + X2^2 = 1/H^2 on a one-dimensional slice coordinate
on_hyp = {kv: sp.simplify(sp.expand(-e[0]**2 + e[1]**2 + e[2]**2 - 1 / H**2).rewrite(sp.exp)) for kv, e in emb.items()}
ok1b = all(v == 0 for v in on_hyp.values())
check("S1 (exact; textbook, de Sitter 1917, Lemaitre 1925, Schrodinger 1956) one de Sitter spacetime with Lambda = 3H^2 solves Friedmann's equation with no matter as a closed universe a = cosh(Ht)/H (k = +1, the law's slices), a flat one a = e^(Ht) (k = 0) and an open one a = sinh(Ht)/H (k = -1). All three lie on the same hyperboloid: they are one spacetime cut three ways, picked by a timelike, null or spacelike direction in the embedding",
      ok1 and ok1b, f"Friedmann and acceleration residuals {[(kv, fr[kv]) for kv in fr]}; on the hyperboloid: {on_hyp}")
results["S1"] = dict(friedmann={str(kv): [str(x) for x in v] for kv, v in fr.items()}, on_hyperboloid={str(kv): str(v) for kv, v in on_hyp.items()})

# ================= S2. the law's dial on each slicing: q(t) = k/a(t)^2, any sign, changing in time
say("== S2. The law's dial on the slices")
qs = {kv: sp.simplify(kv / a**2) for kv, a in slicings.items()}
q_at_rest = sp.simplify(qs[+1].subs(t, 0))
ok2 = q_at_rest == H**2 and qs[0] == 0 and all(sp.simplify(sp.diff(qs[kv], t)) != 0 for kv in (+1, -1))
check("S2 (exact) on each slicing the slices are space of constant curvature, so the law's dial there is q(t) = k/a(t)^2: positive (the law's turning 3-spheres, round 18), zero (plain vector addition) or negative (velocity-like boosts, round 19), all with the same Lambda. Only the closed slicing at its moment of rest has q = Lambda/3 (round 20); everywhere else the slice's dial drifts with time and is not Lambda/3",
      ok2, f"q(t): closed {qs[1]}, flat {qs[0]}, open {qs[-1]}; closed at t = 0: {q_at_rest} = Lambda/3")
results["S2"] = {str(kv): str(v) for kv, v in qs.items()}

# ================= S3. with matter, k is free initial data, independent of Lambda's sign
say("== S3. With matter, the sign of space's dial is free")
a = sp.Function('a')(t)
# Friedmann: (a'/a)^2 = 8 pi G rho/3 - k/a^2 + Lambda/3 with rho = rho0 (a0/a)^3; today's data (a0, H0, rho0) fix k:
H0s, Gs = sp.symbols('H0 G', positive=True)
k_today = sp.solve(sp.Eq(H0s**2, 8 * sp.pi * Gs * rho0 / 3 - kk / a0**2 + Lam / 3), kk)[0]
signs = {}
for lab_, vals in (("closed", {H0s: 1, Gs: 1, rho0: sp.Rational(3, 8) / sp.pi * 2, Lam: 3, a0: 1}),
                   ("flat", {H0s: 1, Gs: 1, rho0: sp.Rational(3, 8) / sp.pi * sp.Rational(1, 2), Lam: sp.Rational(3, 2), a0: 1}),
                   ("open", {H0s: 1, Gs: 1, rho0: sp.Rational(3, 8) / sp.pi * sp.Rational(1, 4), Lam: sp.Rational(3, 2), a0: 1})):
    signs[lab_] = sp.nsimplify(k_today.subs(vals))
ok3 = signs["closed"] > 0 and signs["flat"] == 0 and signs["open"] < 0
check("S3 (exact; textbook Friedmann) with matter, k = a0^2 (8 pi G rho0/3 + Lambda/3 - H0^2) is set by today's density and expansion rate. With the same positive Lambda it can be +, 0 or - (three examples here): nothing ties space's sign to Lambda's. The law's step 'q > 0 closes, so Lambda > 0' (round 20, audit B7) ties them only at a moment of rest, which round 27 showed our universe never had",
      ok3, f"k = {k_today}; examples with Lambda > 0: closed {signs['closed']}, flat {signs['flat']}, open {signs['open']}")
results["S3"] = dict(k_today=str(k_today), examples={kx: str(v) for kx, v in signs.items()})

# ================= S4. audit C4: does the law fix a spin-twist (torsion) coupling different from Einstein-Cartan's?
say("== S4. The law's torsion and Einstein-Cartan (audit C4)")
f_ = sp.symbols('f')
eps = lambda i, j, l: sp.LeviCivita(i, j, l)
# left-invariant frame on the unit 3-sphere: [X_a, X_b] = 2 eps_abc X_c; family nabla_X Y = f [X, Y]/2
br = lambda a_, b_: [2 * eps(a_, b_, c_) for c_ in range(3)]
Tcomp = [[[sp.simplify(f_ * br(a_, b_)[c_] / 2 - f_ * br(b_, a_)[c_] / 2 - br(a_, b_)[c_]) for c_ in range(3)] for b_ in range(3)] for a_ in range(3)]
T_law = [[[Tcomp[a_][b_][c_].subs(f_, 1) for c_ in range(3)] for b_ in range(3)] for a_ in range(3)]
T012 = sp.factor(Tcomp[0][1][2])
law_free = all(T_law[a_][b_][c_] == 0 for a_ in range(3) for b_ in range(3) for c_ in range(3))
# the flat twisted alternatives f = 0, 2 on today's cosmic slices: torsion size 2/a0 per unit length, in energy units hbar c * 2/a0
hbarc_GeV_m = 1.97327e-16
a0_m = 2.51e27                      # round 27 G4: smallest closed curvature radius allowed (Planck + BAO, 2 sigma)
T_cosmic_GeV = hbarc_GeV_m * 2 / a0_m
bound_GeV = 1e-31                   # spin-polarized torsion bounds, Kostelecky, Russell and Tasson 2008
ok4 = law_free and sp.simplify(T012 - 2 * (f_ - 1)) == 0 and T_cosmic_GeV < bound_GeV
check("S4 (exact + numbers; audit C4) the law's own way of carrying directions is the f = 1 member of the family nabla_X Y = f[X,Y]/2, whose torsion is 2(f - 1): zero at the law. So the law's 'twist' gyr is curvature (Levi-Civita holonomy, round 25), not torsion, and the law has no spin-torsion term to fix. It agrees with both GR and Einstein-Cartan outside matter (where Einstein-Cartan torsion vanishes) and says nothing inside spinning matter, where Einstein-Cartan's coupling 8 pi G is set by its action, not by the law. Even the flat twisted members f = 0, 2 on today's cosmic slices would give torsion of only 2 hbar c/a0, about twelve orders below the bounds",
      ok4, f"T^2_01 = {T012}; law (f = 1) torsion all zero: {law_free}; f = 0 or 2 on a0 >= {a0_m:.2e} m: {T_cosmic_GeV:.1e} GeV vs bound ~{bound_GeV:.0e} GeV")
results["S4"] = dict(torsion_component=str(T012), law_torsion_zero=law_free, cosmic_flat_twist_GeV=T_cosmic_GeV, bound_GeV=bound_GeV)

say(f"\n{passed} passed, {failed} failed  ({time.time() - T0:.0f} s)")
results["summary"] = dict(passed=passed, failed=failed)
json.dump(results, open("results.json", "w"), indent=1, default=str)
open("checks_output.txt", "w").write(_log.getvalue())
