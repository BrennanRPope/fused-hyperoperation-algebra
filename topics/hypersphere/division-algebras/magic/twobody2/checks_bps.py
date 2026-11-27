"""Check: the second-order (2PN) velocity terms of Gupta's Lagrangian, at order G^2 v^2, for two slow extremal
same-sign dilaton black holes (coupling a; no static force) against the known moduli-space metric (Shiraishi,
Nucl. Phys. B402 (1993) 399, arXiv:1407.5377, Eqs. 3, 7, 18): relative motion is free motion in the metric
g(r) (dr^2 + r^2 dphi^2), g = 1 + (3 - a^2)/2 * Int F^{2(1-a^2)/(1+a^2)} ..., F = 1 + (1 + a^2) G M / r, which gives
g = 1 + c1 x + c2 x^2 + ..., x = G M/r, c1 = 3 - a^2, c2 = (3 - a^2)(1 - a^2) for any mass ratio.
From our Hamiltonian H = p^2/(2 mu) [1 + A1 x + A2 x^2] + p_r^2/(2 mu) B2 x^2: coordinates are free up to
r -> r + lambda G^2 M^2/r, which leaves 2 A2 + B2 unchanged; the metric predicts A1 = -c1, 2 A2 + B2 = 2 (c1^2 - c2).
Writes checks_bps_output.txt."""
import sys, os; HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import sympy as sp, pn2
R = sp.Rational; Xs = sp.Symbol('X', positive=True); lines = []; ok_all = True
G, m1, m2, q1, q2, a = pn2.G, pn2.m1, pn2.m2, pn2.q1, pn2.q2, pn2.a
L1, L2 = pn2.full_L(); px, py, r = sp.Symbol('px'), sp.Symbol('py'), pn2.r
for av in [0, R(1, 2), R(3, 4), R(4, 3), sp.sqrt(3), sp.sqrt(2)]:
    f = 1 + av**2
    sub = {a: av, pn2.al1: av, pn2.al2: av, pn2.be1: 0, pn2.be2: 0, pn2.bb1: 0, pn2.bb2: 0, pn2.e1: 0, pn2.e2: 0, pn2.z12: 0,
           pn2.d1: av**2/f**2, pn2.d2: av**2/f**2, pn2.gb: -2*av**2/f, pn2.G12: G*f, q1: sp.sqrt(f*G)*m1, q2: sp.sqrt(f*G)*m2}
    sub = {k: (v.subs({m1: Xs, m2: 1 - Xs}) if hasattr(v, 'subs') else v) for k, v in sub.items()}; sub.update({m1: Xs, m2: 1 - Xs})
    H0, H1, H2, k = pn2.hamiltonian(L1, L2, sub); mu = Xs*(1 - Xs)
    def coef(H, n):
        P = sp.Poly(sp.expand(sp.expand(H).coeff(G, n)), px, py)
        return {mon: sp.simplify(c*2*mu*r**n) for mon, c in P.terms() if sum(mon) == 2}
    c1d, c2d = coef(H1, 1), coef(H2, 2)
    A1x, A1y = c1d.get((2, 0), 0), c1d.get((0, 2), 0); A2 = c2d.get((0, 2), 0); B2 = sp.simplify(c2d.get((2, 0), 0) - A2)
    c1 = 3 - av**2; c2 = (3 - av**2)*(1 - av**2)
    ok = sp.simplify(k) == 0 and sp.simplify(A1x + c1) == 0 and sp.simplify(A1y + c1) == 0 and sp.simplify(2*A2 + B2 - 2*(c1**2 - c2)) == 0
    ok_all &= bool(ok)
    lines.append(f"{'PASS' if ok else 'FAIL'}  a = {av}: no static pull ({k}); order G: {A1x}, {A1y} (metric {-c1}); "
                 f"order G^2 invariant 2A2 + B2 = {sp.simplify(2*A2 + B2)} (metric {2*(c1**2 - c2)}), same at every mass ratio")
    print(lines[-1], flush=True)
lines.append("ALL PASS" if ok_all else "SOME FAILED"); print(lines[-1])
open(os.path.join(HERE, 'checks_bps_output.txt'), 'w').write("\n".join(lines) + "\n")
