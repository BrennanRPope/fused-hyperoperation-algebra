"""Checks: the static (v = 0) second-order two-body terms, computed independently (fokker2.py: static field equations,
Fokker action, dimensional regularisation) and compared with Gupta's printed Lagrangian (arXiv:2609.07413, Eq. 7.1).
Writes checks_static_output.txt.  The 3-metric exchange is not computed; its shape is fixed (n_ab n_cd with
n_ab = m_a m_b (1 + alpha_a.alpha_b) - q_a.q_b), so differences of that shape are allowed and anything else is not."""
import sys, os; HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import sympy as sp
import pn2, multi
from fokker2 import lagrangian, at3, emd_body
R = sp.Rational; lines = []; ok_all = True
def rep(name, ok, detail=""):
    global ok_all; ok_all &= bool(ok); lines.append(f"{'PASS' if ok else 'FAIL'}  {name}  {detail}"); print(lines[-1], flush=True)
m1, m2, q1, q2, a = sp.symbols('m1 m2 q1 q2 a')
al1, al2, be1, be2, bp1, bp2 = sp.symbols('al1 al2 be1 be2 bp1 bp2')
params = [m1, m2, q1, q2, a, al1, al2, be1, be2, bp1, bp2]
SW = {m1: m2, m2: m1, q1: q2, q2: q1, al1: al2, al2: al1, be1: be2, be2: be1, bp1: bp2, bp2: bp1}
L2, L3, L4 = lagrangian(1, [[a]], {1: emd_body(m1, q1, al1, be1, bp1), 2: emd_body(m2, q2, al2, be2, bp2)})
l2, l3, l4 = (at3(x, params) for x in (L2, L3, L4))
f = 1 + al1*al2
def paper(corrected):
    pn2.CORRECTED = corrected; h1, h2 = pn2.halves(); pn2.CORRECTED = True
    z = {v: 0 for v in pn2.K.values()}
    sub = {pn2.G: 1, pn2.r: 1, pn2.G12: f, pn2.gb: -2*al1*al2/f, pn2.bb1: be1*al2**2/(2*f**2), pn2.bb2: be2*al1**2/(2*f**2),
           pn2.d2: al2**2/f**2, pn2.e1: bp1*al2**3/f**3, pn2.z12: be1*be2*al1*al2/f**3}
    sym = {pn2.m1: m1, pn2.m2: m2, pn2.q1: q1, pn2.q2: q2, pn2.a: a, pn2.al1: al1, pn2.al2: al2, pn2.be1: be1}
    out = []
    for h in (h1, h2):
        h = h.subs(z).subs(sub).subs(sym, simultaneous=True)
        out.append(sp.expand(sp.cancel(h + h.subs(SW, simultaneous=True))))
    return out
p1, p2 = paper(False); c1, c2 = paper(True)
rep("S1 Newton + Coulomb + scalar at zeroth order", sp.expand(l2 - (f*m1*m2 - q1*q2)) == 0, str(l2))
rep("S2 the independent calculation reproduces every first-order static term of the paper", sp.expand(sp.cancel(l3 - p1)) == 0)
lam = sp.Symbol('lam'); D = sp.expand(sp.cancel(l4 - p2)); Dl = sp.expand(D.subs({m2: lam*m2, q2: lam*q2}, simultaneous=True))
n = lambda A, B: {(1, 1): (1 + al1**2)*m1**2 - q1**2, (2, 2): (1 + al2**2)*m2**2 - q2**2, (1, 2): f*m1*m2 - q1*q2}[(A, B)]
bid = {k: sp.expand(Dl.coeff(lam, k)) for k in range(5)}
ok3 = bid[0] == 0 and bid[4] == 0 and sp.expand(bid[1] + n(1, 1)*n(1, 2)/3) == 0 and sp.expand(bid[3] + n(2, 2)*n(1, 2)/3) == 0
rep("S3 second-order static terms with three powers of one body (and four) agree with the paper up to the 3-metric shape n11 n12", ok3,
    "difference (3,1) = -n11 n12/3")
extra = sp.expand(bid[2] + 2*n(1, 2)**2)
ok4 = sp.expand(extra - 12*(1 - a**2*al1*al2)*m1*m2*q1*q2) == 0
rep("S4 in the both-bodies-twice terms the only difference beyond the 3-metric shape is 12 (1 - a^2 al1 al2) G^2 m1 m2 q1 q2 / r^3 (not of n12^2 or n11 n22 form)",
    ok4, "paper - ours = 2 n12^2 - " + str(sp.factor(extra)))
rep("S5 with that one term corrected, the paper's static terms equal ours up to the 3-metric shape", sp.expand(sp.cancel(c2 - l4 - 2*n(1, 2)**2 - (n(1, 1) + n(2, 2))*n(1, 2)/3)) == 0)
bps = {al1: a, al2: a, be1: 0, be2: 0, bp1: 0, bp2: 0, q1: sp.sqrt(1 + a**2)*m1, q2: sp.sqrt(1 + a**2)*m2}
rep("S6 printed paper: two extremal same-sign black holes feel a static second-order force for every coupling except a = 1",
    sp.factor(sp.expand(p2.subs(bps))) == sp.factor(-12*(1 - a**2)*(1 + a**2)**2*m1**2*m2**2), "static energy term " + str(sp.factor(sp.expand(p2.subs(bps)))))
rep("S7 ours (and the corrected paper): no static force between extremal same-sign black holes, every coupling, first and second order",
    sp.expand(l3.subs(bps)) == 0 and sp.expand(l4.subs(bps)) == 0 and sp.expand(c2.subs(bps)) == 0)
# multi-field (STU): four photons, three scalars
CV = multi.CV
def stu_V(Bd):
    nf = 8
    def entry(idx):
        chi = [i for i in idx if i >= 4]
        if chi: return -Bd.Q[chi[0] - 4] if len(idx) == 1 else 0
        ph = [i - 1 for i in idx if 1 <= i <= 3]
        return Bd.m*sum(wj*sp.prod([cj[s] for s in ph]) for wj, cj in zip(Bd.w, CV))
    T1 = [entry((i,)) for i in range(nf)]; T2 = [[entry((i, j)) for j in range(nf)] for i in range(nf)]
    T3 = [[[entry((i, j, k)) for k in range(nf)] for j in range(nf)] for i in range(nf)]
    return (T1, T2, T3, None)
def nn(A, B): return A.m*B.m*(1 + multi.dot(A.al, B.al)) - sum(x*y for x, y in zip(A.Q, B.Q))
X = sp.Symbol('X'); okB = True; okL = True; det = []
pairs = [((1, 0, 0, 0), (0, 1, 0, 0)), ((1, 0, 0, 0), (0, -1, 0, 0)), ((1, 2, 3, 0), (2, 0, 1, 1)), ((1, 1, 0, 0), (1, 1, 0, 0)), ((0, 0, 0, 2), (-1, -1, -1, 0)),
         ((1, 0, 0, 0), (-1, 0, 0, 0)), ((1, 2, 0, 0), (0, -1, -3, 0)), ((1, 1, 1, 1), (-2, 0, -1, 0)), ((5, 1, 2, 3), (-1, -4, -2, -1))]
for QA, QB in pairs:
    A = multi.Body(QA); B = multi.Body(QB)
    L2m, L3m, L4m = lagrangian(3, [list(c) for c in CV], {1: stu_V(A), 2: stu_V(B)})
    s4 = sp.nsimplify(at3(L4m, [X]) + 2*nn(A, B)**2 + (nn(A, A) + nn(B, B))*nn(A, B)/3)
    s3 = at3(L3m, [X])
    Lf1, Lf2 = multi.full(A, B); sub = {pn2.G: 1, pn2.r: 1, pn2.m1: A.m, pn2.m2: B.m}
    t3 = sp.expand(multi.static(Lf1).subs(sub)); t4 = sp.expand(multi.static(Lf2).subs(sub).subs({multi.g1s: 0, multi.s1s: 0}))
    if nn(A, B) == 0: okB &= (s3 == 0 and s4 == 0)
    okL &= sp.simplify(s3 - t3) == 0 and sp.simplify(s4 - t4) == 0
    det.append(f"{QA}|{QB}: n12={nn(A, B)}, static2={s4}")
rep("S8 four-charge (STU) bodies, charges pointing different ways: no static force at first or second order whenever n12 = 0 (independent calculation)", okB, "; ".join(det))
rep("S9 the lifted multi-field equations (multi.py, corrected term, g1 = s1) equal the independent static calculation for all 9 pairs", okL)
lines.append("ALL PASS" if ok_all else "SOME FAILED"); print(lines[-1])
open(os.path.join(HERE, "checks_static_output.txt"), "w").write("\n".join(lines) + "\n")
