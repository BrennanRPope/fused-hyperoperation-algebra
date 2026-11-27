"""Second-order (2PN) two-body Lagrangian for black holes whose charges point different ways: the STU model
(four photons, three flat scalars, photon couplings e^{-2 c_i.phi} F_i^2 with c_i.c_j = 4 delta_ij - 1).
Built by lifting Gupta's one-photon, one-scalar Einstein-Maxwell-dilaton halves (arXiv:2609.07413, Eqs. 5.9, 7.1;
gupta_halves.txt) term by term: every scalar line ends on a body (alpha, beta, third cumulant) or on a photon line
(c_i), every photon line keeps its index from charge to charge, so each single-field monomial has one multi-field
image -- except the G q1^2 q2^2 terms, where (q1.q2)^2 vs |q1|^2|q2|^2 (and the c_i.c_j pair of them) is a free
split g1 + g2 = 3, s1 + s2 = 1. That split is fixed by demanding no static force between any two BPS bodies.
Bodies: pieces of extremal charge in the four slots, all of one sign (each body BPS on its own).
  m = sum |Q_i| / 2 (G = 1),  w_i = |Q_i|/sum|Q|,  alpha = sum w_i c_i,  beta = sum w_i c_i c_i^T - alpha alpha^T,
  kappa = third cumulant of c over the pieces.   A single piece has alpha.alpha = 3 (Kaluza-Klein, a^2 = 3)."""
import os, sys, sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import pn2
from pn2 import K, G, m1, m2, r, G12, gb, bb1, bb2, d2, e1, z12, L1s, L2s
R = sp.Rational
CV = [sp.Matrix(v) for v in [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]]
g1s, s1s = sp.symbols('g1 s1')

class Body:
    def __init__(self, Q, c=CV):
        Q = [sp.nsimplify(x) for x in Q]; sg = {sp.sign(x) for x in Q if x != 0}
        assert len(sg) == 1, "each body must be BPS on its own (one sign)"
        self.Q, self.c = Q, c; tot = sum(abs(x) for x in Q); self.m = tot/2
        self.w = [abs(x)/tot for x in Q]
        dim = c[0].shape[0]
        self.al = sum((wi*ci for wi, ci in zip(self.w, c)), sp.zeros(dim, 1))
        self.be = sum((wi*ci*ci.T for wi, ci in zip(self.w, c)), sp.zeros(dim, dim)) - self.al*self.al.T
    def kappa(self, x):  # third cumulant contracted three times with x
        return sum(wi*((ci - self.al).T*x)[0]**3 for wi, ci in zip(self.w, self.c))
def dot(x, y): return (x.T*y)[0]

def static_and_G(A, B, gg1=g1s, ss1=s1s, Gs=G):
    """The G-proportional (charge-scalar) part of the half ell(A, B): Gupta's terms with a or q, lifted.
    A plays body 1 (mass m1), B body 2 (mass m2)."""
    c = A.c; n = len(c); q1, q2 = A.Q, B.Q
    ca1 = [dot(ci, A.al) for ci in c]; ca2 = [dot(ci, B.al) for ci in c]; cb1a2 = [dot(ci, A.be*B.al) for ci in c]
    s12 = dot(A.al, B.al); a11 = dot(A.al, A.al); cc = [[dot(ci, cj) for cj in c] for ci in c]
    I = range(n); v = K
    l1 = Gs*m1*sum((1 + ca1[i])*(q1[i]*q2[i] - q2[i]**2/2) for i in I)/r**2
    S1 = Gs**2*m1**2*sum(-R(2, 3)*ca1[i]**2*q1[i]*q2[i] + ca1[i]**2*q2[i]**2 - R(4, 3)*ca1[i]*q1[i]*q2[i] + 2*ca1[i]*q2[i]**2
                         - R(1, 3)*a11*q1[i]*q2[i] - q1[i]*q2[i] + q2[i]**2 for i in I)/r**3
    S2 = Gs**2*m1*m2*sum(R(1, 3)*q1[i]**2*ca1[i]*ca2[i] + 5*q1[i]*q2[i]*ca1[i]*ca2[i] - q1[i]*q2[i]*(ca1[i]*s12 + cb1a2[i])
                         + R(1, 2)*q2[i]**2*(ca1[i]*s12 + cb1a2[i]) + R(1, 3)*q1[i]**2*ca1[i] - 3*q1[i]*q2[i]*ca1[i] + R(5, 6)*q2[i]**2*ca1[i]
                         + R(1, 6)*s12*q1[i]**2 - 3*s12*q1[i]*q2[i] + R(1, 2)*q1[i]**2 - 10*q1[i]*q2[i]
                         + (6*q1[i]*q2[i]*(1 - ca1[i]*ca2[i]) if pn2.CORRECTED else 0) for i in I)/r**3
    S3 = Gs*m1*(R(1, 2)*v['v2v2']*sum(q2[i]**2*ca1[i] for i in I)
                + v['nv1']**2*sum(q1[i]*q2[i]*(1 + ca1[i]) - R(1, 2)*q2[i]**2 for i in I)
                + v['nv1']*v['nv2']*sum(-2*q1[i]*q2[i]*(1 + ca1[i]) + q2[i]**2*(2 + ca1[i]) for i in I)
                - v['nv2']**2*sum(q2[i]**2*(1 + R(1, 2)*ca1[i]) for i in I)
                - R(1, 2)*v['v1v1']*sum(q1[i]*q2[i]*(1 + ca1[i]) for i in I)
                - R(1, 4)*v['v1v1']*sum(q2[i]**2*(1 - ca1[i]) for i in I)
                + R(1, 2)*v['v1v2']*sum(q2[i]**2*(1 - ca1[i]) for i in I))/r**2
    nq1 = sum(x**2 for x in q1); nq2 = sum(x**2 for x in q2); qq = sum(x*y for x, y in zip(q1, q2))
    Q22 = (gg1*nq1*nq2 + (3 - gg1)*qq**2 + ss1*sum(q1[i]**2*q2[j]**2*cc[i][j] for i in I for j in I)
           + (1 - ss1)*sum(q1[i]*q2[i]*q1[j]*q2[j]*cc[i][j] for i in I for j in I))
    S4 = Gs*(-R(1, 3)*sum(q1[i]**2*cc[i][j]*q1[j]*q2[j] for i in I for j in I) + Q22/2)/r**3
    return l1, S1 + S2 + S3 + S4

def gravity_part(A, B):
    """Gupta's terms without G (gravity-scalar block in G12, gb, bb, d, e, z; kinetic terms; the q1 q2 kinematic blocks),
    with the single-field combinations replaced by their multi-field contractions."""
    l1, l2 = pn2.halves()
    s12 = dot(A.al, B.al); f = sp.Symbol('f'); qq = sum(x*y for x, y in zip(A.Q, B.Q))   # f = 1 + s12 kept symbolic (can be 0)
    sub = {G12: f, gb: -2*s12/f, bb1: dot(B.al, A.be*B.al)/(2*f**2), bb2: dot(A.al, B.be*A.al)/(2*f**2),
           d2: dot(B.al, B.al)/f**2, e1: A.kappa(B.al)/f**3, z12: dot(A.al, B.be*A.be*B.al)/f**3,
           pn2.q1: 1, pn2.q2: qq}
    lam = {L1s: 7 + 10*gb + R(7, 2)*gb**2 - 2*d2 + 2*bb1, L2s: 8 + 10*gb + R(7, 2)*gb**2 - 2*d2 + 4*bb1}
    out = []
    for ex in (l1, l2):
        ex = sp.expand(ex.subs(lam).subs(G, 0).subs(sub, simultaneous=True))
        assert not any(fa.is_Pow and fa.base == f and fa.exp < 0 for t in sp.Add.make_args(ex) for fa in sp.Mul.make_args(t))
        out.append(sp.expand(ex.subs(f, 1 + s12)))
    return out

SW = {m1: m2, m2: m1, K['v1v1']: K['v2v2'], K['v2v2']: K['v1v1'], K['nv1']: -K['nv2'], K['nv2']: -K['nv1'], K['na1']: -K['na2'],
      K['na2']: -K['na1'], K['v1a2']: K['a1v2'], K['a1v2']: K['v1a2'], K['v1a1']: K['v2a2'], K['v2a2']: K['v1a1']}
def half(A, B, gg1=g1s, ss1=s1s):
    a1, a2 = static_and_G(A, B, gg1, ss1, Gs=1); b1, b2 = gravity_part(A, B)
    return sp.expand(a1 + b1), sp.expand(a2 + b2)
def full(A, B, gg1=g1s, ss1=s1s):
    h1, h2 = half(A, B, gg1, ss1); i1, i2 = half(B, A, gg1, ss1)
    return h1 + i1.subs(SW, simultaneous=True), h2 + i2.subs(SW, simultaneous=True)
def static(ex):
    z = {v: 0 for v in K.values()}; return sp.expand(ex.subs(z))
def pull(A, B):
    return sp.nsimplify((1 + dot(A.al, B.al))*A.m*B.m - sum(x*y for x, y in zip(A.Q, B.Q)))

# The G q1^2 q2^2 split enters only as (g1 - s1)(|q1|^2|q2|^2 - (q1.q2)^2) (c_i.c_j = 4 delta_ij - 1). With the paper's
# m1 m2 q1 q2 term corrected (pn2.CORRECTION) and g1 = s1, the lift equals the independent static calculation
# (fokker2.py) for every pair tested; g1 = s1 = 0 below.
G1, S1 = sp.Integer(0), sp.Integer(0)
def turns(A, B, gg1=None, ss1=None):
    gg1 = G1 if gg1 is None else gg1; ss1 = S1 if ss1 is None else ss1
    L1, L2 = full(A, B, gg1, ss1)
    k = pull(A, B); qq = sum(x*y for x, y in zip(A.Q, B.Q))
    sub = {G: 1, m1: A.m, m2: B.m, G12: (k + qq)/(A.m*B.m), pn2.q1: 1, pn2.q2: qq}
    H0, H1, H2, kk = pn2.hamiltonian(L1, L2, sub)
    assert sp.simplify(kk - k) == 0
    return pn2.turn_coeffs(H0, H1, H2, kk, A.m*B.m/(A.m + B.m))
def body(ws, ss, m):
    """Body of mass m from piece weights ws and signs ss (all signs equal)."""
    return Body([2*m*sv*wv for wv, sv in zip(ws, ss)])
