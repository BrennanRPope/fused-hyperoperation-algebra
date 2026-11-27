"""Two comparable extremal four-charge black holes (STU), each charge carrying its own sign.
Probe limit 'B around A': exact background of A, B as point particle made of its charge pieces.
Two-body turn at first order: T = X_A T(B around A) + X_B T(A around B)   (pn1.py: linear in X, no nu)."""
import sys; sys.path.insert(0, '..')
import sympy as sp
from misalign import W_stu
from precess import T_series
def probe(kA, sA, kB, sB):
    mB = sum(kB); ws = [sp.S(k)/mB for k in kB]; rel = [a*b for a, b in zip(sA, sB)]
    W = W_stu([sp.S(k) for k in kA], ws, rel)
    return T_series(W)
def pair_rule(kA, sA, kB, sB):
    mB = sum(kB); w = [sp.S(k)/mB for k in kB]; x = [w[j]*(1 - sA[j]*sB[j]) for j in range(4)]
    K = sum(kA[j]*x[j] for j in range(4))
    N = 2*sum(kA[i]*kA[j]*(x[i]*(1 - w[j]) + x[j]*(1 - w[i]) + x[i]*x[j]) for i in range(4) for j in range(i + 1, 4))
    return sp.simplify(N/K**2), K
def two_body(kA, sA, kB, sB, exact=True):
    f = probe if exact else pair_rule
    TA, KA = f(kA, sA, kB, sB); TB, KB = f(kB, sB, kA, sA)
    mA, mB = sum(kA), sum(kB); XA = sp.S(mA)/(mA + mB)
    sym = sp.simplify(KA/mA - KB/mB)          # Newton's third law: pull per source mass is shared
    return sp.nsimplify(sp.simplify(XA*TA + (1 - XA)*TB)), TA, TB, sym
if __name__ == "__main__":
    P, M = 1, -1
    cases = {
        "1 charge vs opposite 1 charge":            ([2, 0, 0, 0], [P]*4, [2, 0, 0, 0], [M]*4),
        "2 charges (+,+) vs (-,-), equal":          ([1, 1, 0, 0], [P]*4, [1, 1, 0, 0], [M]*4),
        "2 charges (+,+) vs (-,-), unequal masses": ([3, 3, 0, 0], [P]*4, [1, 1, 0, 0], [M]*4),
        "2 charges (+,+) vs (-,+)":                 ([1, 1, 0, 0], [P]*4, [1, 1, 0, 0], [M, P, P, P]),
        "2 charges (1,2) vs 1 charge opposite":     ([1, 2, 0, 0], [P]*4, [3, 0, 0, 0], [M]*4),
        "4 charges vs 4 opposite, equal":           ([1, 1, 1, 1], [P]*4, [1, 1, 1, 1], [M]*4),
        "4 charges vs 1 opposite, equal mass":      ([1, 1, 1, 1], [P]*4, [4, 0, 0, 0], [M]*4),
        "3 charges vs 2 charges, mixed signs":      ([1, 2, 3, 0], [P, M, P, P], [0, 2, 1, 3], [P, P, M, M]),
    }
    for name, (kA, sA, kB, sB) in cases.items():
        T, TA, TB, sym = two_body(kA, sA, kB, sB)
        Tp = two_body(kA, sA, kB, sB, exact=False)[0]
        print(f"{name:45s} T={str(T):8s} (B round A {TA}, A round B {TB}; pair rule {Tp}; 3rd law {sym})")
