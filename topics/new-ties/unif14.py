"""Round 14 (DRAFT): excursions off the line, for every stretch at once.

A nonlinear excursion starts at a linear neighbour map T(y) = alpha y (see nb14.py) and continues with words P (h side) and
Q (g side): T = Q^{-1} o alpha o P, so a(alpha) = alpha s_P/s_Q and b(alpha) = (alpha t_P - t_Q)/s_Q are linear in alpha.
Instead of fixing alpha, keep the closed interval I of alpha in [1/20, 20] for which every step so far passed the real
cut (T(J) meets J); the peel rule (h if a >= 1) splits I where needed. Grammar before the excursion is over-approximated
(last letters l or m on each side), which can only add branches.
The excursion comes back to the line (b = 0) exactly at alpha* = t_Q/t_P, if that lies in I. Such an alpha* must be a
{2,5,19}-unit, since a stretch is one. This script counts the search nodes by length and lists the alpha* values that are
units. Usage: python3 unif14.py [maxlen=11] (about 15 s; each extra letter costs about 3.5 times more). Writes unif14.json.
"""
import sys, json
from fractions import Fraction as F
from collections import Counter
L = {'l': (F(19, 200), F(0)), 'm': (F(1, 20), F(0)), 'N': (F(1, 19), F(1, 19)), 'K': (F(1, 10), F(1, 10))}
NXT = {'l': 'lmNK', 'm': 'mN', 'N': 'lmNK', 'K': 'lmNK'}
J = F(1, 9)
def cut(lo, hi, sP, tP, sQ, tQ):
    if tP > 0: hi = min(hi, (tQ + sQ * J) / tP)
    elif tQ + sQ * J < 0: return None
    lo = max(lo, tQ / (tP + sP * J))
    return (lo, hi) if lo <= hi else None
def sunit(q):
    for n in (q.numerator, q.denominator):
        for p in (2, 5, 19):
            while n % p == 0: n //= p
        if n != 1: return False
    return True
def v(q, p):
    k = 0; a, b = q.numerator, q.denominator
    while a % p == 0: a //= p; k += 1
    while b % p == 0: b //= p; k -= 1
    return k
def run(maxlen):
    stack = []
    for lh in 'lm':
        for lg in 'lm':
            for x in NXT[lh]:
                if x in 'NK':
                    sP, tP = L[x]; c = cut(F(1), F(20), sP, tP, F(1), F(0))
                    if c: stack.append((c[0], c[1], sP, tP, F(1), F(0), x, lg, x, ''))
            for y in NXT[lg]:
                if y in 'NK':
                    s, t = L[y]; c = cut(F(1, 20), F(1), F(1), F(0), s, t)
                    if c: stack.append((c[0], c[1], F(1), F(0), s, t, lh, y, '', y))
    nodes = Counter(); ret_all = Counter(); ret_unit = Counter(); alphas = Counter(); first = {}
    while stack:
        lo, hi, sP, tP, sQ, tQ, lh, lg, wP, wQ = stack.pop()
        dep = len(wP) + len(wQ); nodes[dep] += 1
        if tP != 0 and dep > 1:
            az = tQ / tP
            if lo <= az <= hi:
                ret_all[dep] += 1
                if sunit(az):
                    ret_unit[dep] += 1; alphas[az] += 1; first[az] = min(first.get(az, 99), dep)
        if dep >= maxlen: continue
        c0 = sQ / sP
        if hi >= c0:
            for x in NXT[lh]:
                s, t = L[x]; nsP, ntP = sP * s, sP * t + tP
                r = cut(max(lo, c0), hi, nsP, ntP, sQ, tQ)
                if r: stack.append((r[0], r[1], nsP, ntP, sQ, tQ, x, lg, wP + x, wQ))
        if lo < c0:
            for y in NXT[lg]:
                s, t = L[y]; nsQ, ntQ = sQ * s, sQ * t + tQ
                r = cut(lo, min(hi, c0), sP, tP, nsQ, ntQ)
                if r: stack.append((r[0], r[1], sP, tP, nsQ, ntQ, lh, y, wP, wQ + y))
    return nodes, ret_all, ret_unit, alphas, first
if __name__ == "__main__":
    maxlen = int(sys.argv[1]) if len(sys.argv) > 1 else 11
    nodes, ret_all, ret_unit, alphas, first = run(maxlen)
    out = dict(maxlen=maxlen, nodes_by_length=dict(sorted(nodes.items())), returns_by_length=dict(sorted(ret_all.items())),
               unit_returns_by_length=dict(sorted(ret_unit.items())),
               unit_alphas={str(a): dict(count=c, v2=v(a, 2), v5=v(a, 5), v19=v(a, 19), first_length=first[a]) for a, c in alphas.items()})
    json.dump(out, open("unif14.json", "w"), indent=1)
    print("nodes by length", dict(sorted(nodes.items())))
    print("returns to the line (any alpha*)", sum(ret_all.values()), " with alpha* a {2,5,19}-unit", sum(ret_unit.values()))
    print("distinct unit alpha*:", sorted(alphas, key=float), " 19-adic valuations", sorted(set(v(a, 19) for a in alphas)))
