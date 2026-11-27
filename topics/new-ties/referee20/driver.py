"""Bound rounds for kappa(P) = kappa(Q), P with p items, Q with q items (p, q in 1..3), P != Q.

Round(Hs) (Hs >= gap sum of both words):
  A  Arch(0,0): 3-log LLL -> min(H1(P), H1(Q)) <= X1      (vanishing case: Lemma N, part a)
  for each orientation (the word with H1 <= X1 is called P):
  B  Arch(1,0): for every head h = x_1 of P (gap <= X1) a 4-log LLL on log K(N h)
     -> min(H2(P), H1(Q)) <= X2(h)                           (non-vanishing: K(N h) != 1)
     B1: P has items 1,2 bounded;  B2: P item 1 and Q item 1 bounded.
  then a 2-adic tree from the back (exact discrete logs, Lemma N part b for vanishing) until one word
  has all items bounded ("short"); if too big, a Matveev cascade instead.
  A short word bounds the other one by size (derived_bound), using S <= SMAX.
"""
import math, sys, time, json
from engine import *

def gapsum(t):
    return sum(gap(it) for it in t)

class Budget(Exception):
    pass

def cascade(p, q, fP, fQ, YP, YQ, B):
    """Matveev-only bound for the rest of the tree.  Returns max candidate for new Hs."""
    if fP >= p:
        return max(YP, derived_bound(p, YP, q))
    if fQ >= q:
        return max(YQ, derived_bound(q, YQ, p))
    if fP == 0 and fQ == 0:
        X = matveev_h(3, None, B)
    else:
        Aa = max(logTmax(fP, YP) if fP else 0.0, logTmax(fQ, YQ) if fQ else 0.0)
        X = matveev_h(4, Aa, B)
    return max(cascade(p, q, fP + 1, fQ, YP + X, YQ, B), cascade(p, q, fP, fQ + 1, YP, YQ + X, B))

def cascade2(p, q, fP, fQ, YP, YQ, cP, cQ, tsP, tsQ, B):
    """Matveev cascade where P is complete once fP >= cP (its other items, a known-bounded tail with gap
    sum <= tsP, are already bounded), likewise Q."""
    if fP >= cP:
        box = YP + tsP
        return max(box, derived_bound(p, box, q))
    if fQ >= cQ:
        box = YQ + tsQ
        return max(box, derived_bound(q, box, p))
    if fP == 0 and fQ == 0:
        X = matveev_h(3, None, B)
    else:
        Aa = max(logTmax(fP, YP) if fP else 0.0, logTmax(fQ, YQ) if fQ else 0.0)
        X = matveev_h(4, Aa, B)
    return max(cascade2(p, q, fP + 1, fQ, YP + X, YQ, cP, cQ, tsP, tsQ, B),
               cascade2(p, q, fP, fQ + 1, YP, YQ + X, cP, cQ, tsP, tsQ, B))

def adic_tree(p, q, fP, fQ, FP, FQ, Mb, budget, vlog, B=None, maxdepth=99):
    maxP = maxQ = -1
    maxC = -1          # cascade fallbacks (already combined with derived bounds)
    stack = [((), ())]
    nodes = 0
    while stack:
        tP, tQ = stack.pop()
        if fP + len(tP) >= p:
            maxP = max(maxP, FP + gapsum(tP)); continue
        if fQ + len(tQ) >= q:
            maxQ = max(maxQ, FQ + gapsum(tQ)); continue
        nodes += 1
        if nodes > budget:
            raise Budget()
        bP, bQ = adic_bounds(p, q, tP, tQ, Mb, vlog)
        itP, itQ = items_g2(bP), items_g2(bQ)
        deep = len(tP) + len(tQ) >= maxdepth
        if itP:
            if fP + len(tP) + 1 >= p:     # adding one item completes P: record directly
                maxP = max(maxP, FP + gapsum(tP) + max(gap(x) for x in itP))
            elif deep:                    # P gets one more bounded tail item; finish by Matveev
                maxC = max(maxC, cascade2(p, q, fP, fQ, FP, FQ, p - len(tP) - 1, q - len(tQ),
                                          gapsum(tP) + max(gap(x) for x in itP), gapsum(tQ), B))
            else:
                for x in itP:
                    stack.append(((x,) + tP, tQ))
        if itQ:
            if fQ + len(tQ) + 1 >= q:
                maxQ = max(maxQ, FQ + gapsum(tQ) + max(gap(y) for y in itQ))
            elif deep:
                maxC = max(maxC, cascade2(p, q, fP, fQ, FP, FQ, p - len(tP), q - len(tQ) - 1,
                                          gapsum(tP), gapsum(tQ) + max(gap(y) for y in itQ), B))
            else:
                for y in itQ:
                    stack.append((tP, (y,) + tQ))
    return maxP, maxQ, nodes, maxC

def subtree(p, q, fP, fQ, FP, FQ, Hs, B, budget, out, vlog, tag):
    Mb = Hs + 3
    # enumerate tails fully only if the item lists are small; else stop at depth 1 and finish by Matveev
    V0 = 2 + int(math.log2(2 * Mb + 1)) + 3
    maxdepth = 99 if len(items_g2(V0)) <= 800 else 1
    try:
        mP, mQ, nodes, mC = adic_tree(p, q, fP, fQ, FP, FQ, Mb, budget, vlog, B=B, maxdepth=maxdepth)
        cands = [mC] if mC >= 0 else []
        if mP >= 0: cands.append(max(mP, derived_bound(p, mP, q)))
        if mQ >= 0: cands.append(max(mQ, derived_bound(q, mQ, p)))
        c = max(cands)
        out.append('   %s: 2-adic tree (depth %s), %d nodes; short P box sum <= %d, short Q box sum <= %d, Matveev fallbacks -> %s; -> %d'
                   % (tag, 'full' if maxdepth == 99 else maxdepth, nodes, mP, mQ, fmt(mC) if mC >= 0 else 'none', c))
        return c
    except Budget:
        c = cascade(p, q, fP, fQ, FP, FQ, B)
        out.append('   %s: 2-adic tree too big (> %d nodes); Matveev cascade -> %d' % (tag, budget, c))
        return c

def one_round(p, q, Hs, budget=400000, maxheads=200000, verbose=True):
    t0 = time.time()
    out = []
    B = 3 * Hs + 10
    X1, e1 = lll3(B)
    out.append('round: Hs=%s  X0=B=%s' % (fmt(Hs), fmt(B)))
    out.append('   A: 3-log LLL (C=10^%d): min(H1,H1\') <= %d' % (e1, X1))
    vlog = []
    cand = []
    orients = [(p, q)] if p == q else [(p, q), (q, p)]
    info = dict(X1=X1)
    for (pp, qq) in orients:
        if pp == 1:
            c = max(X1, derived_bound(1, X1, qq))
            out.append('   [P=%d items] P short at once: box %d -> %d' % (pp, X1, c))
            cand.append(c)
            continue
        heads = items_upto(X1, first=True)
        if len(heads) > maxheads:
            c = cascade(pp, qq, 1, 0, X1, 0, B)
            out.append('   [P=%d items] %d heads: too many, Matveev cascade from (1,0) -> %d' % (pp, len(heads), c))
            cand.append(c)
            continue
        X2 = {}
        es = set()
        for h in heads:
            K = strip(Tval([h]))
            assert K != 1
            X2[h], e = lll4(K, B)
            es.add(e)
        F2 = max(gap(h) + X2[h] for h in heads)
        X2max = max(X2.values())
        out.append('   [P=%d items, Q=%d] B: %d heads, 4-log LLL each (C=10^%s): X2max=%d, max(H1+X2)=%d'
                   % (pp, qq, len(heads), sorted(es), X2max, F2))
        info['X2max_%d%d' % (pp, qq)] = X2max
        info['F2_%d%d' % (pp, qq)] = F2
        if pp == 2:
            c = max(F2, derived_bound(2, F2, qq))
            out.append('   B1: P short, box %d -> %d' % (F2, c))
        else:
            c = subtree(pp, qq, 2, 0, F2, 0, Hs, B, budget, out, vlog, 'B1 (P items 1,2 bounded)')
        cand.append(c)
        if qq == 1:
            c = max(X2max, derived_bound(1, X2max, pp))
            out.append('   B2: Q short, box %d -> %d' % (X2max, c))
        else:
            c = subtree(pp, qq, 1, 1, X1, X2max, Hs, B, budget, out, vlog, 'B2 (P item 1, Q item 1 bounded)')
        cand.append(c)
    newHs = max(cand)
    vans = sorted(set((v[6], v[1] == v[2]) for v in vlog))
    out.append('   2-adic vanishing cases met: %d, resolutions %s' % (len(vlog), vans))
    out.append('   new Hs <= %d   (%.1fs)' % (newHs, time.time() - t0))
    if verbose:
        print('\n'.join(out), flush=True)
    return newHs, out, info

def fmt(x):
    return ('%.3e' % x) if x > 10 ** 6 else str(x)

def apriori(p, q):
    f = lambda H: cascade(p, q, 0, 0, 0, 0, 3 * H + 10)
    H = 10 ** 300
    for _ in range(200):
        nH = f(H)
        if nH >= H:
            break
        H = nH
    H0 = H
    # f(x) < x for all x >= H0+1, checked on a doubling grid up to 1e3000 (beyond, f(x) is a
    # polynomial in log x, hence f(2x) < x trivially).
    # f is nondecreasing; for x in [x_i, 2 x_i] we have f(x) <= f(2 x_i) < x_i <= x.
    x = H0 + 1
    while f(2 * x) >= x:
        x = x + x // 1000 + 1
    H0 = x
    while x < 10 ** 3000:
        assert f(2 * x) < x, x
        x = 2 * x
    return H0

def run(p, q, budget=400000, label=''):
    print('==== kappa(P)=kappa(Q): %d items vs %d items (%d vs %d letters) %s' % (p, q, p + 1, q + 1, label), flush=True)
    H0 = apriori(p, q)
    print('Matveev a priori (cascade): Hs <= %s' % fmt(H0), flush=True)
    Hs = H0
    hist = [H0]
    while True:
        nH, out, info = one_round(p, q, Hs, budget=budget)
        if nH >= Hs:
            break
        Hs = nH
        hist.append(Hs)
    print('rounds:', [fmt(h) for h in hist], flush=True)
    return Hs, hist

if __name__ == '__main__':
    p, q = int(sys.argv[1]), int(sys.argv[2])
    run(p, q)
