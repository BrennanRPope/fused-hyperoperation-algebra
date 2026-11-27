"""Round 21 exploration (DRAFT): what "the powers" can and cannot say about every length.

Part A  the four letters' multipliers and their powers of 2 (Brennan's 1, 2, 4, 8), on the collision surface.
Part B  vanishing sums of products of 2, 5 and 19 whose two sides are each 10-lacunary: how many primitive ones.
Part C  the same, but each side "word-shaped" (consecutive ratios are products of word ratios): only 1 + 19 = 20.
Part E  the same balances with only the 2- and 5-sizes rule (they never grow toward the small end): 319 to 10^8.
Part D  Lemma 21.4: for canonical words N u, the pair (canonical form of u, t_u) determines u.

Usage: python3 explore21.py [A|B|C|D|E ...]   (no argument: all).  Writes explore21.json.
Part B peaks near 3.1 GB and Part C near 2.6 GB (Python dicts of all sums); both well under 11.5 GB.
"""
import itertools, json, math, sys, time, collections, resource
from fractions import Fraction as F

OUT = {}

def v(n, p):
    n = abs(n); k = 0
    while n and n % p == 0:
        n //= p; k += 1
    return k

def vq(q, p):
    q = F(q)
    return v(q.numerator, p) - v(q.denominator, p)

# ----------------------------------------------------------------------------- Part A
def surface(kappa, mu):
    kappa, mu = F(kappa), F(mu)
    return dict(nu=kappa / (2 - kappa), kappa=kappa, mu=mu, lam=mu * (2 - kappa))

def two_parts(br):
    """2-part of each letter's stretch 1/s, in the order N (nu), K (kappa), m (mu), l (lambda)."""
    return [2 ** (-vq(br[k], 2)) for k in ('nu', 'kappa', 'mu', 'lam')]

def partA():
    res = {}
    b1 = surface(F(1, 10), F(1, 20))           # section 13.4
    b3 = surface(F(3, 40), F(3, 80))           # section 16 (the 3-adic break)
    assert b1['lam'] == F(19, 200) and b1['nu'] == F(1, 19)
    assert b3['lam'] == F(231, 3200) and b3['nu'] == F(3, 77)
    res['break_13_4'] = dict(values={k: str(x) for k, x in b1.items()}, two_parts=two_parts(b1),
                             five_parts=[5 ** (-vq(b1[k], 5)) for k in ('nu', 'kappa', 'mu', 'lam')])
    res['break_3adic'] = dict(values={k: str(x) for k, x in b3.items()}, two_parts=two_parts(b3))
    # a grid of surface points: kappa = 1/(2^a q), mu = 1/(2^b q'), q, q' odd
    grid = []
    for a in range(0, 5):
        for b in range(0, 5):
            for q in (1, 3, 5, 7, 9):
                for q2 in (1, 3, 7):
                    kap = F(1, 2 ** a * q); m = F(1, 2 ** b * q2)
                    if kap >= 1 or m >= 1:
                        continue
                    br = surface(kap, m)
                    e = [-vq(br[k], 2) for k in ('nu', 'kappa', 'mu', 'lam')]
                    grid.append(dict(a=a, b=b, q=q, q2=q2, exps=e))
    ok_formula = all(g['exps'] == ([0, g['a'], g['b'], g['a'] + g['b']] if g['a'] >= 1 else
                                   [-vq(F(1, g['q']) / (2 - F(1, g['q'])), 2), 0, g['b'], g['b'] + (-vq(2 - F(1, g['q']), 2))])
                     for g in grid)
    one248 = sorted({(g['a'], g['b']) for g in grid if g['exps'] == [0, 1, 2, 3]})
    res['grid_points'] = len(grid)
    res['formula_0_a_b_a+b_holds_for_a>=1'] = ok_formula
    res['(a,b)_giving_1_2_4_8'] = one248
    OUT['A'] = res
    print('A:', res['break_13_4']['two_parts'], res['break_3adic']['two_parts'], len(grid), ok_formula, one248)

# ----------------------------------------------------------------------------- Parts B and C
def units(X):
    return sorted(2 ** a * 5 ** b * 19 ** c for a in range(int(math.log2(X)) + 1)
                  for b in range(int(math.log(X, 5)) + 1) for c in range(int(math.log(X, 19)) + 1)
                  if 2 ** a * 5 ** b * 19 ** c <= X)

def word_shaped(big, small):
    """small/big = 19^L / (19^n19 10^nK 200^L 20^M) for some L, M, n19, nK >= 0 with n19 + nK >= 1:
    a product of word ratios 19^L/(r 200^L 20^M), r in {19, 10} (the K-only-with-M = 0 rule is dropped: a relaxation)."""
    a = v(small, 2) - v(big, 2); b = v(small, 5) - v(big, 5); c = v(small, 19) - v(big, 19)
    if small * big == 0:
        return False
    LM = b - a
    if LM < 0:
        return False
    for L in range(LM + 1):
        M = LM - L; nK = -b - 2 * L - M; n19 = L - c
        if nK >= 0 and n19 >= 0 and nK + n19 >= 1:
            return True
    return False

def primitive_relations(sums):
    rels = set()
    for s, L in sums.items():
        if len(L) < 2:
            continue
        for A, B in itertools.combinations(L, 2):
            if set(A) & set(B):
                continue
            bad = False
            for i in range(1, len(A) + 1):
                for a in itertools.combinations(A, i):
                    sa = sum(a)
                    for j in range(1, len(B) + 1):
                        if (i, j) == (len(A), len(B)):
                            continue
                        for b in itertools.combinations(B, j):
                            if sa == sum(b):
                                bad = True; break
                        if bad: break
                    if bad: break
                if bad: break
            if bad:
                continue
            g = 0
            for x in A + B:
                g = math.gcd(g, x)
            rels.add(tuple(sorted([tuple(x // g for x in A), tuple(x // g for x in B)])))
    return rels

def search(X, kmax, shaped, rat=10):
    U = units(X)
    if shaped:
        nxt = {u: [w for w in U if w < u and word_shaped(u, w)] for u in U}
    else:
        nxt = {u: [w for w in U if w * rat <= u] for u in U}
    sums = collections.defaultdict(list)
    cnt = 0
    stack = [((u,), u) for u in U]
    while stack:
        t, s = stack.pop()
        sums[s].append(t); cnt += 1
        if len(t) < kmax:
            for w in nxt[t[-1]]:
                stack.append((t + (w,), s + w))
    rels = primitive_relations(sums)
    return len(U), cnt, rels

def side_ok(S):
    return all(word_shaped(S[i], S[i + 1]) for i in range(len(S) - 1))

def partB(X=10 ** 9, kmax=3):
    t = time.time()
    nU, cnt, rels = search(X, kmax, shaped=False)
    shaped = sorted(r for r in rels if side_ok(r[0]) and side_ok(r[1]))
    ex = sorted(rels, key=lambda r: (len(r[0]) + len(r[1]), max(r[0][0], r[1][0])))[:25]
    OUT['B'] = dict(X=X, kmax=kmax, units=nU, lacunary_sums=cnt, primitive=len(rels),
                    word_shaped_among_them=[list(map(list, r)) for r in shaped],
                    smallest_examples=[list(map(list, r)) for r in ex], seconds=round(time.time() - t))
    print('B:', nU, cnt, len(rels), shaped, ex[:8])

def partC(X=10 ** 15, kmax=4):
    t = time.time()
    nU, cnt, rels = search(X, kmax, shaped=True)
    OUT['C'] = dict(X=X, kmax=kmax, units=nU, word_shaped_sums=cnt,
                    primitive=[list(map(list, r)) for r in sorted(rels)], seconds=round(time.time() - t))
    print('C:', nU, cnt, sorted(rels))

def mono(S):
    """descending terms whose 2- and 5-sizes never grow toward the small end (true of every word's terms)."""
    return all(v(S[i + 1], 2) <= v(S[i], 2) and v(S[i + 1], 5) <= v(S[i], 5) for i in range(len(S) - 1))

def partE(X=10 ** 8, kmax=3):
    t = time.time()
    nU, cnt, rels = search(X, kmax, shaped=False)
    m = sorted((r for r in rels if mono(r[0]) and mono(r[1])), key=lambda r: (max(r[0][0], r[1][0]), r))
    OUT['E'] = dict(X=X, kmax=kmax, primitive=len(rels), sizes_rule_only=len(m),
                    smallest=[list(map(list, r)) for r in m[:15]], seconds=round(time.time() - t))
    print('E:', len(rels), len(m), m[:6])

# ----------------------------------------------------------------------------- Part D
S = {'l': F(19, 200), 'm': F(1, 20), 'N': F(1, 19), 'K': F(1, 10)}
D = {'l': 0, 'm': 0, 'N': 1, 'K': 1}

def st(w):
    """the map w(y) = s y + t, with xw = x o w."""
    s, t = F(1), F(0)
    for x in reversed(w):
        s, t = S[x] * s, S[x] * t + S[x] * D[x]
    return s, t

def canon(w):
    """section 19.3: drop leading l, m; first letter N; second letter N -> m; mK -> lN, ml -> lm; drop trailing l, m."""
    w = w.lstrip('lm')
    if not w:
        return ''
    w = 'N' + w[1:]
    if len(w) > 1 and w[1] == 'N':
        w = 'Nm' + w[2:]
    while True:
        w2 = w.replace('mK', 'lN').replace('ml', 'lm')
        if w2 == w:
            break
        w = w2
    return w.rstrip('lm')

def canonical_words(nmax):
    out = []
    for n in range(1, nmax + 1):
        for tup in itertools.product('lmNK', repeat=n - 1):
            w = 'N' + ''.join(tup)
            if w[-1] in 'NK' and canon(w) == w:
                out.append(w)
    return out

def partD(nmax=11):
    t = time.time()
    W = canonical_words(nmax)
    seen = {}
    clash = []
    for w in W:
        u = w[1:]
        if not u:
            continue
        key = (canon(u), st(u)[1])
        if key in seen and seen[key] != u:
            clash.append((seen[key], u))
        seen[key] = u
    OUT['D'] = dict(nmax=nmax, canonical_words=len(W), suffixes=len(seen), clashes=clash[:10],
                    seconds=round(time.time() - t))
    print('D:', len(W), len(seen), clash[:5])

if __name__ == '__main__':
    parts = sys.argv[1:] or ['A', 'B', 'C', 'D', 'E']
    try:
        OUT.update(json.load(open('explore21.json')))
    except Exception:
        pass
    for p in parts:
        {'A': partA, 'B': partB, 'C': partC, 'D': partD, 'E': partE}[p]()
        json.dump(OUT, open('explore21.json', 'w'), indent=1, default=str)
    print('peak MB', round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
