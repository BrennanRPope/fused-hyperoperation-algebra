"""Round 21 checks (DRAFT). Run from this folder: python3 checks21.py   (about 15 seconds, peak about 0.2 GB; reads explore21.json)."""
import itertools, json, random, time
from fractions import Fraction as F
import explore21 as E

out, res = [], {}
def rep(name, ok, msg):
    out.append(f"{name}: {'PASS' if ok else 'FAIL'}  {msg}"); res[name] = dict(ok=ok, msg=msg); print(out[-1])

# Q1: the four letters' powers of 2
b1 = E.surface(F(1, 10), F(1, 20)); b3 = E.surface(F(3, 40), F(3, 80))
t1, t3 = E.two_parts(b1), E.two_parts(b3)
grid_ok, only, npts = True, set(), 0
for a in range(1, 7):
    for b in range(0, 7):
        for q in (1, 3, 5, 7, 11):
            for q2 in (1, 3, 7, 9):
                kap, mu = F(1, 2 ** a * q), F(1, 2 ** b * q2)
                if mu >= 1: continue
                br = E.surface(kap, mu); npts += 1
                e = [-E.vq(br[k], 2) for k in ('nu', 'kappa', 'mu', 'lam')]
                grid_ok &= (e == [0, a, b, a + b])
                if e == [0, 1, 2, 3]: only.add((a, b))
rep('Q1', t1 == [1, 2, 4, 8] and t3 == [1, 8, 16, 128] and grid_ok and only == {(1, 2)},
    f"2-parts of N, K, m, l: {t1} at the 13.4 break, {t3} at the 3-adic break; on the surface the exponents are "
    f"(0, a, b, a+b) on all {npts} grid points; 1, 2, 4, 8 only for (a, b) = {sorted(only)}")

# Q2: nothing changes with length: N N w ~ N m w at every length, and no clash among canonical words to 11 letters
random.seed(21)
okp = True
for n in range(0, 25):
    for _ in range(20):
        w = ''.join(random.choice('lmNK') for _ in range(n)) + random.choice('NK')
        s1, t1_ = E.st('NN' + w); s2, t2_ = E.st('Nm' + w)
        okp &= (t1_ / s1 == t2_ / s2)
def fp(w):
    s, t = E.st(w); q = t / s
    n, d = q.numerator, q.denominator
    for p in (2, 5, 19):
        while n % p == 0: n //= p
        while d % p == 0: d //= p
    return F(n, d)
W = E.canonical_words(11)
bylen, seen, clash = {}, {}, 0
for w in W:
    f = fp(w)
    if f in seen: clash += 1
    seen[f] = w
    bylen[len(w)] = bylen.get(len(w), 0) + 1
rep('Q2', okp and clash == 0,
    f"psi(NNw) = psi(Nmw) for 20 random w at each length 1 to 25; {len(W):,} canonical words up to 11 letters, "
    f"{len(seen):,} fingerprints, no clash; per length {bylen}")

# Q3: primitive lacunary relations are plentiful without the word shape
ex = [((19, 1), (20,)), ((76, 5), (80, 1)), ((125, 4), (128, 1)), ((361, 20), (380, 1)), ((380, 19, 1), (400,))]
ok3 = all(sum(A) == sum(B) for A, B in ex)
ok3 &= [E.side_ok(A) and E.side_ok(B) for A, B in ex] == [True, False, False, False, False]
nU, cnt, rels = E.search(10 ** 6, 3, shaped=False)
shaped = [r for r in rels if E.side_ok(r[0]) and E.side_ok(r[1])]
monos = [r for r in rels if E.mono(r[0]) and E.mono(r[1])]
J = json.load(open('explore21.json'))
ok3 &= J['B']['primitive'] == 18852 and len(rels) > 100 and shaped == [((19, 1), (20,))]
ok3 &= len(monos) > 1 and J['E']['sizes_rule_only'] == 319 and ((380, 19, 1), (400,)) in monos
rep('Q3', ok3, f"5 sample relations hold and only 19 + 1 = 20 is word-shaped; to 10^6 with 3 terms a side: {len(rels):,} "
    f"primitive lacunary relations, {len(monos)} obeying the 2- and 5-sizes rule, 1 word-shaped; explore21.json: "
    f"{J['B']['primitive']:,} to 10^9, and {J['E']['sizes_rule_only']} obeying the sizes rule to 10^8")

# Q4: word-shaped relations
nU, cnt, rels = E.search(10 ** 12, 3, shaped=True)
okc = J['C']['primitive'] == [[[19, 1], [20]]] and J['C']['X'] == 10 ** 15 and J['C']['kmax'] == 4
rep('Q4', sorted(rels) == [((19, 1), (20,))] and okc,
    f"word-shaped primitive relations to 10^12 with 3 terms a side: {sorted(rels)} ({cnt:,} sums); explore21.json to "
    f"10^15 with 4 a side: {J['C']['primitive']} ({J['C']['word_shaped_sums']:,} sums)")

# Q5: the facts Lemma 21.2 uses
def taus(w):
    """the terms of psi(w): for each N or K at position i, the product of 1/s over the letters after it."""
    out = []
    for i, x in enumerate(w):
        if x in 'NK':
            r = F(1)
            for y in w[i + 1:]: r /= E.S[y]
            out.append(r)
    return out
only19 = []
for n in range(0, 9):
    for mid in itertools.product('lm', repeat=n):
        for y in 'NK':
            r = F(1)
            for x in mid + (y,): r *= E.S[x]
            if r == F(1, 19): only19.append(''.join(mid) + y)
okshape = True
for w in random.sample(W, 3000):
    T = taus(w); g = 1
    for x in T: g *= x.denominator
    Ti = [int(x * g) for x in T]
    okshape &= all(E.word_shaped(Ti[i], Ti[i + 1]) for i in range(len(Ti) - 1))
    okshape &= all(Ti[i] >= 10 * Ti[i + 1] for i in range(len(Ti) - 1))
rep('Q5', only19 == ['N'] and max(E.S.values()) == F(1, 10) and okshape,
    f"a ratio 1/19 between neighbouring terms comes only from {only19} (runs up to 8 l/m); every ratio is at most 1/10; "
    f"the terms of 3,000 canonical words are word-shaped and 10-lacunary")

# Q6: Lemma 21.3: the moves equation and the injectivity check
lam, mu = F(19, 200), F(1, 20)
sols = sorted((a, b, c, d) for a in range(-20, 21) for b in range(-20, 21) for c in (-1, 0, 1) for d in (-1, 0, 1)
              if lam ** a * mu ** b * F(19, 10) ** c * F(20, 19) ** d == 1)
seen6, cl6 = {}, 0
for w in E.canonical_words(10):
    u = w[1:]
    if not u: continue
    k = (E.canon(u), E.st(u)[1])
    if k in seen6 and seen6[k] != u: cl6 += 1
    seen6[k] = u
rep('Q6', sols == [(-1, 1, 1, 0), (0, 0, 0, 0), (1, -1, -1, 0)] and cl6 == 0,
    f"lam^a mu^b (19/10)^c (20/19)^d = 1 only for {sols}; (canonical form of u, t_u) determines u for all "
    f"{len(seen6):,} canonical words N u up to 10 letters")

out.append(f"{sum(r['ok'] for r in res.values())}/{len(res)} checks pass")
print(out[-1])
open('checks21_output.txt', 'w').write('\n'.join(out) + '\n')
json.dump(res, open('checks21.json', 'w'), indent=1)
