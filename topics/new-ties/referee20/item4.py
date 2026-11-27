"""Item 4: 2- and 5-parts of T for reduced words with 2, 3, 4 letters (n = 1, 2, 3 items).
Exact: an item is 'big' if its multiplier mu = r*200^L*20^M is 0 mod 2^7*5^3 (then T_k = 19^(A_k) mod
2^7 5^3, with A_k mod ord(19) arbitrary -> we try every residue); otherwise 'small' (finitely many:
L<=2, M<=3).  Words whose items are all small are enumerated exactly with big integers."""
import itertools
MOD = 2 ** 7 * 5 ** 3
def vp(x, p, cap=99):
    if x == 0: return cap
    v = 0
    while x % p == 0 and v < cap:
        x //= p; v += 1
    return v
o19 = 1
while pow(19, o19, MOD) != 1: o19 += 1
print('ord(19 mod 16000) =', o19)
def mu(L, M, y): return (19 if y == 'N' else 10) * 200 ** L * 20 ** M
small = []
for L in range(0, 6):
    for M in range(0, 8):
        for y in 'NK':
            if y == 'K' and M: continue
            if mu(L, M, y) % MOD != 0:
                small.append((L, M, y))
# completeness of the small list: every item outside it has mu = 0 mod MOD
for L in range(0, 30):
    for M in range(0, 30):
        for y in 'NK':
            if y == 'K' and M: continue
            if (L, M, y) not in small: assert mu(L, M, y) % MOD == 0
assert max(s[0] for s in small) <= 2 and max(s[1] for s in small) <= 3
print('small items:', len(small))
def Texact(items):
    Tv, A = 1, 0
    for (L, M, y) in items:
        A += L; Tv = 19 ** A + mu(L, M, y) * Tv
    return Tv
for n in (1, 2, 3):
    pairs = {}
    # A: all items small, exact
    for combo in itertools.product(small, repeat=n):
        if combo[0] == (0, 0, 'N'): continue
        Tv = Texact(combo)
        e = (vp(Tv, 2), vp(Tv, 5))
        pairs.setdefault(e, []).append(combo)
    exactA = dict(pairs)
    # B: last big item followed by j small items
    capped = set()
    for j in range(0, n):
        for a in range(o19):
            for combo in itertools.product(small, repeat=j):
                Tv, A = pow(19, a, MOD), a
                for (L, M, y) in combo:
                    A += L; Tv = (pow(19, A, MOD) + mu(L, M, y) * Tv) % MOD
                e = (vp(Tv, 2, 7), vp(Tv, 5, 3))
                if e[0] >= 7 or e[1] >= 3: capped.add(e)
                pairs.setdefault(e, []).append(('big', a, combo))
    vals = {e: 2 ** e[0] * 5 ** e[1] for e in pairs}
    mx = max(vals.values())
    print('n=%d items (%d letters): (e2,e5) pairs reached:' % (n, n + 1), sorted(pairs), ' capped:', capped)
    print('   max 2^e2 5^e5 = %d, reached by (all-small words only?)' % mx,
          [(e, len(pairs[e]), all(x[0] != 'big' for x in pairs[e])) for e in pairs if vals[e] == mx])
    if n == 3:
        for e in pairs:
            if vals[e] == mx:
                print('   words with max:', exactA.get(e), ' big-item paths:', [x for x in pairs[e] if x[0] == 'big'][:5])
    print('   max e2 = %d, max e5 = %d' % (max(e[0] for e in pairs), max(e[1] for e in pairs)))
