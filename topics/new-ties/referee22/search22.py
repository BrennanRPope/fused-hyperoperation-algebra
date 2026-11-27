"""Item 2: search for primitive balances whose two sides are word-shaped.
Usage: python3 search22.py K BOUND   (K numbers a side at most; all numbers <= BOUND, gcd of all = 1)
Independent code; word-ratio test derived by hand (see core22.is_ratio_product) and
cross-checked here against a direct enumeration of products of word ratios."""
import sys, time, json
from fractions import Fraction as F
from itertools import combinations
from collections import defaultdict
sys.path.insert(0, '.')
from core22 import is_ratio_product

K = int(sys.argv[1]); BOUND = int(sys.argv[2])
t0 = time.time()

# --- cross-check the ratio test against brute-force products of word ratios ---
wr = set()
for L in range(0, 5):
    for M in range(0, 9):
        for r in (19, 10):
            wr.add(F(19 ** L, r * 200 ** L * 20 ** M))
prods = set(wr); frontier = set(wr)
for _ in range(3):
    new = set()
    for a in frontier:
        for b in wr:
            c = a * b
            if c.denominator < 10 ** 30: new.add(c)
    frontier = new - prods; prods |= new
# test: every product found passes; and random 2,5,19-rationals with small exponents agree
for q in prods:
    assert is_ratio_product(q), q
import random
random.seed(5)
mism = 0
for _ in range(20000):
    e2, e5, e19 = random.randint(-14, 2), random.randint(-9, 2), random.randint(-4, 4)
    q = F(2) ** e2 * F(5) ** e5 * F(19) ** e19
    # brute: is q in products of <=4 word ratios with L<=4, M<=6? only decidable if small
    if is_ratio_product(q) and q.denominator < 10 ** 12:
        # find decomposition: L,p,q from formulas, must lie in brute set if k<=4 and small
        k = 2 * (-e5) - (-e2) - e19
        if k <= 4 and q not in prods:
            mism += 1
    if not is_ratio_product(q) and q in prods:
        mism += 1
print('ratio-test cross-check mismatches:', mism, flush=True)

# --- numbers 2^a 5^b 19^c <= BOUND ---
nums = []
a = 0
while 2 ** a <= BOUND:
    b = 0
    while 2 ** a * 5 ** b <= BOUND:
        c = 0
        while 2 ** a * 5 ** b * 19 ** c <= BOUND:
            nums.append((2 ** a * 5 ** b * 19 ** c, (a, b, c))); c += 1
        b += 1
    a += 1
nums.sort(reverse=True)
idx = {n: i for i, (n, e) in enumerate(nums)}
print('numbers:', len(nums), flush=True)


def ok(ex, ey):
    # ratio y/x with exponents ey-ex; same test as is_ratio_product, integer form
    e2, e5, e19 = ey[0] - ex[0], ey[1] - ex[1], ey[2] - ex[2]
    A, B = -e2, -e5
    if 2 * B - A - e19 < 1: return False
    return max(0, e19) <= min(A - B, 2 * B - A)

succ = [[] for _ in nums]
for i, (x, ex) in enumerate(nums):
    for j in range(i + 1, len(nums)):
        y, ey = nums[j]
        if y * 10 > x: continue
        if ok(ex, ey): succ[i].append(j)
print('edges:', sum(map(len, succ)), flush=True)

# --- all word-shaped lists of size 1..K, keyed by sum ---
bysum = {}
cnt = 0
def enc(path):
    v = 0
    for k in path: v = v * 8192 + (k + 1)
    return v
def dec(v):
    out = []
    while v: out.append(v % 8192 - 1); v //= 8192
    return tuple(reversed(out))
def rec(path, s):
    global cnt
    e = enc(path); cnt += 1
    old = bysum.get(s)
    if old is None: bysum[s] = e
    elif isinstance(old, list): old.append(e)
    else: bysum[s] = [old, e]
    if len(path) == K: return
    for j in succ[path[-1]]:
        path.append(j); rec(path, s + nums[j][0]); path.pop()
for i in range(len(nums)):
    rec([i], nums[i][0])
print('word-shaped lists:', cnt, ' distinct sums:', len(bysum), ' t=%.1fs' % (time.time() - t0), flush=True)


def primitive(left, right):
    allv = [(0, v) for v in left] + [(1, v) for v in right]
    n = len(allv)
    for mask in range(1, (1 << n) - 1):
        sl = sr = 0
        for k in range(n):
            if mask >> k & 1:
                if allv[k][0] == 0: sl += allv[k][1]
                else: sr += allv[k][1]
        if sl == sr: return False
    return True

found = []
for s, lists in bysum.items():
    lists = [dec(v) for v in lists] if isinstance(lists, list) else [dec(lists)]
    if len(lists) < 2 and not (len(lists) == 1 and len(lists[0]) == 1):
        continue
    # pairs (unordered), including a list with itself only for single numbers
    for i in range(len(lists)):
        for j in range(i, len(lists)):
            A_, B_ = lists[i], lists[j]
            if i == j and len(A_) > 1: continue
            L = [nums[k][0] for k in A_]; R = [nums[k][0] for k in B_]
            exps = [nums[k][1] for k in A_ + B_]
            if min(e[0] for e in exps) or min(e[1] for e in exps) or min(e[2] for e in exps):
                continue  # gcd > 1
            if not primitive(L, R): continue
            found.append((sorted(L, reverse=True), sorted(R, reverse=True)))
print('primitive word-shaped balances with gcd 1:', len(found))
for f in found: print('  ', f)
expected = [([1], [1]), ([19, 1], [20]), ([20], [19, 1])]
others = [f for f in found if f not in expected]
print('examples other than {c}={c} and {19c,c}={20c}:', len(others), ' total time %.1fs' % (time.time() - t0))
json.dump({'K': K, 'BOUND': BOUND, 'numbers': len(nums), 'lists': cnt, 'found': found, 'others': others},
          open('search22_K%d_%d.json' % (K, len(str(BOUND)) - 1), 'w'))
