"""Item 2 controls. (a) Independent slow search, numbers <= 10^7, up to 3 a side, using the
Fraction-based ratio test and itertools (no shared code with search22.py besides core22).
(b) Positive control: same search with 'word-shaped' loosened to 'each ratio <= 1/10', to show
that the balance/primitivity code does find other primitive balances when they exist."""
import sys
from fractions import Fraction as F
from itertools import combinations
from math import gcd
from functools import reduce
from collections import defaultdict
sys.path.insert(0, '.')
from core22 import is_ratio_product

B = 10 ** 7
nums = sorted({2 ** a * 5 ** b * 19 ** c for a in range(25) for b in range(11) for c in range(6)
               if 2 ** a * 5 ** b * 19 ** c <= B}, reverse=True)

def primitive(Lt, Rt):
    allv = [(0, v) for v in Lt] + [(1, v) for v in Rt]
    n = len(allv)
    for mask in range(1, (1 << n) - 1):
        sl = sum(v for k, (sd, v) in enumerate(allv) if mask >> k & 1 and sd == 0)
        sr = sum(v for k, (sd, v) in enumerate(allv) if mask >> k & 1 and sd == 1)
        if sl == sr: return False
    return True

def run(shaped):
    lists = defaultdict(list)
    for k in (1, 2, 3):
        for c in combinations(nums, k):      # combinations of a descending list are descending
            if all(shaped(F(c[i + 1], c[i])) for i in range(k - 1)):
                lists[sum(c)].append(c)
    out = set()
    for s, ls in lists.items():
        for i in range(len(ls)):
            for j in range(i, len(ls)):
                a, b = ls[i], ls[j]
                if i == j and len(a) > 1: continue
                if reduce(gcd, a + b) != 1: continue
                if primitive(a, b): out.add(tuple(sorted((a, b))))
    return sorted(out)

ratio_cache = {}
def ws(q):
    if q not in ratio_cache: ratio_cache[q] = is_ratio_product(q)
    return ratio_cache[q]
print('numbers <= 10^7:', len(nums))
r = run(ws)
print('(a) word-shaped, <=3 a side, <=10^7: primitive gcd-1 balances:', r)
r2 = run(lambda q: q <= F(1, 10))
print('(b) loosened to ratios <= 1/10: primitive gcd-1 balances found:', len(r2), '; first few:', r2[:8])
bad = [x for x in r2 if not all(all(ws(F(s[i + 1], s[i])) for i in range(len(s) - 1)) for s in x)]
print('    of which not word-shaped on some side:', len(bad))
