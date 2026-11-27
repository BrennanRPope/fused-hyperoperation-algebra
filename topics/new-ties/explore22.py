"""Round 22 exploration (DRAFT): Conjecture 21.2 for small sizes.

Coordinates (Proposition 22.2): every product of 2, 5 and 19 is uniquely 2^j (19/10)^E (1/20)^F, and a product of word
ratios is exactly such a number with j >= 1, F >= j, -j <= E <= F - j (j = how many word ratios). Call these the cone.

Part B  sides of 2 against 1 or 2: a balance a1 + a2 = b1 (+ b2) with a2 = rho a1, b2 = sigma b1 exists exactly when the
        fingerprint strip(n + d) of 1 + rho = (d + n)/d is 1, or equal for rho != sigma. Search the cone to 10^300.
Part C  sides with 3 numbers against 1, 2 or 3: fingerprints of 1 + rho1 + rho1 rho2, every collision tested for primitivity.

Usage: python3 explore22.py [B|C] [digits]   Writes explore22.json. Part B to 10^300: about 25 s, peak about 2 GB.
"""
import itertools, json, math, sys, time, resource

def strip(x):
    for p in (2, 5, 19):
        while x % p == 0:
            x //= p
    return x

def exps(j, E, F):
    """exponents of 2, 5, 19 in 2^j (19/10)^E (1/20)^F"""
    return j - E - 2 * F, -E - F, E

def nd(e2, e5, e19):
    n = 2 ** max(e2, 0) * 5 ** max(e5, 0) * 19 ** max(e19, 0)
    d = 2 ** max(-e2, 0) * 5 ** max(-e5, 0) * 19 ** max(-e19, 0)
    return n, d

def cone(X):
    """all cone elements rho = n/d (reduced) with max(n, d) <= X, as (n, d, j, E, F)"""
    out = []
    Fmax = int(math.log(X) / math.log(20 / 1.9)) + 3
    for F in range(1, Fmax + 1):
        for j in range(1, F + 1):
            for E in range(-j, F - j + 1):
                n, d = nd(*exps(j, E, F))
                if max(n, d) <= X:
                    out.append((n, d, j, E, F))
    return out

def primitive(A, B):
    """True if no nonempty proper subfamily of the numbers of A and B balances."""
    for i in range(0, len(A) + 1):
        for a in itertools.combinations(A, i):
            for k in range(0, len(B) + 1):
                if (i, k) in ((0, 0), (len(A), len(B))):
                    continue
                for b in itertools.combinations(B, k):
                    if sum(a) == sum(b):
                        return False
    return True

OUT = {}

def partB(digits=300):
    t = time.time()
    C = cone(10 ** digits)
    seen, ones, coll = {}, [], []
    for n, d, j, E, F in C:
        k = strip(n + d)
        if k == 1:
            ones.append([n, d, j, E, F])
        if k in seen:
            coll.append([seen[k], [n, d]])
        else:
            seen[k] = [n, d]
    OUT['B'] = dict(digits=digits, cone_elements=len(C), unit_sums=ones, collisions=coll, seconds=round(time.time() - t))
    print('B:', len(C), ones, len(coll))

def three_sums(X):
    """1 + r1 + r1 r2 over cone pairs with every number of the cleared balance at most X: (key, side as integers)"""
    C = cone(X)
    out = []
    for (n1, d1, *_r) in C:
        for (n2, d2, *_s) in C:
            if d1 * d2 > X:
                continue
            side = (d1 * d2, n1 * d2, n1 * n2)          # a1 = d1 d2, a2 = n1 d2, a3 = n1 n2 (then divide by gcd)
            g = math.gcd(math.gcd(*side[:2]), side[2])
            side = tuple(x // g for x in side)
            out.append((strip(sum(side)), side))
    return out, C

def partC(digits=30):
    t = time.time()
    X = 10 ** digits
    S3, C = three_sums(X)
    S2 = [(strip(n + d), (d, n)) for (n, d, *_r) in C]
    by = {}
    for k, side in S2 + S3:
        by.setdefault(k, []).append(side)
    found, nonprim, pairs = [], 0, 0
    for k, sides in by.items():
        lists = list(sides) + ([(1,)] if k == 1 else [])
        for A, B in itertools.combinations(lists, 2):
            if len(A) < 3 and len(B) < 3:
                continue                                  # 2 against 2 or 1 is Part B
            pairs += 1
            sA, sB = sum(A), sum(B)
            L = sA * sB // math.gcd(sA, sB)
            A2 = tuple(x * (L // sA) for x in A); B2 = tuple(x * (L // sB) for x in B)
            g = 0
            for x in A2 + B2:
                g = math.gcd(g, x)
            A2 = tuple(x // g for x in A2); B2 = tuple(x // g for x in B2)
            if set(A2) & set(B2) or not primitive(A2, B2):
                nonprim += 1
                continue
            found.append([list(A2), list(B2)])
    OUT['C'] = dict(digits=digits, cone_elements=len(C), three_sides=len(S3), equal_fingerprint_pairs=pairs,
                    not_primitive=nonprim, primitive_found=found[:20], n_primitive=len(found),
                    seconds=round(time.time() - t))
    print('C:', len(C), len(S3), pairs, nonprim, len(found), found[:5])

if __name__ == '__main__':
    try:
        OUT.update(json.load(open('explore22.json')))
    except Exception:
        pass
    part = sys.argv[1] if len(sys.argv) > 1 else 'B'
    digits = int(sys.argv[2]) if len(sys.argv) > 2 else (300 if part == 'B' else 30)
    {'B': partB, 'C': partC}[part](digits)
    json.dump(OUT, open('explore22.json', 'w'), indent=1, default=str)
    print('peak MB', round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
