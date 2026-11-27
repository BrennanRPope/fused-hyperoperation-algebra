"""Items 1 and 3 (Lemmas F2, F3) and the 5-letter claim (item 9): exact 19-adic reachability.

T_0 = 1, T_k = 19^(A_k) + r_k P_k T_(k-1), P_k = 200^L_k 20^M_k, r = 19 (N) or 10 (K, only with M_k = 0).
We compute the exact set of residues T_n mod 19^E reachable by reduced words with n items, for every
choice of L_k, M_k >= 0.  This is a finite computation because:
  * 19^(A) mod 19^E only depends on min(A, E);
  * for an N item, M ranges over all integers >= 0, so P = 10^L 20^(L+M) ranges over the whole coset
    10^L <20>, and <20> mod 19^E is exactly the set of units = 1 mod 19 (order 19^(E-1), checked);
    so the multiplier set is {19 u : u unit, u = 10^L mod 19}: only L mod 18 and min(L, E) matter;
  * for a K item (M = 0) with L < E the multiplier is the single unit 10*200^L; for L >= E the
    additive term vanishes and the multipliers 10*200^L (L >= E) run over ALL units, because 200 is a
    primitive root mod 19^E (checked below).
  * the rule "first item is not a bare N" changes nothing mod 19^E: M = 19^(E-1) gives 20^M = 1.
Everything is done on boolean arrays indexed by residues, with the coset multiplications done exactly
by classes (v19, unit part mod 19), which is exact because a full coset times any residue x = 19^v w is
the full class {19^(v+1) w' : w' = c w mod 19}.
"""
import numpy as np, sys, time

def setup(E):
    M = 19 ** E
    assert pow(200, 18 * 19 ** (E - 2), M) != 1 and pow(200, 6 * 19 ** (E - 1), M) != 1 \
        and pow(200, 9 * 19 ** (E - 1), M) != 1, '200 not primitive'
    # order of 20 is 19^(E-1)
    assert pow(20, 19 ** (E - 1), M) == 1 and pow(20, 19 ** (E - 2), M) != 1
    x = np.arange(M, dtype=np.int64)
    val = np.zeros(M, dtype=np.int64)
    w = x.copy()
    val[0] = E
    nz = x != 0
    y = x.copy()
    for v in range(1, E):
        m = nz & (y % 19 == 0)
        val[m] += 1
        y[m] //= 19
    unitres = np.where(nz, y % 19, 0)
    classes = {}
    for v in range(E):
        for d in range(1, 19):
            classes[(v, d)] = np.nonzero(nz & (val == v) & (unitres == d))[0]
    return M, val, unitres, classes

def run(E, n, verbose=True):
    M, val, unitres, classes = setup(E)
    pw = [19 ** a % M for a in range(E)] + [0]
    # state: dict A(capped) -> bool array
    st = {0: np.zeros(M, dtype=bool)}
    st[0][1] = True
    for k in range(n):
        new = {}
        def add(A, arr):
            if A in new:
                new[A] |= arr
            else:
                new[A] = arr.copy()
        for A, arr in st.items():
            idx = np.nonzero(arr)[0]
            pairs = set(zip(val[idx].tolist(), unitres[idx].tolist()))
            vs = set(val[idx].tolist())
            # ---- N items: multiplier coset c = 10^L mod 19, additive 19^min(A+L,E)
            Ls = [(L, pow(10, L, 19), min(A + L, E)) for L in range(E)]
            Ls += [(None, c, E) for c in range(1, 19)]          # L >= E, any residue class mod 18
            for L, c, Ap in Ls:
                out = np.zeros(M, dtype=bool)
                for (v, d) in pairs:
                    if v + 1 >= E:
                        out[0] = True
                    else:
                        out[classes[(v + 1, c * d % 19)]] = True
                out = np.roll(out, pw[Ap])          # x -> x + 19^Ap
                add(Ap, out)
            # ---- K items, L < E: exact unit multiplier
            for L in range(E):
                mu = 10 * pow(200, L, M) % M
                Ap = min(A + L, E)
                out = np.zeros(M, dtype=bool)
                out[(mu * idx + pw[Ap]) % M] = True
                add(Ap, out)
            # ---- K items, L >= E: all units, additive 0
            out = np.zeros(M, dtype=bool)
            for v in vs:
                if v >= E:
                    out[0] = True
                else:
                    for d in range(1, 19):
                        out[classes[(v, d)]] = True
            add(E, out)
        st = new
    tot = np.zeros(M, dtype=bool)
    for arr in st.values():
        tot |= arr
    vals = val[np.nonzero(tot)[0]]
    hist = {int(v): int((vals == v).sum()) for v in sorted(set(vals.tolist()))}
    return hist

if __name__ == '__main__':
    t = time.time()
    # sanity: compare with brute force on small words
    def T(items):
        Tv, A = 1, 0
        for (L, Mm, y) in items:
            A += L
            Tv = 19 ** A + (19 if y == 'N' else 10) * 200 ** L * 20 ** Mm * Tv
        return Tv
    def vp(x, p):
        v = 0
        while x % p == 0:
            x //= p; v += 1
        return v
    import itertools
    its = [(L, Mm, y) for L in range(6) for Mm in range(6) for y in 'NK' if not (y == 'K' and Mm)]
    for n, E in ((1, 3), (2, 4), (3, 5)):
        hist = run(E, n)
        mx = max(v for v in hist if v < E) if any(v < E for v in hist) else None
        print('n=%d items (%d letters), mod 19^%d: reachable v19 classes (v: #residues) %s' % (n, n + 1, E, hist))
        zero = E in hist
        print('   zero mod 19^%d reachable: %s  -> max v19(T_%d) = %s' % (E, zero, n, ('>= %d' % E) if zero else mx))
    # brute-force cross-check of maxima on small words
    for n in (1, 2, 3):
        best = 0
        cnt = 0
        for combo in itertools.product(its, repeat=n):
            if combo[0] == (0, 0, 'N'):
                continue
            if n == 3 and sum(c[0] + c[1] for c in combo) > 9:
                continue
            best = max(best, vp(T(list(combo)), 19)); cnt += 1
        print('brute force n=%d (%d words, small L,M): max v19 = %d' % (n, cnt, best))
    print('time %.1fs' % (time.time() - t))
