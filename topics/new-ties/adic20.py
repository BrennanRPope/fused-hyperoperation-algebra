"""Round 20 (DRAFT): the 19-adic step for five N/K letters.

Usage: python3 adic20.py [all|check]   (all: writes adic20.json, a few minutes; check: quick recomputation, under a minute)

What it does (README section 24, notes in adic20_notes.md):
  1. Lists the 12 five-letter tie patterns whose leading units cancel modulo 19 and shows each one has the same shape:
     W = N x1 x2 x3 N with the last item a bare-gap N (L4 = 0) carrying m^M, and
         T(W) = 19^A3 * (1 + 20^M * u),   u = T(N x1 x2 x3) / 19^(A3 - 1),   u = -1 mod 19,
     so g = v19(T(W)) = A3 + 1 + v19(M - mu(u)), where mu(u) is the 19-adic number with 20^mu = -1/u.
     The free parameters are M and the m-counts of the N-items among x1, x2, x3 (the L's and letters are fixed).
  2. Computes mu(u) modulo 19^k exactly, digit by digit (an exact 19-adic discrete logarithm to base 20), checked by pow.
  3. The digit argument: for 0 <= M <= B, the largest v19(M - mu) is the largest k with (mu mod 19^k) <= B, reached at
     M = mu mod 19^k. Exact maximum of g over all M <= B, for every front with free gap sum <= G.
  4. The two-logarithm 19-adic bound of Bugeaud and Laurent (1996), as restated in arXiv 2401.06555, Lemma 2.2, for all M.
"""
import sys, json, math, time, itertools, random

P19 = 19
LOG19 = math.log(19.0)
LOG20 = math.log(20.0)


def v(a, p):
    a = abs(a); k = 0
    while a % p == 0: a //= p; k += 1
    return k


def P_of(L, M): return 200 ** L * 20 ** M


def T_of(items):
    T = 1; A = 0
    for (r, L, M) in items:
        A += L; T = 19 ** A + r * P_of(L, M) * T
    return T


def name(items):
    return "N" + "".join("l" * L + "m" * M + ("N" if r == 19 else "K") for (r, L, M) in items)


# ---------------------------------------------------------------- 1. the cancelling tie patterns
def casework(n, Lmax):
    """Same casework as checks19.py Q1: term j of T has 19-adic size A_j + #(N after j) and leading unit
    prod_(i > j) c_i 10^(L_i) mod 19 (c = 1 for N, 10 for K). Returns (number of tie patterns, cancelling patterns)."""
    vr = lambda r: 1 if r == 19 else 0
    cr = lambda r: 1 if r == 19 else 10
    ties = 0; bad = []
    for rs in itertools.product((19, 10), repeat=n):
        for Ls in itertools.product(range(Lmax + 1), repeat=n):
            A = [0] + list(itertools.accumulate(Ls))
            t = [A[j] + sum(vr(rs[i]) for i in range(j, n)) for j in range(n + 1)]
            u = []
            for j in range(n + 1):
                x = 1
                for i in range(j, n): x = x * cr(rs[i]) * pow(10, Ls[i], 19) % 19
                u.append(x)
            g = min(t); tied = [j for j in range(n + 1) if t[j] == g]
            if len(tied) > 1:
                ties += 1
                if sum(u[j] for j in tied) % 19 == 0:
                    bad.append(dict(rs=rs, Ls=Ls, tmin=g, tied=tied))
    return ties, bad


def pattern_front(pat, Ms):
    """front items x1, x2, x3 of a pattern with free m-counts Ms for its N-items (K-items have M = 0)."""
    it = []; k = 0
    for r, L in zip(pat["rs"][:3], pat["Ls"][:3]):
        if r == 19: it.append((19, L, Ms[k])); k += 1
        else: it.append((10, L, 0))
    return tuple(it)


def n_free(pat): return sum(1 for r in pat["rs"][:3] if r == 19)


def fronts(pat, G):
    """all fronts of a pattern whose free m-counts add up to at most G"""
    k = n_free(pat)
    for Ms in itertools.product(range(G + 1), repeat=k):
        if sum(Ms) <= G: yield Ms, pattern_front(pat, Ms)


def u_of(front):
    A3 = sum(L for (_, L, _) in front)
    T3 = T_of(front)
    assert v(T3, 19) == A3 - 1, (front, v(T3, 19), A3)
    u = T3 // 19 ** (A3 - 1)
    assert u % 19 == 18
    a = u
    for q in (2, 5): 
        while a % q == 0: a //= q
    assert a > 1   # u has a prime factor other than 2, 5, so 20 and -u are multiplicatively independent
    return A3, u


# ---------------------------------------------------------------- 2. mu(u) modulo 19^k, exactly
def mu_digits(u, K):
    """mu mod 19^K, where 20^mu = w := -1/u in Z_19 (w = 1 mod 19). Digit j: z = w * 20^(-x) = 1 mod 19^(j+1);
    20^(19^j) = 1 + c_j 19^(j+1) mod 19^(j+2) with c_j a unit, so digit d = ((z - 1)/19^(j+1)) / c_j mod 19."""
    mod = 19 ** (K + 1)
    w = (-pow(u, -1, mod)) % mod
    assert w % 19 == 1
    z = w; x = 0
    base = 20 % mod          # 20^(19^j) mod 19^(K+1)
    binv = pow(20, -1, mod)  # its inverse
    for j in range(K):
        pj1 = 19 ** (j + 1)
        assert (z - 1) % pj1 == 0
        c = ((base - 1) // pj1) % 19
        d = (((z - 1) // pj1) % 19) * pow(c, -1, 19) % 19
        if d:
            x += d * 19 ** j
            z = z * pow(binv, d, mod) % mod
        base = pow(base, 19, mod); binv = pow(binv, 19, mod)
    assert pow(20, x, mod) * pow(w, -1, mod) % (19 ** (K + 1)) == 1   # 20^x = w mod 19^(K+1)
    return x


def max_v_upto(mu, K, B):
    """max over 0 <= M <= B of v19(M - mu), given mu mod 19^K (needs 19^K > B so the answer is < K).
    Returns (kmax, M at which it is reached, exact v19 there)."""
    assert 19 ** K > B
    kmax = 0
    for k in range(K + 1):
        if mu % 19 ** k <= B: kmax = k
    M = mu % 19 ** kmax
    vv = v(M - mu, 19) if M != mu else K  # exact v19(M - mu) is the first nonzero digit of mu above kmax
    vv = min(vv, K)
    assert vv >= kmax and vv < K
    # by maximality vv == kmax unless mu mod 19^(kmax+1) == M, which would contradict kmax being the largest
    assert vv == kmax
    return kmax, M


def K_for(B): return int(math.log(B) / LOG19) + 8


# ---------------------------------------------------------------- 4. Bugeaud-Laurent, two logarithms, p = 19
def bl_bound(u, M):
    """v19(gamma1^b1 gamma2^b2 - 1) <= 24 p g / ((p - 1) (log p)^4) * E^2 * D^4 * h'(gamma1) h'(gamma2), with
    gamma1 = 20, b1 = M, gamma2 = -u, b2 = 1, D = 1, g = 1 (both are 1 mod 19), h'(gamma) >= max(h(gamma), log p),
    E' = b1/h'(gamma2) + b2/h'(gamma1), E = max(log E' + log log p + 0.4, 10, 10 log p).
    Returns the bound on v19(1 + 20^M u) (= v19(20^M (-u) - 1)) as a float (rounded up)."""
    h1 = max(LOG20, LOG19); h2 = max(math.log(abs(u)), LOG19)
    Ep = M / h2 + 1 / h1
    E = max(math.log(Ep) + math.log(LOG19) + 0.4, 10.0, 10 * LOG19)
    C = 24 * 19 * 1 / (18 * LOG19 ** 4)
    return C * E * E * h1 * h2


# The original statement, as quoted in some later papers, has p^f - 1 in place of g (here 18 in place of 1): a
# conservative reading multiplies every bound below by 18. Both restatements we read (arXiv 2401.06555 Lemma 2.2 and
# arXiv 2606.18500 Theorem 2.3) have g; we did not see the 1996 paper itself.
CONSERVATIVE_FACTOR = 18


def bl_bound_g(A3, u, M):
    return A3 + bl_bound(u, M)


# ---------------------------------------------------------------- the runs
def identity_check(pats, nrand=300, seed=20):
    """g(W) = A3 + 1 + v19(M - mu(u)) on random words of every pattern (T computed exactly)"""
    rnd = random.Random(seed); n = 0
    for pat in pats:
        for _ in range(nrand // len(pats) + 1):
            Ms = tuple(rnd.randint(0, 12) for _ in range(n_free(pat)))
            front = pattern_front(pat, Ms)
            A3, u = u_of(front)
            M = rnd.choice([rnd.randint(0, 60), rnd.randint(0, 3000)])
            W = front + ((19, 0, M),)
            g = v(T_of(W), 19)
            mu = mu_digits(u, 12)
            pred = A3 + 1 + v(M - mu, 19) if v(M - mu, 19) < 12 else None
            assert pred is not None and g == pred, (name(W), g, pred)
            # the identity T = 19^A3 (1 + 20^M u)
            assert T_of(W) == 19 ** A3 * (1 + 20 ** M * u)
            n += 1
    return n


def brute_check(u, A3, Bsmall, K=20):
    """max g over M <= Bsmall by brute force (T modulo 19^(A3+K)), against the digit argument"""
    mod = 19 ** (A3 + K)
    best = (-1, None); p = 1
    base = 19 ** A3
    for M in range(Bsmall + 1):
        T = (base + base * p * u) % mod
        g = v(T, 19) if T else A3 + K
        assert g < A3 + K
        if g > best[0]: best = (g, M)
        p = p * 20 % mod
    mu = mu_digits(u, K)
    k, Mk = max_v_upto(mu, K, Bsmall)
    return best, (A3 + 1 + k, Mk)


def main(mode):
    t0 = time.time(); out = {}
    ties5, bad5 = casework(4, 5)
    out["patterns"] = dict(ties=ties5, cancel=len(bad5), list=[])
    for pat in bad5:
        assert pat["tied"] == [0, 1, 4] and pat["rs"][3] == 19 and pat["Ls"][3] == 0
        assert sum(pat["Ls"]) == pat["tmin"]
        front0 = pattern_front(pat, (0,) * n_free(pat))
        A3, u0 = u_of(front0)
        out["patterns"]["list"].append(dict(
            word=name(front0) + "m^M N" if n_free(pat) == 0 else name(front0) + " m^M N  (free: m-counts of the N-items)",
            rs=pat["rs"], Ls=pat["Ls"], A3=A3, free_front_gaps=n_free(pat), u_at_zero=str(u0),
            form="g = %d + v19(M - mu(u)),  20^mu = -1/u" % (A3 + 1)))
        print(out["patterns"]["list"][-1], flush=True)
    nid = identity_check(bad5, 300 if mode == "all" else 60)
    out["identity_words"] = nid
    print("identity checked on", nid, "random words", flush=True)

    # the README family N K lK K m^M N
    fam = next(p for p in bad5 if p["rs"] == (10, 10, 10, 19))
    A3, u = u_of(pattern_front(fam, ()))
    mu = mu_digits(u, 30)
    out["family"] = dict(word="N K l K K m^M N", A3=A3, u=str(u), u_mod_361=u % 361,
                         mu_mod_19_10=mu % 19 ** 10, mu_digits_low_first=[(mu // 19 ** j) % 19 for j in range(30)])
    gk = {M: v(T_of(pattern_front(fam, ()) + ((19, 0, M),)), 19) for M in (0, 38, 5092)}
    out["family"]["g_at_0_38_5092"] = [gk[0], gk[38], gk[5092]]
    # record M where g first reaches each value: M = mu mod 19^k
    out["family"]["first_M_for_g"] = {str(A3 + 1 + k): int(mu % 19 ** k) for k in range(0, 12)}
    print("family", out["family"], flush=True)
    bf = brute_check(u, A3, 10 ** 5 if mode == "all" else 20000)
    out["family"]["brute_force"] = dict(B=10 ** 5 if mode == "all" else 20000, brute=bf[0], digits=bf[1])
    assert bf[0][0] == bf[1][0], bf
    print("brute force vs digits", bf, flush=True)

    # exact max g over M <= B for every front, every pattern
    G = 30 if mode == "all" else 8
    Bs = [10 ** 2, 10 ** 100, 10 ** 200, 10 ** 300]
    Kmax = K_for(Bs[-1])
    res = {str(b): dict(max_g=0, at=None, hist={}) for b in [math.floor(math.log10(x)) for x in Bs]}
    nfront = 0
    blmax = {}
    for pat in bad5:
        for Ms, front in fronts(pat, G):
            A3, u = u_of(front)
            mu = mu_digits(u, Kmax)
            nfront += 1
            for B in Bs:
                k, Mk = max_v_upto(mu, K_for(B), B)
                g = A3 + 1 + k
                key = str(math.floor(math.log10(B)))
                r = res[key]
                r["hist"][str(g - 1 - math.ceil(math.log(B) / LOG19) - A3)] = r["hist"].get(
                    str(g - 1 - math.ceil(math.log(B) / LOG19) - A3), 0) + 1
                if g > r["max_g"]: r["max_g"] = g; r["at"] = dict(word=name(front) + " m^%s N" % Mk if B < 10 ** 30
                                                                     else name(front) + " m^M N", M_digits=len(str(Mk)), A3=A3)
            # Bugeaud-Laurent bound at the same B's, for comparison (all M <= B)
            for B in Bs:
                key = str(math.floor(math.log10(B)))
                bb = bl_bound_g(A3, u, B)
                blmax[key] = max(blmax.get(key, 0), bb)
    out["reduction"] = dict(G=G, fronts=nfront, by_B=res,
                            note="hist keys: k_max - ceil(log_19 B), i.e. how far the exact maximum is above the "
                                 "generic log_19 B; max_g is A3 + 1 + k_max")
    out["bugeaud_laurent_max_over_fronts"] = {k: round(x, 1) for k, x in blmax.items()}
    # BL bound as a function of the front height, and of M, for the family
    out["bl_family"] = {("log10 M = %d" % e): round(bl_bound_g(A3_f := u_of(pattern_front(fam, ()))[0],
                                                               u_of(pattern_front(fam, ()))[1], 10 ** e), 1)
                        for e in (0, 10, 12, 13, 20, 50, 100, 200, 300)}
    pat0 = next(p for p in bad5 if p["rs"] == (19, 19, 19, 19) and p["Ls"] == (1, 2, 1, 0))
    out["bl_vs_front"] = {}
    for s in (0, 10, 30, 100, 300, 1000):
        A3s, us = u_of(pattern_front(pat0, (s, 0, 0)))
        b = bl_bound_g(A3s, us, 10 ** 100)
        out["bl_vs_front"][str(s)] = dict(log_u=round(math.log(us), 1), bl_g_at_M_1e100=round(b, 1),
                                          per_unit_log_u=round((b - A3s) / math.log(us), 1),
                                          exact_max_g_M_upto_1e100=A3s + 1 + max_v_upto(mu_digits(us, K_for(10 ** 100)),
                                                                                     K_for(10 ** 100), 10 ** 100)[0])
    out["time_s"] = round(time.time() - t0, 1)
    print(json.dumps({k: out[k] for k in ("reduction", "bugeaud_laurent_max_over_fronts", "bl_family", "bl_vs_front")},
                     indent=1), flush=True)
    return out


def check():
    """quick recomputation of the key claims (under a minute): returns list of (name, ok, detail)"""
    r = []
    ties5, bad5 = casework(4, 5)
    ok = len(bad5) == 12 and all(p["tied"] == [0, 1, 4] and p["rs"][3] == 19 and p["Ls"][3] == 0 for p in bad5)
    r.append(("12 cancelling patterns, all of shape N x1 x2 x3 m^M N with tied terms 0, 1, 4", ok, "%d of %d" % (len(bad5), ties5)))
    n = identity_check(bad5, 48)
    r.append(("g = A3 + 1 + v19(M - mu(u)) and T = 19^A3 (1 + 20^M u)", True, "%d random words" % n))
    fam = next(p for p in bad5 if p["rs"] == (10, 10, 10, 19))
    A3, u = u_of(pattern_front(fam, ()))
    bf = brute_check(u, A3, 20000)
    r.append(("digit argument = brute force, family, M <= 20000", bf[0][0] == bf[1][0], str(bf)))
    mu = mu_digits(u, K_for(10 ** 300))
    ks = {e: max_v_upto(mu, K_for(10 ** e), 10 ** e)[0] for e in (100, 200, 300)}
    r.append(("family: exact max g for M <= 10^100, 10^200, 10^300", True,
              str({e: A3 + 1 + k for e, k in ks.items()})))
    return r, dict(A3=A3, u=u, ks=ks)


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "all"
    if mode == "check":
        rr, _ = check()
        for x in rr: print(("PASS " if x[1] else "FAIL ") + x[0] + ": " + x[2])
    else:
        out = main(mode)
        json.dump(out, open("adic20.json", "w"), indent=1, default=str)
        print("time", out["time_s"])
