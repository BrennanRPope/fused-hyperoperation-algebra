"""Round 16 (DRAFT): Theorem 16.2, Conjecture 15.3 for every pair of words with at most two N/K letters. Computer-assisted,
by the same method as baker16.py plus one 2-adic step.

A reduced word with two N/K letters is N l^L m^M x (x = K only when M = 0, since mK = lN). Its psi is T/19^L with
T = 19^L + r P, P = 2^(3L+2M) 5^(2L+M), r = 19 (x = N) or 10 (x = K); its fingerprint is K = T / 19^g, g = v19(T)
(T is prime to 10 except for NN, whose fingerprint is 1, settled by Theorem 15.5). Suppose two different words share K.
Step 1: K = 19^(L-g) + rP/19^g = 19^(L'-g') + r'P'/19^(g'). If D := (L-g) - (L'-g') = 0 the second terms agree, and
        unique factorisation gives the same word or the pair mK ~ lN, which is not reduced. So D != 0.
Step 2 (2-adic): K = 19^(L-g) mod 2^mu with mu = v2(rP), and likewise for the other word, so 19^D = 1 mod 2^min(mu, mu'),
        hence min(mu, mu') <= v2(19^|D| - 1) <= 2 + log2 |D|. The word with the smaller mu ("small") has
        H_s = L_s + M_s <= mu_s / 2 <= 1 + log2(|D|)/2, and |D| <= 3 H + 2 with H the larger of the two H's.
Step 3: for the other ("big") word, Lambda = log(19^g K / (r P)) = log(1 + 19^L/(rP)) lies in (0, 10^-H / 10): a linear form
        in log 19, log 2, log 5 and log K with coefficients at most 3H + 3. Matveev (n = 4) bounds H; for each small word,
        LLL in dimension 4 cuts the bound; then every pair below the bound is compared exactly.
Usage: python3 pairs16.py   Writes pairs16.json.
"""
import json, math, time
from decimal import Decimal
from baker16 import matveev_constant, largest_ok, reduce_once, LOG
def strip(a):
    for p in (2, 5, 19):
        while a % p == 0: a //= p
    return a
def words(Hmax):
    for H in range(Hmax + 1):
        for L in range(H + 1):
            M = H - L
            yield (19, L, M)
            if M == 0: yield (10, L, M)
def T_of(w):
    r, L, M = w
    return 19 ** L + r * 2 ** (3 * L + 2 * M) * 5 ** (2 * L + M)
def name(w):
    r, L, M = w
    return "N" + "l" * L + "m" * M + ("N" if r == 19 else "K")
if __name__ == "__main__":
    t0 = time.time(); out = {}
    C3 = matveev_constant(n=4, alphas=(19, 2, 5))          # times A4 = max(log K, 0.16), supplied below
    logK_bound = lambda H: 3.64 + 5.30 * (1 + math.log2(3 * H + 2) / 2)
    H0 = largest_ok(lambda H: H * math.log(10) + math.log(10) < C3 * (1 + math.log(3 * H + 3)) * logK_bound(H), 3, 10 ** 20)
    Hs_max = math.floor(1 + math.log2(3 * H0 + 2) / 2)
    out["matveev"] = dict(constant_without_A4=C3, H_big_bound=H0, H_small_bound=Hs_max)
    print(out["matveev"], flush=True)
    per_small = []; Hall = max(Hs_max, 2)
    for s in words(Hs_max):
        K = strip(T_of(s))
        if K == 1: continue
        alphas = (19, 2, 5, K); logs = [LOG[19], LOG[2], LOG[5], Decimal(K).ln()]
        cur = H0; rounds = 0
        while True:
            X0 = 3 * cur + 3; C = 10 ** (4 * len(str(X0)) + 8)
            lam, b1sq = reduce_once(X0, C, alphas, logs); tries = 0
            while (lam is None or lam <= 0) and tries < 4:
                C *= 10 ** 6; tries += 1; lam, b1sq = reduce_once(X0, C, alphas, logs)
            if lam is None or lam <= 0: break                     # keep the bound we have
            new = largest_ok(lambda H: Decimal(H + 1) * Decimal(10).ln() < -lam.ln(), 0, cur)   # Lambda < 10^-(H+1)
            rounds += 1
            if new >= cur or rounds > 6: break
            cur = new
        per_small.append(dict(word=name(s), K=K, final_H_bound=cur, rounds=rounds))
        Hall = max(Hall, cur)
    out["small_words"] = len(per_small); out["largest_final_bound"] = max(p["final_H_bound"] for p in per_small)
    out["worst_small_words"] = sorted(per_small, key=lambda p: -p["final_H_bound"])[:5]
    print("small words", len(per_small), "largest final bound", out["largest_final_bound"], round(time.time() - t0, 1), "s", flush=True)
    seen = {}; clashes = []; n = 0
    for w in words(Hall):
        n += 1; K = strip(T_of(w))
        if K in seen: clashes.append((name(seen[K]), name(w)))
        else: seen[K] = w
    out["exact_check"] = dict(H_up_to=Hall, words=n, clashes=clashes)
    out["seconds"] = round(time.time() - t0, 1)
    json.dump(out, open("pairs16.json", "w"), indent=1, default=str)
    print(out["exact_check"], out["seconds"])
