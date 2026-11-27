# Referee report, round 22 (DRAFT)

*Saved by the round-22 thread from the referee's final message (the referee could not write this file itself); the
content is the referee's, lightly condensed in layout only. The referee worked from BRIEF.md alone and did not open
README section 25, explore21.*, checks21.* or figures21.py. Standard library only; peak memory 1.45 GB.*

| Item | Verdict |
|---|---|
| 1. Prop 21.1 | **Holds.** The hypothesis a >= 1 can be weakened to a >= 0. |
| 2. Conj 21.2 search | **Reproduced, no other example.** |
| 3. Lemma 21.3 | **Holds, with wording fixes.** |
| 4. Lemma 21.4 | **Holds.** Hand proof plus a check of every canonical word up to length 12. |
| 5. Thm 21.5 | **Holds.** The induction closes. It is still conditional on 21.2. |

## 0. Facts taken as proved, spot-checked (check_basics_log.txt)

- Theorem 14.1 and psi(xw) = psi(w) + d_x/s_w: no failures on all 87,381 words of length 0 to 8.
- Proposition 15.2: no failures on 3,000 random tails.
- Canonical form: keeps the fingerprint and is idempotent (all 349,525 words up to length 9); the rewriting mK -> lN,
  ml -> lm is order-independent (left sides cannot overlap; random orders on 20,000 words agreed); a word is canonical
  exactly when it starts with N, its second letter is not N, it has no mK and no ml, and it ends in N or K (lengths 1-10).
- Each term of psi(w) is at most 1/10 of the one before; no canonical word up to length 13 has its first two terms in
  ratio 1/19. phi(s) = v2(s) - 2 v5(s) - v19(s) equals the number of N/K letters (20,000 words).

## Item 1: holds

If v2(1/kappa) = a >= 0 then v2(kappa) = -a < 1 = v2(2), so v2(2 - kappa) = -a, v2(nu) = 0, v2(lambda) = -(a + b);
the valuations of 1/nu, 1/kappa, 1/mu, 1/lambda are 0, a, b, a + b. kappa = 1/10, mu = 1/20 gives 1, 2, 4, 8, and these
four values are s_N, s_K, s_m, s_l. Check: 138,283 random cases with a >= 0, no failures.
Fixes (wording): a >= 0 is enough (for a <= -1 the pattern fails: 60,524 random failures; for a <= -2,
v2(1/nu) = a + 1); say "2-adic valuation", not "factors of 2" (b can be 0 or negative, 1/lambda is not an integer);
kappa != 2 is needed, automatic when a >= 0.

## Item 2: the search reproduced (search22.py and logs)

Word-ratio test, by hand: a product of k >= 1 word ratios is 19^(L-p) / (10^q 200^L 20^M) with p + q = k; so
2^e2 5^e5 19^e19 qualifies exactly when, with A = -e2, B = -e5: 2B - A - e19 >= 1 and some integer L has
max(0, e19) <= L <= min(A - B, 2B - A). Cross-checked against brute force and 20,000 random triples.

| Per side | Bound | Numbers | Lists | Primitive balances found |
|---|---|---|---|---|
| 3 | 10^12 | 1,402 | 311,128 | only {1}={1} and {19,1}={20} |
| 4 | 10^12 | 1,402 | 1,365,589 | same two |
| 5 | 10^12 | 1,402 | 3,600,687 | same two |
| 3 | 10^15 | 2,600 | 1,656,668 | same two |
| 4 | 10^15 | 2,600 | 11,722,261 | same two (the full claimed search) |
| 6 | 10^10 | 853 | 679,682 | same two |

Controls (search22_control.log): an independent slow search (3 a side, to 10^7) finds the same two; loosening to "each
ratio <= 1/10" finds 14,448 primitive balances (e.g. {76,5}={80,1}, {125,4}={128,1}), 14,446 not word-shaped. Evidence
only; 21.2 is still open.

## Item 3: holds, with wording fixes

Proof as stated, with: take P a balancing subfamily containing the largest number x that is minimal by inclusion; P is
primitive (a proper balancing Q inside P cannot contain x, and then P \ Q is a smaller balancing part containing x).
Case {20c} vs {19c, c}: 19c must be S tau'_0 (else S tau'_0 >= 190c > x), c the next term (1/19 > 1/100), and a
neighbour ratio 1/19 forces L = M = 0, r = 19, i.e. W' starts N N, impossible for a canonical word.
Checks (check_lemmas_log_12.txt): 4,001 random balances, the minimal part containing the largest number always primitive.
The canonical hypothesis is needed: NN and Nm both have psi = 20, minimal part {19, 1} = {20}, leading terms 19 and 20 differ.
Wording: say "a minimal (by inclusion) balancing subfamily containing the largest number" (not unique) and state it is
primitive; W, W' must be nonempty for tau_0 to exist (psi = 0 forces both empty).

## Item 4: holds

Hand proof: u = p x v with p = l^i m^j, x in {N, K} (j >= 1 forces x = N); t_u = sigma (t_v + 1), sigma = s_p s_x,
phi(sigma) = 1. canon(u) = N c with c = v if v does not start with N; if v = N v1: v1 empty gives c empty,
v1 = l^a K v2 gives c = l^(a+1) N v2, v1 = l^a m^b N v2 gives c = l^a m^(b+1) N v2. So c has at most two preimages
(v = c, v = N v1), and t_(N v1) + 1 = (20/19)(t_c + 1). If v = v', sigma = sigma' and the valuations leave only
(i + 1, j - 1, N) as another candidate, which needs m^j K, ruled out. If v != v', sigma' = (20/19) sigma, so
phi(sigma') = 2 != 1. u empty: t_u = 0 < t_u'.
Exhaustive: all 443,520 canonical W up to length 12 grouped by (canon(u), t_u): 443,520 groups, none with two u
(canon(u) alone gives 129,905 classes, largest 140). Without canonicity, u = lN and u' = mK are the same map (9,199 such
groups up to length 8).

## Item 5: holds (conditional on 21.2)

The induction on (k, k') closes: canonical forms keep the fingerprint and never raise the N/K count; empty case by
psi = 0; W = N u, tau_0 = 1/s_u, Lemma 21.3 (21.2 needed only for sides of at most k and k' numbers) gives
S'/s_u = S/s_u', so S' psi(u) = S psi(u') and t_u = t_u'; induction gives canon(u) = canon(u'); Lemma 21.4 gives u = u'.
Checks: 15.3 directly on all 87,381 words up to length 8; all 443,520 canonical words up to length 12 have different
fingerprints; tau_0 = 1/s_u and psi(Nu) = tau_0 + psi(u) on 5,000 words.
Wording: say "one side at most k numbers and the other at most k'"; the searches cover bounded sizes only, so 21.5
reduces 15.3 to 21.2 and does not make 15.3 a theorem even for k, k' <= 4.
