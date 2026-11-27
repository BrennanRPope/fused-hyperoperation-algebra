# Referee brief, round 20: recheck Theorems 18.1 and 19.1 to 19.3 from their statements alone

You are an independent referee. Redo the proofs below from these definitions and statements, with your own code.
Do NOT open these files (they hold the proofs you are checking): README.md sections 22 and 23, pairs18.py, pairs18.json,
pairs19.py, pairs19.json, checks18*, checks19*, figures18.py, figures19.py, figures 56 to 63. You may reuse the earlier
referee's tools in referee18/ (they were written from definitions only) and anything in the Python standard library,
mpmath and python-flint (both installed). Install nothing.

## Definitions

Letters are maps y -> s y + t:  l = (19/200, 0), m = (1/20, 0), N = (1/19, 1/19), K = (1/10, 1/10).
A word is a composition of letters; its map is w(y) = s_w y + t_w; psi(w) = t_w / s_w.
Equivalently psi(wx) = psi(w)/s_x + d_x and psi(xw) = psi(w) + d_x/s_w, with d = 1 for N and K and 0 for l and m.
The fingerprint kappa(w) is psi(w) written as a reduced fraction with every factor 2, 5 and 19 removed from numerator
and denominator.

A reduced word with n+1 N/K letters is N x_1 ... x_n, each item x_k = l^L_k m^M_k y_k with y_k in {N, K}; K is allowed
only when M_k = 0 (because mK = lN as maps); and the first item is not a bare N (no word starts N N). Its gaps are
L_k + M_k.

Convenience claim (check it yourself against psi): with r_k = 19 for y_k = N and 10 for K, P_k = 200^L_k 20^M_k and
A_k = L_1 + ... + L_k, put T_0 = 1 and T_k = 19^(A_k) + r_k P_k T_(k-1). Then kappa(word) = T_n with all factors 2, 5,
19 removed. Write T_n = S K with S = 2^e2 5^e5 19^g.

## Results taken as already proved (and already independently rechecked)

Theorem 15.5: among words with at most two N/K letters, psi is a product of 2s, 5s, 19s only for N, K, NN, KN (with any
leading l, m). Theorem 16.1: no reduced word with three N/K letters has fingerprint 1. Theorem 16.2: different reduced
two-letter words have different fingerprints. Theorem 17.1: a reduced three-letter word never shares a fingerprint with a
reduced two-letter word.

## Statements to recheck

1. Lemma F2 (claimed by hand): for every reduced three-letter word, v19(T_2) <= 2.
2. Theorem 18.1 (claimed, computer-assisted): two different reduced three-letter words never share a fingerprint.
   Claimed final bound on every gap after lattice rounds: 44, then an exact comparison of 226,965 words with no clash.
3. Lemma F3 (claimed by hand plus finite casework): for every reduced four-letter word, v19(T_3) <= 3.
4. Claimed (by hand): for a reduced four-letter word, 2^e2 5^e5 <= 80, reached only by N K m N N; so S <= 80 * 19^3.
5. Lemma N (claimed): for two different reduced words with at most four N/K letters and the same fingerprint, no linear
   form in logarithms used in the proof vanishes. (Check whichever vanishing cases your own proof needs to exclude.)
6. Theorem 19.1 (claimed, computer-assisted): no reduced four-letter word has fingerprint 1, and none shares a fingerprint
   with a reduced two-letter word. Claimed: Matveev start about 2.1e29; rounds 369, 61, 50, 50.
7. Theorem 19.2 (claimed): no reduced four-letter word shares a fingerprint with a reduced three-letter word.
   Claimed: start 2.1e62; rounds 4.0e18, 269, 68, 59, 55, 55.
8. Theorem 19.3 (claimed): no two different reduced four-letter words share a fingerprint.
   Claimed: start 4.8e112; rounds 4.0e40, 1.0e23, 371, 74, 66, 61, 59, 59; exact comparison of 106,624,096 four-letter
   words (gap sum <= 59) against the short boxes, no clash.
9. Claimed (by hand) for five letters: for W = N K l K K m^M N, T = 19 + 19 * 20^M * T(N K l K K), T(N K l K K) = -1
   modulo 361, so v19(T) = 2 + v19(M - mu) for a 19-adic mu, unbounded (3, 4, 5 at M = 0, 38, 5092).

## What to do

Build your own proofs (Matveev's theorem in the form of Bugeaud, Mignotte and Siksek 2006, Theorem 9.4; LLL with an exact
check of reducedness; de Weger-style reduction; exact 2-adic or 5-adic logarithms; any case split you like). Your bounds
need not match ours; what matters is whether each statement holds. Priorities, in order: items 1, 3, 4, 9 (cheap);
item 2; items 6 and 7; item 8. If item 8 is too heavy to finish, redo its Matveev start and as many rounds as you can,
and an exact check up to whatever bound you reach, and say plainly how far you got.
Keep peak memory under 4 GB. Work only in referee20/. Put scripts and logs there, and a REPORT.md with: for each item,
holds / fails / not finished; your Matveev start, your rounds, your exact check, next to ours; any gap in the statements
or in what they would need (for example a vanishing case or a missing hypothesis).
