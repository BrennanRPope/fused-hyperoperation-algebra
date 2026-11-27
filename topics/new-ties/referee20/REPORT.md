# Referee 20: Theorems 18.1 and 19.1 to 19.3, rechecked from the statements alone

Written by an independent Claude session (a subagent started in round 20) from BRIEF.md, without opening README sections 22
and 23, pairs18/19, checks18/19 or figures18/19. Saved here by the thread session, which could not let the subagent write
this file itself; the content is the subagent's report, unchanged except for layout.

**All nine items hold.** No counterexample, no false step, nothing left unfinished.

| # | Verdict | Referee's numbers | Claimed |
|---|---|---|---|
| 1 | holds | v19(T_2) <= 2, and 2 is reached | <= 2 |
| 2 (18.1) | holds | start 6.4e45; rounds 422, 42, 33; exact check of 80,291 words, no clash | every gap <= 44; 226,965 words |
| 3 | holds | v19(T_3) <= 3, and 3 is reached | <= 3 |
| 4 | holds | max 2^e2 5^e5 = 80, only N K m N N; e2 <= 4, e5 <= 1 | same |
| 5 | holds for every form the referee's proof uses | each vanishing case forces W = W' or is impossible | as stated |
| 6 (19.1) | holds | 4 vs 2: start 6.4e45; rounds 558, 50, 38 | start 2.1e29; rounds 369, 61, 50, 50 |
| 7 (19.2) | holds | start 2.8e62; rounds 765, 59, 47, 45 | start 2.1e62; rounds 4.0e18, 269, 68, 59, 55, 55 |
| 8 (19.3) | holds | start 1.55e79; rounds 1.57e19, 263, 61, 51, 49; exact check of 39,039,480 words, no clash | start 4.8e112; rounds 4.0e40, 1.0e23, 371, 74, 66, 61, 59, 59; 106,624,096 words |
| 9 | holds | first M with v19 = 3, 4, 5, 6: 0, 38, 5092, 25669 | 3, 4, 5 at M = 0, 38, 5092 |

The referee's bounds are on the gap sum of each word, not on each gap, so the word counts are not comparable. Extra check:
the exact comparison at the claimed bound 59 gives exactly 106,624,096 four-letter words; with all two- and three-letter
words to 60 included, no clash, and only the word N has fingerprint 1 (exact_59.log).

## Gaps and remarks

- Order of proof (item 5): the vanishing argument leans on the earlier theorems; 18.1 is used inside 19.3 and 17.1 inside
  19.2, so 18.1 has to be proved first.
- Non-reduced back parts (item 5): the argument also needs the rewrite N N z -> N (z + m) (and mK = lN). It keeps the
  fingerprint, but the statements never mention it.
- Fingerprint 1 in 19.1 needs no Baker theory: item 4 and F3 give T = S <= 548,720, which forces gap sum <= 2.
- Everything depends on the S bounds: items 1, 3 and 4 feed the size argument, the 2-adic cases and the exact check.

## Methods

- Items 1 and 3 (v19dp.py): the exact set of residues of T_n modulo 19^E over all L, M >= 0. It is finite because for an
  N item P runs over a whole coset 10^L <20>, and <20> is exactly the units = 1 mod 19; 200 is a primitive root mod 19^E,
  so K items with L >= E give every unit; 19^A matters only through min(A, E). 0 mod 19^(n+2) is never reached for n = 1,
  2, 3 items, and v19 = n is reached. Brute force on small words agrees.
- Item 4 (item4.py): an item whose multiplier is 0 mod 16000 resets T to 19^A (every A mod 800 tried); words of the 9 other
  items enumerated exactly. Four letters: (e2, e5) in {(0,0), (1,1), (2,1), (3,1), (4,1)}; three letters: max 40 (N m N N);
  two letters: no 2 or 5. So S <= 80 * 19^3, 40 * 19^2, 19.
- Item 9 (item9.py): T(NKlKK) = 220209 = -1 mod 361; mu computed 19-adically (v19(mu) = 1); v19(T) = 2 + v19(M - mu)
  verified for every M < 200,000.
- Items 2, 6, 7, 8 (engine.py, driver.py, exactcheck.py): each round a 3-log LLL step (min(H1, H1') <= X1); for each first
  item (43 to 13,693) a 4-log LLL step on log K(N x1) (min(H2, H1') <= X2); a 2-adic tree from the back using
  T = 19^(A_b) tau mod 2^(v2(R_b) + g2_b), bounding v2(19^m Q - 1) exactly for |m| <= Hs + 3; then the size step
  T(Q) <= T(P) SMAX. Each LLL output (python-flint) is checked: rows in the lattice, |det| unchanged, lambda1^2 >= min of
  the Gram-Schmidt norms. Matveev in the form of Bugeaud, Mignotte and Siksek (2006, Theorem 9.4).
- Lemma N (lemmaN.py): (a) the 3-log form vanishes only if the back-part fingerprints are equal, so the back parts are
  identical (rho = 1, same item, W = W') or of the NN type (rho = 20/19), which contradicts v2 - 2 v5 = v19; (b) the 4-log
  forms with log(K_a/K'_a) are nonzero since K_a != K'_a; (c) the 2-adic forms: vanishing makes the front fingerprints
  equal, so identical fronts, and rho = 1 with identical tails gives W = W'. Up to 1,635 vanishing cases per round were met;
  every one came out "equal" or "impossible".

Peak memory under about 3 GB. Scripts need mpmath (engine) and python-flint; they are not part of the repository's checks.
Logs: t181_rounds.log, t191_rounds.log, t192_rounds.log, t193_rounds.log, t181_exact.log, exact_all.log, exact_59.log,
v19dp.log, item4.log, item9.log, lemmaN.log, tests.log, test_exact.log. Figure: fig_rounds.png.
