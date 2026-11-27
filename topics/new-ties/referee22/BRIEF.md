# Referee brief, round 22 (DRAFT): recheck four statements from their statements alone

You are an independent referee. Check each statement below with your own reasoning and your own code. Do NOT open
README.md section 25, explore21.py, explore21.json, checks21.py, checks21_output.txt, checks21.json or figures21.py in
/mnt/project-files/math-theory/new-ties/ (they hold the original proofs and checks). Everything else may be read.
Write your scripts and logs in /mnt/project-files/math-theory/new-ties/referee22/ and your verdict in REPORT.md there.
Python 3 with fractions is enough; do not install packages. Keep any script's peak memory under 4 GB.

## Definitions

Letters are maps y -> s y + t on the line:  l = (19/200, 0),  m = (1/20, 0),  N = (1/19, 1/19),  K = (1/10, 1/10).
A word is a string of letters; as a map, xw means x o w (apply w first). For a word w write w(y) = s_w y + t_w and
psi(w) = t_w / s_w. The fingerprint of w is psi(w) with every factor 2, 5 and 19 removed from numerator and denominator.

Taken as proved (you may use them, and may spot-check them):
- Theorem 14.1: psi(w) is the sum, over the positions i of N and K in w, of the product of 1/s_x over the letters x after i.
  Also psi(xw) = psi(w) + d_x / s_w, with d = 1 for N and K and 0 for l and m.
- Proposition 15.2: psi(Kw) = psi(Nw); psi(NNw) = psi(KNw) = psi(Nmw) = psi(Kmw); mK = lN and ml = lm as maps; leading
  l and m do not change psi; a trailing l multiplies psi by 200/19 and a trailing m by 20.
- Canonical form of a word: drop leading l and m; make the first letter N; if the second letter is N, make it m; rewrite
  mK -> lN and ml -> lm until nothing changes; drop trailing l and m. It keeps the fingerprint. A word is canonical if it
  equals its canonical form.
- Conjecture 15.3 (open): two words with the same fingerprint have the same canonical form.

Collision surface (for item 1): multipliers kappa, mu > 0 rational, nu = kappa/(2 - kappa), lambda = mu (2 - kappa).

A word ratio is 19^L / (r 200^L 20^M) with integers L, M >= 0 and r in {19, 10}. A list of positive integers is
word-shaped if, sorted largest first, every ratio (next / previous) is a product of one or more word ratios. A balance is
two lists of products of powers of 2, 5 and 19 with equal sums; it is primitive if no nonempty proper subfamily of all its
numbers (taken from both sides) balances, i.e. sum of the chosen left numbers = sum of the chosen right numbers.

## Statements to check

1. Proposition 21.1. If 1/kappa has a >= 1 factors of 2 and 1/mu has b factors of 2 (2-adic sizes), then 1/nu, 1/kappa,
   1/mu, 1/lambda have 0, a, b, a + b factors of 2. With kappa = 1/10, mu = 1/20 this gives 1, 2, 4, 8.
2. Conjecture 21.2 (open; check the search only). In a primitive balance whose two sides are each word-shaped, either each
   side is one number (the same) or the sides are {19c, c} and {20c}. Claimed search: no other example with up to 4 numbers
   a side and all numbers (after dividing by their gcd) at most 10^15. Redo a smaller search yourself (for example 3 a side,
   up to 10^12) with your own code.
3. Lemma 21.3. Let W, W' be canonical words with S' psi(W) = S psi(W') for products S, S' of 2, 5 and 19. If Conjecture
   21.2 holds for the smallest balancing part (of the balance S' (terms of psi(W)) = S (terms of psi(W')), cleared of
   denominators) that contains the largest term, then the leading terms agree: S' tau_0 = S tau'_0, where tau_0 is the
   largest term of psi(W).
4. Lemma 21.4. Let W = N u and W' = N u' be canonical words with t_u = t_(u') and the same canonical form of u and u'.
   Then u = u'.
5. Theorem 21.5. If Conjecture 21.2 holds, Conjecture 15.3 holds for all words; more precisely, 21.2 for lists of at most
   k and k' numbers gives 15.3 for every pair of words with at most k and k' letters N or K.

For each item say: holds / holds with a fix (give the fix) / fails (give a counterexample or the broken step). Look hard
for gaps: hidden assumptions about canonical words, the order of letters, empty words, trailing letters, the case where a
word has one N/K letter, and whether the induction in 5 really closes. Prove or disprove by hand, and back each with an
exhaustive or random check where possible (for example on all canonical words up to 10 or 11 letters).
