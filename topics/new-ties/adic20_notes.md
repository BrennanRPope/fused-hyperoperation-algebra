# Round 20 (DRAFT): notes on the 19-adic step for five letters

Files: adic20.py (`python3 adic20.py all` writes adic20.json, about 4 minutes, peak about 13 MB; `python3 adic20.py check`
recomputes the key claims in under a second), adic20.json.

## 1. All 12 cancelling patterns have one shape (by hand, checked)

The five-letter casework (four items, as in checks19.py Q1) has 3,838 tie patterns; 12 cancel modulo 19. In all 12 the tied
terms are j = 0, 1, 4, the last item is a bare N with L4 = 0 (it may carry m^M), and L1 + L2 + L3 = A3 equals the number of
N letters, so the L's and letters of the front are fixed by the pattern and the only free gaps are M = M4 and the m-counts
of the N-items among x1, x2, x3 (0 to 3 of them). Writing W = N x1 x2 x3 m^M N,

  T(W) = 19^A3 + 19 · 20^M · T3,  v19(T3) = A3 - 1,  u := T3 / 19^(A3 - 1) = -1 mod 19,

so T(W) = 19^A3 (1 + 20^M u) exactly, and since 20^x = 1 + 19 Z_19 with v19(20^x - 1) = 1 + v19(x),

  **g(W) = A3 + 1 + v19(M - mu(u)),  where mu(u) in Z_19 has 20^mu = -1/u.**

The linear form is Lambda = 20^M · (-u) - 1, i.e. M log_19(20) + log_19(-u), two 19-adic logarithms with one integer
exponent; u depends only on the front. (The README family N K l K K m^M N is the pattern with no free front gap:
A3 = 1, u = 220,209 = -1 mod 361, so mu = 0 mod 19 and g = 2 + v19(M - mu).)

| Pattern (front, then m^M N) | A3 | free front gaps | u at zero front gaps |
|---|---|---|---|
| N lN llN lN | 4 | 3 | 1,608,003,819 |
| N lN lllN N | 4 | 3 | 1,608,000,380 |
| N lN llN K | 3 | 2 | 80,400,209 |
| N lN lK lN | 3 | 2 | 80,403,819 |
| N lN llK N | 3 | 2 | 80,400,380 |
| N lN lK K | 2 | 1 | 4,020,209 |
| N K llN lN | 3 | 2 | 88,003,819 |
| N K lllN N | 3 | 2 | 88,000,380 |
| N K llN K | 2 | 1 | 4,400,209 |
| N K lK lN | 2 | 1 | 4,403,819 |
| N K llK N | 2 | 1 | 4,400,380 |
| N K lK K | 1 | 0 | 220,209 |

Every other five-letter word has g = the smallest term size <= 4 (the casework). Checked: the identity on 312 random words
(exact T); u has a prime factor other than 2 and 5 for every front used (so 20 and -u are multiplicatively independent).

## 2. The two-logarithm 19-adic bound (published; constants from restatements)

Bugeaud and Laurent (1996, J. Number Theory 61, 311-342), as restated in arXiv 2401.06555 (Lemma 2.2, "Corollary 1 of
[BL]") and arXiv 2606.18500 (Theorem 2.3, "Corollary 2 of [BL]"), which agree: for multiplicatively independent eta1, eta2
with v_p(eta_i) = 0 and positive integers b1, b2,

  v_p(eta1^b1 eta2^b2 - 1) <= 24 p g / ((p - 1)(log p)^4) · E^2 · D^4 · h'(eta1) h'(eta2),

with h'(eta) >= max(h(eta), log p / D), E' = b1/h'(eta2) + b2/h'(eta1), E = max(log E' + log log p + 0.4, 10, 10 log p),
g the least integer with v_p(eta_i^g - 1) > 0 for both. **We did not see the 1996 paper itself.** Some citations give
p^f - 1 in place of g; that reading multiplies every bound below by 18.

Here p = 19, eta1 = 20, b1 = M, eta2 = -u, b2 = 1, D = 1, g = 1: the constant is 0.337, and 10 log 19 = 29.44 is the
largest term of E until M > about 1.4 × 10^12 · log|u|. So

  g(W) <= A3 + 0.337 · E^2 · log 20 · log|u|  ~  A3 + 875 · log|u|  (M up to 10^13),  about 52,000 · log|u| at M = 10^100.

**Linear in the height of the front** (log|u| grows like 3 + (free front gaps) · log 20), **polylogarithmic in M.**
Family N K l K K: g <= 10,771 for every M <= 10^13, <= 652,696 for M <= 10^100. Against the size step this is useless as it
stands (it lets S grow like 19^(c · height)); its use is as the a priori input below.

## 3. The 19-adic reduction for a fixed front (exact, computed)

For a fixed front, mu(u) mod 19^k is computed exactly, digit by digit (20^(19^j) = 1 + c_j 19^(j+1), c_j a unit), and
checked by pow. Then **for 0 <= M <= B the largest v19(M - mu) is the largest k with (mu mod 19^k) <= B, reached at
M = mu mod 19^k** (and no larger value is possible: a larger one would make mu mod 19^(k+1) <= B). Checked against brute
force for the family up to M = 10^5 (largest g = 6 at M = 25,669) and against g = 3, 4, 5 at M = 0, 38, 5092. The family's
records: g >= 4, 5, 6, 7, 8 first at M = 38; 5,092; 25,669; 1,459,200; 43,552,883.

All 12 patterns, every front with free front gaps adding to at most 30 (13,517 fronts):

| B | exact max g over all M <= B | where | log_19 B | exceedance over ceil(log_19 B): count of fronts |
|---|---|---|---|---|
| 10^2 | 9 | N lm^4N llm^11N lm^12N m^43 N (A3 = 4) | 1.56 | -1: 9,715; 0: 3,612; 1: 187; 2: 3 |
| 10^100 | 87 | N K lllm^23N m^5N m^M N (A3 = 3) | 78.2 | -1: 12,207; 0: 1,230; 1: 75; 2: 4; 4: 1 |
| 10^200 | 164 | (A3 = 4) | 156.4 | -1: 11,198; 0: 2,196; 1: 116; 2: 7 |
| 10^300 | 242 | (A3 = 4) | 234.6 | -1: 9,309; 0: 4,024; 1: 178; 2: 6 |

So in practice g <= A3 + 1 + log_19 B + (at most 4 here), as for a random 19-adic mu.

## 4. Feasibility for all five-letter words (judgement)

- The 19-adic layer is feasible and cheap. Every unbounded case is a two-logarithm form with b2 = 1, so Bugeaud-Laurent
  suffices; Yu's many-logarithm theorem is not needed.
- A priori bounds: g enters the front-step forms only as the coefficient of log 19 (and the back step only through |m|),
  so g <= C · X · (log X)^2 with C about 5 × 10^4 (× 18 conservatively) replaces 9X + 20; Matveev's ell = 1 + log(coeff)
  grows by about 11 + 2 log log X, so X* stays finite. My rough guess for 5 against 5 is 10^150 to 10^250 (extrapolated
  from 2.1e29, 2.1e62, 4.8e112 for 4v2 to 4v4; not computed).
- Reduction: once a round's front steps fix the front, the digit argument gives the exact g for M up to the round's bound,
  costing milliseconds per front, so S(other) <= 2^e2 5^e5 19^(A3 + 1 + log_19 B + small): the size step loses only about
  log B, and closes as before.
- What remains is the five-letter version of everything in round 19: the 2- and 5-part bound for five letters (by hand,
  not done), the trees 5v1 to 5v5 with a third middle gap (more exact back steps over two- and three-item back parts), and
  the exact check. Counting words: four items with gap sum <= 59, 65, 70 is 1.0 × 10^10, 2.0 × 10^10, 3.4 × 10^10, against
  1.07 × 10^8 three-item words at 59 in round 19, so the exact check alone is 100 to 300 times round 19's (hours in this
  Python, and a short-box table that may need care with memory). Large, but of the same kind.

## 5. What can be stated now

**Proposition 20.1 (by hand).** Every reduced five-letter word whose 19-part is not the smallest term size lies in one of
the 12 patterns above, and for those g(W) = A3 + 1 + v19(M - mu(u)), with u = T3 / 19^(A3 - 1) = -1 mod 19.

**Proposition 20.2 (published bound applied).** For every such word, g(W) <= A3 + 0.337 · E^2 · log 20 · max(log|u|, log 19)
(E as above, b1 = M); × 18 under the conservative reading.

**Proposition 20.3 (computed, exact).** For every word of the 12 patterns whose free front gaps add up to at most 30 and
every M <= 10^300, g(W) <= 242; for M <= 10^100, g <= 87; for M <= 100, g <= 9; each maximum is reached (table above). For
N K l K K m^M N alone: g <= 80, 159, 236 for M <= 10^100, 10^200, 10^300.

Who: Claude (round 20, a worker of the thread) found the 12-pattern shape, the identity, the digit argument and the
computations; Bugeaud-Laurent is published (constants from the two restatements named above). Not checked by a person.
