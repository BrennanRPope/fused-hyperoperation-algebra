# The power set with one dial per element (Nov 25, 2026)

Task: Brennan's "mini dials" idea from the voice call: each element gets its own dial for how much it takes part. Steered during the work to (1) merging symmetry, then complement symmetry, and whether they match; (2) the primes as an easier test bed (squarefree numbers, zeta); (3) other known sets, looking for common dials.

Page: published as an artifact (see PAGES-INDEX.md, "power set dials"); source index.html here.

## Files
- checks.py / checks_output.txt: 34 of 34 pass (fractions, brute-force OR-convolution, zeta/Mobius uniqueness, complement, Bernoulli simulation).
- primes_check.py / primes_output.txt: sieve to 10^7; zeta(s)/zeta(2s) at s=2,3; growth exponent of three halves of the primes.
- sets_check.py / sets_output.txt: naturals by size (Li2), positive rationals zeta(s)^2/zeta(2s), Gaussian prime ideals.
- lean/PowersetDials.lean: 48 theorems, merged and compiled in the shared library as Topics/PowersetDials (library now 2783 theorems, 110 modules).
- index.html: the page.

## Findings
1. Count: sum over subsets of prod of dials = prod(1+t) = 1 + (t1 ✳ ... ✳ tn).
2. Merging (OR-convolution, A∪B) combines dials elementwise by the circle operation ∘₁ (called mult-add before 1.0) (Lean conv_W). Respecting disjoint union forces one dial per element. Unique nonneg square root of P(S): dial sqrt2-1 per element. P^alpha = dial 2^alpha - 1, P^a join P^b = P^(a+b).
3. Halves form a box; partner(t) = (1-t)/(1+t), corners = splits P(A)xP(S\A), centre sqrt2-1. Angle t = tan(theta): halves add to 45 deg (halves only; thirds don't).
4. Complement A -> S\A turns t into 1/t, swaps union-merge and intersection-merge, fixes only t=1 (full P(S)). Union half sqrt2-1 vs intersection half sqrt2+1 (probabilities 0.293/0.707). No p is a half both ways.
5. Primes: prod(1+p^-s) = zeta(s)/zeta(2s) (s>1). Participation c_p with average alpha gives counts ~ X (log X)^(alpha-1) (Selberg-Delange). Corner (p=1 mod 4), plain 1/2, merging half all measure alpha = 1/2. Tiny dials linearise, so halves agree; participation complement fixes 1/2 = merging half.
6. Other sets: naturals by size (prod(1+q^n)): corner and merging halves agree (pi^2/24), plain 1/2 gives 54.5% (needs c=0.455). Integers: N is an exact half. Rationals: alpha=2, N and squarefree are complementary halves. Number fields: same as primes (Chebotarev). Reals: Poisson intensity halves exactly. Algebraic numbers: no canonical dial without degree+height.
7. No dial gives an iteration half step (f(0)=1 forced).

## Verdicts
- Half by joining: settled, known math (p-biased measures, OR-convolution, Poisson superposition). Not a lead.
- Primes/exponent alpha as a common dial: modest lead, a yardstick for any fractional power set; predicts nothing beyond Selberg-Delange yet.
- Iteration-sense half: dials can't reach it; still open without choice (previous page).
