# Numerical analysis through the circle operation (Nov 25, 2026)

## Files
- index.html: the page (self-contained; figure data from page_data.json inlined).
- checks.py -> checks_output.txt: 41 checks on real float64 (math.fma, log1p, expm1) against exact
  Fractions and 60-digit mpmath. All pass. Also writes page_data.json.
- page_data.json: data behind the page's four charts.
- lean/Numerical.lean: 32 theorems, compiles clean (see lean/README.md).

## Findings
1. The lift (1+a)(1+b) − 1 loses about log10(1/|a|) digits in float64 (4 of 16 left at |a| = 1e-12);
   a + b + ab and fma(a, b, a + b) keep all 16. Mechanism: the lift is a translation that carries the
   identity 0 onto the coarse float grid near 1; the way back amplifies a rounding by (1 + a)/a. Known
   (catastrophic cancellation; Kahan, Goldberg 1991).
2. Near the inverse (a ✳ b ≈ 0) plain a + b + ab loses up to ten digits, but fma(a, b, a + b) is
   correctly rounded whenever 1/2 ≤ 1 + a ≤ 2, because the inverse then lies in Sterbenz's window
   and a + b is exact (248/248; outside the octave 572/952 were not). Ours as stated; a direct
   consequence of Sterbenz 1974 plus single-rounding fma.
3. log1p/expm1 are the rung chart (E at t = ±1) that turns ✳ into + while fixing 0 with slope 1.
   One routine covers every dial: x1 ✳k … ✳k xn = expm1(Σ log1p(k x))/k. Compound interest each
   second: naive keeps 7 digits, log1p/expm1 all 16; "at least one of 100,000 chances of 1e-12":
   naive off by 1.6e-7, dial fold exact. Known trick (Goldberg 1991), packaged per dial value.
4. Relative errors compose by ✳ (products) and ⊟ (quotients); Higham's θn is exactly the ✳-fold of
   the step errors (checked to the bit), bounded by (1+u)^n − 1 ≤ nu/(1 − nu) (proved). Known content.
5. Newton's division-free inverse is v ← v ✳ (−(u ✳ v)); the residual goes e ↦ −e², the corner of
   (1+e)(1−e). Goldschmidt: (−e) ✳ e ✳ e² ✳ … = −e^(2^n). Horner: a chain of circle operations sums the
   2^n − 1 cube corners in n fma steps. Classical, restated.

## Verdict
A small practical lead, not a research one: a tested two-function kit (multadd via fma, dialfold via
log1p/expm1, with the wall handled) would give every digit to any code in the project that compounds
rates, combines small probabilities or folds the cube. It predicts nothing new mathematically.
