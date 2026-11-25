# Functional analysis through the circle operation (Nov 25, 2026)

Finding: for operators A ∘₁ B = (I + A)(I + B) − I; the way back exists exactly when −1 is not in the spectrum, and
det(I + K) − 1 is the circle-operation fold of the eigenvalues (the string's determinant is the sine product, zero each time an
eigenvalue crosses −1). At k = −1 the same fold is the chance of seeing at least one point of a determinantal process.

## Files
- index.html: the page (self-contained, all figures drawn as SVG by script; two sliders).
- checks.py -> checks_output.txt: 48 of 48 checks pass (numpy, scipy, sympy; Nyström quadrature; 40000-sample CUE Monte Carlo).
  results.json: Volterra norm roots, Neumann misses, CUE kernel eigenvalues and count distributions, void chances.
- lean/FunctionalAnalysis.lean: 35 theorems, compiles clean, no sorry, standard axioms; operator laws over noncommutative rings.
  See lean/README.md.

## Findings
1. Perlis (1942) quasi-regularity r + s + rs = 0 is exactly circle-operation invertibility; Jacobson's circle a + b − ab is the k = −1
   dial. Quasi-inverse (I + A)^-1 − I, reverse-order rule, resolvent identity A' − B' = (I + A')(B − A)(I + B') (Hilbert 1904).
   Known.
2. Neumann partial sums miss by exactly one term: A ∘₁ s_n = −(−A)^{n+1} (Lean, any ring), so nilpotents have finite
   quasi-inverses. Volterra V: ||V^n|| ≤ 1/n!, spectral radius 0, every λV stays off the wall; quasi-inverse kernel
   −λ e^{−λ(x−t)}. Known (Neumann 1877, Volterra).
3. det and trace are blind to order: det(I + A ∘₁ B) = det(I+A)det(I+B) for noncommuting A, B. Two nilpotents N, tN' combine by the circle operation ∘₁ (called mult-add before 1.0) to
   eigenvalues riding |1 + z| = 1 (the complex page's circle), a sixth circle root of 0 at t = 1, never reaching the wall.
   Small worked example, ours; the facts are classical.
4. Fredholm: string Green kernel min(x,y) − xy, det(I − z²G) = sin z / z (checked by Nyström); zeros where an eigenvalue of −z²G
   hits −1; sign flips per crossing; cube layers e1 = 1/6, e2 = 1/120. Lidskii, Fredholm. Known.
5. Determinantal points: P(at least one point in a region) = λ1 ∘₋₁ λ2 ∘₋₁ ... = 1 − det(I − K); counts = independent coin flips
   (HKPV 2006), checked against CUE(8) Monte Carlo. Shirai–Takahashi α-processes: void = ∏(1 + αλ)^(−1/α) = (1 + kF_k)^(−1/k):
   the dial runs fermions (k=−1) → Poisson (k=0) → bosons (k=+1). Known math; reading α as the dial is ours.

## Verdict
Yes, one lead: point processes on the dial. Everything else is classical functional analysis seen through the lift. Determinantal
and permanental processes use the whole dial with physical meaning at each named value and give checkable statistics; open
questions: what the landscape bend d and the fractional rungs mean as point processes.

Links used: linear algebra 6MmBPSMK179Q7ZRFTjz3L3. The complex-analysis page (same session, not yet published) is not linked; link it once it has an id.
