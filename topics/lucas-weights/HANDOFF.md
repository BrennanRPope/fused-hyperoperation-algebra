# Lucas corner weights U_m(k, d): handoff (Nov 25, 2026)

Stopped partway at the coordinator's request: Brennan asked for the corner-weights work to be merged into one stream, led by the size/pair-weights task (corner weights as principal minors of one matrix). This file holds everything this thread found. Nothing here is pushed.

## Verdict so far

**The weights are known, and the corner formula is folklore-level: a two-line corollary of classical facts.** It is still worth keeping as a tool, but on its own it is not a research lead. It fits straight into the merged stream's picture: the Lucas weights are what you get when the "one matrix" is a 2×2 matrix N with trace k and determinant d.

## The statement (ours as a written formula, not as a discovery)

Plane law F(x, y) = (x + y + k·xy)/(1 − d·xy). U₀ = 0, U₁ = 1, U_m = k·U_(m−1) − d·U_(m−2). For a chain x₁ ⊕ … ⊕ xₙ, with eₘ the elementary symmetric polynomials of the xᵢ:

  x₁ ⊕ … ⊕ xₙ = Σₘ Uₘ·eₘ / (1 − d·Σₘ U_(m−1)·eₘ)

A rank-m cube corner gets weight Uₘ on top and −d·U_(m−1) on the bottom.

## Why it's true (the whole proof)

1. Adding a in the law is the Möbius map with matrix I + a·N, where N = [[k, 1], [−d, 0]] (trace k, det d). All these matrices commute.
2. Equivalently, the law is the circle operation (called mult-add before 1.0) in the number system ε² = kε − d: (1 + aε)(1 + bε) = (1 − d·ab)·(1 + F(a, b)·ε).
3. So a chain is ∏(I + aᵢN) = Σ eₘ Nᵐ (Vieta: one term per cube corner).
4. Cayley–Hamilton: Nᵐ = Uₘ·N − d·U_(m−1)·I (classical; e.g. Greene, "Lucas sequences and traces of matrix products").
5. Reading off the slope (apply the matrix to 0) gives the formula.

With equal entries this reduces to the library's existing `Patterns.iter_lucas` (same recursion, binomial eₘ).

## Known vs ours

| Piece | Status | Source |
|---|---|---|
| Uₘ(k, d) themselves | Known: Lucas sequences (Lucas 1878); = Dickson polynomials of the 2nd kind E_(m−1)(k, d); = d^((m−1)/2)·Chebyshev U_(m−1)(k/2√d); bivariate Fibonacci polynomials | Yamagishi, Fib. Quart. 57(3) 2019, https://www.mathstat.dal.ca/FQ/Papers/57-3/yamagishi06142019.pdf |
| Nᵐ = UₘN − dU_(m−1)I | Known (Cayley–Hamilton) | Greene, https://www.d.umn.edu/~jgreene/papers/Lucas_Traces.pdf |
| tan / tanh of a sum of n terms as odd/even eₘ | Known | ProofWiki "Tangent of Sum of Series of Angles"; MathPages https://www.mathpages.com/home/kmath735/kmath735.htm |
| Classification of the plane laws by k² − 4d | Known (Coleman & McGuinness 1991, Pacific J. Math 147) | https://msp.org/pjm/1991/147-1/pjm-v147-n1-p02-s.pdf |
| ε² = kε − d number systems (complex / dual / double) | Known (Yaglom, *Complex Numbers in Geometry*); already in the project's why-i/cayley-klein page (curvature −Δ/4 = d − k²/4) | |
| General (k, d) corner formula, written out | Not found in print in searches on Nov 25 (5 searches: tangent generalizations, rational formal group n-fold sums, Fibonacci tangent analogues, Möbius/commuting matrices, hyperbolic Fibonacci functions). Call it "ours as stated, probably folklore". | |
| Possibly relevant, unread | Belbachir, "Inverse trigonometric and hyperbolic identities, series and products for Lucas sequences", arXiv 2610.24667 (Oct 2026). Abstract only: one two-term identity; full PDF not checked. | https://arxiv.org/abs/2610.24667 |

## Dial-plane map (worked out, figure not drawn yet)

Three regions, split by the parabola k² = 4d:
- **Inside (d > k²/4), elliptic:** weights oscillate like Chebyshev U (rotations). Layer n goes silent (Uₙ = 0) when d = 1, k = 2cos(πj/n): (0,1) n=2, (1,1) n=3, (√2,1) n=4, (golden φ,1) n=5 (matches `Patterns.lucas_five_golden`).
- **On the parabola, parabolic:** Uₘ = m·(k/2)^(m−1), like a derivative; the law is conjugate to plain addition.
- **Outside, hyperbolic:** weights grow like (larger root)^m; law conjugate to multiplication. All d < 0 is hyperbolic.
- **d = 0 line:** Uₘ = k^(m−1), the single-dial circle operation.

Named points, U₁..U₈ (checked, `checks_output.txt`):

| Point | Name | Weights |
|---|---|---|
| (0, 0) | plain addition | 1, 0, 0, 0, … |
| (1, 0) | circle operation | 1, 1, 1, 1, … |
| (−1, 0) | inclusion–exclusion | 1, −1, 1, −1, … |
| (0, 1) | tangent of a sum | 1, 0, −1, 0, 1, … |
| (0, −1) | tanh / velocity addition | 1, 0, 1, 0, … |
| (1, −1) | Fibonacci | 1, 1, 2, 3, 5, 8, 13, 21 |
| (2, −1) | Pell | 1, 2, 5, 12, 29, 70, 169, 408 |
| (1, −2) | Jacobsthal | 1, 1, 3, 5, 11, 21, 43, 85 |
| (3, 2) | Mersenne 2^m − 1 | 1, 3, 7, 15, 31, 63, 127, 255 |
| (2, 1) | parabola, Uₘ = m | 1, 2, 3, 4, … |
| (1, 1) | period 6 | 1, 1, 0, −1, −1, 0, 1, 1 |

## Checks

`checks.py` (sympy): 16 symbolic checks pass. Corner formula for chains of 1–5 fully symbolic in k, d, xᵢ; the pair recursion equals the law; Nᵐ = UₘN − dU_(m−1)I for m = 1..10; the named-point table.

## Lean: planned, not written

No Lean toolchain here (the elan installer was refused by the platform's safety check; nothing else tried). Plan that should compile in core Lean 4.37 with `grind`, reusing `Patterns.lucas k d` (same recursion, so no twin):
- `corners w l`: Σ wₘ eₘ(l), defined by `corners w [] = w 0`, `corners w (a :: l) = corners w l + a * corners (shift w) l`. Use named `shift` / `lin` / `smul` defs so rewrites stay syntactic.
- `corners_lin`, `corners_smul` (linearity in w, induction on l generalizing w).
- `chainPair`: (0, 1) then (A, B) ↦ (A + a(B + kA), B − d·a·A). Main theorem: chainPair = (corners U, corners V) with V 0 = 1, V(m+1) = −d·U m, using shift U = lin V k U and shift V = smul (−d) U pointwise.
- `chain_cross`: if x·B = A and y·(1 − d·x·a) = x + a + k·x·a then y·B' = A' (no division needed).
- `npow_lucas`: 2×2 power N^(m+1) = (U(m+2), U(m+1), −d·U(m+1), −d·U m).
- `plane_is_circ_eps`: (1, a)·(1, b) = (1 − d·ab, a + b + k·ab) in ε² = kε − d.
- Optional bridge: corners w l = Σ_(m ≤ |l|) wₘ·eₘ(l) with explicit `esym`.

## Open questions for the merged stream

1. In the principal-minors picture, is the Lucas case exactly "the matrix is a polynomial in one 2×2 N"? If so, the whole (k, d) plane is the rank-2 slice of the general matrix family.
2. Read Belbachir 2610.24667 in full to rule it out as a prior statement.
3. Figures still to draw: dial plane with regions and named points; the 3-cube with two numbers per corner (Uₘ top, −dU_(m−1) bottom).
