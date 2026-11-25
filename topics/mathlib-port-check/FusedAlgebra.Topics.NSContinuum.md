# FusedAlgebra.Topics.NSContinuum: review

Theorems (18): SAME 10, GENERAL 6, EQUIV 1, DROPPED-OK 1. No WEAKER / DIFFERENT / MISSING / VACUOUS. No low-confidence items.

- `mpow_one_nat`, DROPPED-OK to `MultAdd.Cube.mpow_one_add_one`: `mpow 1 n + 1 = 2 ^ n`, with numerals defaulting to ℕ. Same statement.
- `annulus_graded` (EQUIV): bsum becomes a range sum, layer becomes Nat.choose, and the factor 1^(k+1) is dropped.
- onRay and its laws went from ℤ to a comm ring; onRay_trans needs IsDomain, which holds for ℤ. The decide checks use the same numbers (3,25,2,20), (3,20,2,16), ¬(3,20,3,25).
- `planar_stretching_zero`: the pointwise ∀ i form is now the function equation `G *ᵥ Pi.single 2 w = 0`. Pi.single 2 w is `if k = 2 then w else 0`, and mulVec has the same G i j convention.

Route changes worth knowing:
- `annulus_graded`: the original used the user's own binomial theorem (CubeProofs.binomial with bsum). The port uses `Cube.mpow_eq_sum_choose`, i.e. mathlib `add_pow`.
- `mpow_one_nat`: the user's induction is replaced by Cube's one_add_mpow.

New port theorems with no original: 0.
