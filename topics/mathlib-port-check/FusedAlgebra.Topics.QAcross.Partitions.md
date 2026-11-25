# FusedAlgebra.Topics.QAcross.Partitions

Theorems (13): SAME 13.

No WEAKER / DIFFERENT / MISSING / VACUOUS / low-confidence items. All truncation lengths (30, 51), Fin bounds (10, 7, 5, 13/3), the explicit 30-coefficient Euler list, the 51 partition numbers and the ten r4 prime values are identical. Changed definitions (`isSq` -> `IsSquare`, `sigmaStar` -> Finset sum over `Nat.divisors`, `pentCoeff` Bool -> Prop) compute the same values.

Note (not a statement issue): `triple_product_30` and `jacobi_four_squares_30` now make the kernel reduce `IsSquare` on ℕ and a `Finset` sum over `Nat.divisors`; whether `decide +kernel` still goes through is a build question, not a meaning question.

Route changes:
- `lagrange_below_30`: was `decide +kernel`; now a corollary of the new general `r4_pos`.

New port theorems with no original: 3 (`one_add_multAddList_neg_pow`, `r4_prime`, `r4_pos`). As the porter says, `r4_pos` is positivity of Jacobi's formula, not Lagrange's theorem (no link between `r4` and actual representation counts is proved, in either version).
