# FusedAlgebra.Topics.QAcross.Knots

Theorems (19): SAME 18, GENERAL 1.

No WEAKER / DIFFERENT / MISSING / VACUOUS / low-confidence items.

- All 14 kernel checks: PD codes, n, m, writhes and every coefficient list / pair compared mechanically and by eye; identical. Still `decide +kernel`, no native_decide.
- Changed definitions (pzeros, trim, numA/numB, evalOne, negOnePow -> (-1)^n) are the same functions.
- `state_layers`: GENERAL (all n, k instead of n, k < 6; `Nat.choose` instead of `Cube.layer`, which has the identical Pascal recursion).

Route changes:
- `state_layers`: finite decide became an induction proof (strengthening).
- Local rules (`loop_value_q`, Reidemeister I/II, `writhe_cancels`): `grind` became `linear_combination`; same identities.

New port theorems with no original: 10 (toPolynomial_* bridge lemmas, length_allStates).
