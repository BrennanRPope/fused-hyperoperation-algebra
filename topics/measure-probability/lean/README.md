# Lean: MeasureProbability.lean

- 36 theorems, core Lean 4.37.0-rc1, no mathlib, no `sorry`.
- Compiles clean: `LEAN_PATH=/tmp/lib/build lean lean/MeasureProbability.lean` (no errors, no warnings).
- Axioms checked with `#print axioms` on lam_corners, dual_weight, no_norm_past_wall, sure_atom_exists, odds_or, indep_dial, dial_gap_cov, yor, lam_equal_atoms, norm_two_sign: only propext, Classical.choice, Quot.sound.
- Imports `FusedAlgebra.Base` and `FusedAlgebra.Topics.Statistics` (reuses its finite expectation `Ex`).
- Intended module path: `FusedAlgebra/Topics/MeasureProbability.lean`.

Parts: (1) Sugeno λ-measure as the dial fold `lam`: disjoint union, lift, scale, Möbius corners, 3-atom expansion, λ = 0 additive, λ = −1 union; (2) two-atom normalization, its sign, multiplicative inclusion–exclusion, duality with weight `multAddInv λ`, involution, q* = 1/q, equal atoms = g·[m]_q; (3) the wall: `wall_blocks`, `wall_list`, `no_norm_past_wall`, sure atom at λ = −1; (4) odds of an "or" and odds against an "and" mult-add; (5) independent product distributions, `indep_dial` (E[X ✳ₖ Y] = EX ✳ₖ EY), `dial_gap_cov` (gap = k·Cov), `martingale_step`; (6) Yor's formula, discrete.
