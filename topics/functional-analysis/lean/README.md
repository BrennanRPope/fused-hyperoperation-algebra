# FunctionalAnalysis.lean

- 35 theorems (including one record-extensionality helper and two power lemmas), core Lean 4.37.0-rc1, no mathlib, no `sorry`.
- Compile status: compiles clean against the shared library build (`LEAN_PATH=/tmp/lib/build`, imports `FusedAlgebra.Base` only).
  `#print axioms` on resolvent_identity, nilpotent_quasi_inverse, gap_fold, op_dial_lift: propext, Classical.choice, Quot.sound
  only (nilpotent_pair_sixth_root: propext only). The #print lines were removed afterwards.
- Intended module path: `FusedAlgebra/Topics/FunctionalAnalysis.lean` (namespace `FunctionalAnalysis`). Add after Base in build.sh.
- Operator laws are proved over a general noncommutative `[Lean.Grind.Ring R]` with the operator circle operation ∘₁ (called mult-add before 1.0) `opCirc a b = a + b + a*b`
  (local notation ◯ₒ, because the library's ◯ needs a commutative semiring); `opCirc_comm_eq` shows it is ◯ when R is commutative.
  Note: in a noncommutative ring `grind` does not use equational hypotheses, so those proofs rewrite with the hypothesis first
  and call `grind` only on pure identities.

Sections: lift and wall (op_lift, op_assoc, wall_left/right, wall_no_quasi_inverse); quasi-inverses (perlis_iff,
quasi_inverse_unique, qi_of_inverse, jacobson_circle, qi_reverse, resolvent_identity); Neumann series (neumann_miss,
neumann_miss_right, nilpotent_quasi_inverse, square_zero_qi, cube_zero_qi); commuting, conjugation, characters, dial
(op_comm_of_commute, op_conj, character_op, op_dial_lift); noncommuting nilpotent pair on 2×2 records (tr_circ_swap,
det_circ_swap, nilpotent_pair, nilpotent_pair_noncomm, nilpotent_pair_sixth_root); gap probabilities (foldW_lift, gap_fold,
fredholm_fold).
Not in Lean (needs analysis): norms, spectra of infinite operators, Volterra, Fredholm determinants of integral operators,
random matrices. Those are checked in checks.py.
