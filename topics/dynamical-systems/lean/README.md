# Dynamics.lean

- 31 theorems, core Lean 4.37.0-rc1 (no mathlib), no `sorry`.
- Compiles clean against the shared library build (/tmp/lib/build copy):
  `LEAN_PATH=/tmp/lib/build lean lean/Dynamics.lean`. `#print axioms` on iter_msq_lift,
  logistic_fold, cheb_joukowski, pflow_ruler: only propext, Classical.choice, Quot.sound.
- Imports FusedAlgebra.Base only.
- Intended module path: FusedAlgebra/Topics/Dynamics.lean.

Contents: iter; msq_lift, iter_msq_lift (1 + fⁿ(x) = (1+x)^(2ⁿ)), pow_mul_R, mpow_mpow,
dial_square_units, dial_power_units, parity_logistic (k = −2 is logistic r = 2), failure_square,
wall_superattracting, identity_repelling, iter_wall, iter_zero, periodic_iff, wall_fixed_mpow,
chain_circ, jets (msqJet, iter_jet_fst, iter_jet_ruler: slopes double in the log ruler),
flow_time_law, Chebyshev (cheb_joukowski, logistic_T2, cheb_comp); fields: fold_eq,
sum_is_neg_product, inv_inv, fold_inv, msq_inv, logistic_fold (4X(1−X) = X(w ∘₁ w)),
joukowski_circ, dial_wall_superattracting, pflow_add, pflow_ruler (parabolic model flow and
its exact Fatou ruler −2/x).
