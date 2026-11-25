# Numerical.lean

- 32 theorems, core Lean 4.37.0-rc1 (no mathlib), no `sorry`, axioms only propext, Classical.choice, Quot.sound.
- Compiles clean against the shared library build: `LEAN_PATH=/tmp/lib/build lean lean/Numerical.lean`.
- Imports `FusedAlgebra.Base` only. Namespace `Numerical`, the library's scoped notation `◯` for `circ` (`open scoped CircleOp`).
- Intended module path: `FusedAlgebra/Topics/Numerical.lean`.

Main results: `relerr_mul`, `relerr_mul3`, `relerr_div` (relative errors compose by ∘₁ and ⊟),
`theta_is_fold`, `chain_eq`, `theta_bound` (Higham's Lemma 3.1 as a circle-operation fold), `gamma_lift`,
`gamma_bound` ((1+u)^n − 1 ≤ nu/(1 − nu), multiplied out), `lift_amplifier`, `inverse_sensitivity`,
`sterbenz_window_pos/neg` (the inverse is within a factor 2 of −a when 1/2 ≤ 1 + a ≤ 2),
`newton_residual`, `newton_iter_residual` (residual −(−e₀)^(2ⁿ)), `gold_eq` (Goldschmidt),
`hornerCirc_eq`, `horner_cube`, `mpow_horner`, `bernoulli`.
Floating point itself is not modelled; the float64 behaviour is checked in ../checks.py.
