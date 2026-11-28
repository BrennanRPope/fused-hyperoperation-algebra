# Calculus from the circle operation (Nov 25, 2026)

> **Names:** some words or symbols on this page (⊕, mult-add) also mean something else in standard math, physics or computing, or inside this project. [NAMING-COLLISIONS.md](../../NAMING-COLLISIONS.md) says which meaning is in play.

Page: https://claude.ai/artifact/VBRCmFrQVwQiXjMxZkY5SY (source index.html here)

- lean/FusedCalculus.lean: 29 theorems, core Lean 4.37.0-rc1, no mathlib, no sorry, standard axioms only
  (checked against the library's Base and Cube). Namespace MultAdd.Calculus.
  For the library: copy to FusedAlgebra/Topics/Calculus.lean (imports FusedAlgebra.Base and FusedAlgebra.Cube),
  add FusedAlgebra.Topics.Calculus to build.sh and FusedAlgebra.lean. Defines its own mpow_mul (MultAdd.Calculus.mpow_mul).
- checks.py: 12 symbolic/numeric checks of the continuous claims (limits), all pass.

Findings: discrete product rule = corner square; relative rates over step h combine by the circle operation ∘₁ (called mult-add before 1.0) with
weight h = time-scale calculus (Hilger 1988, Bohner-Peterson 2001: p ⊕ q = p+q+μpq). Output-lifted derivative
f'/(1+f) = geometric calculus (Grossman-Katz 1972), input-lifted = Jackson q-derivative (q = 1+h), both = bigeometric.
Fundamental theorem = Volterra product integral. Newton forward formula (cube layers) -> Taylor in the limit.
