# KnotTheory.lean

- 45 theorems, core Lean 4.37.0-rc1 (no mathlib), no `sorry`.
- Compiles clean against the shared library build (/tmp/lib/build copy):
  `LEAN_PATH=/tmp/lib/build lean lean/KnotTheory.lean`. `#print axioms` on tangle_figEight,
  twist_formula, conway_stack, qnum_mul: only propext, Classical.choice, Quot.sound.
- Imports FusedAlgebra.Base, FusedAlgebra.Topics.QAcross.QDial, FusedAlgebra.Topics.QAcross.Knots.
- Intended module path: FusedAlgebra/Topics/KnotTheory.lean (after QAcross.Knots in build order).

Contents: tangle algebra (tsum, vsum, closures) on any `Lean.Grind.CommRing`; the lift law
closeD_tsum (connected sum) and walls; twist_formula ([n] = Aⁿ[0] + Aⁿ⁻²[n]_q[∞], q = −A⁻⁴);
walls_merge / loop_zero / q_one (A⁴ = −1); potts_Q; qnum_mul ([pq]_t = [p]_t[q]_(t^p)),
alexander_T2; field section: ratio_dial (sum is the d-dial in g/f), lift_is_closing, stack_dial,
twist_step_y (y ↦ qy + 1), stack_step_y, twists_on_circle, twists_close, conway_twist,
conway_add, conway_stack; integer continued fractions (rational_fractions, SL(2,ℤ) steps);
Laurent-bracket kernel computations showing the tangle algebra reproduces the PD-code brackets
of QAcross.Knots (trefoil, figure-eight twice, cinquefoil).
