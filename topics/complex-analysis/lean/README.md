# ComplexAnalysis.lean

> **Names:** some words or symbols on this page (wall) also mean something else in standard math, physics or computing, or inside this project. [NAMING-COLLISIONS.md](../../../NAMING-COLLISIONS.md) says which meaning is in play.

- 39 theorems (35 mathematical plus 3 record-extensionality helpers and one list lemma), core Lean 4.37.0-rc1, no mathlib, no `sorry`.
- Compile status: compiles clean against the shared library build (`LEAN_PATH=/tmp/lib/build`, imports `FusedAlgebra.Base` only).
  `#print axioms` on fold_bound, trans_compose, fold_wall_iff_int, cayley_tan: propext, Classical.choice, Quot.sound only
  (gauss_root_steps: propext only). The #print lines were removed afterwards.
- Intended module path: `FusedAlgebra/Topics/ComplexAnalysis.lean` (namespace `ComplexAnalysis`). Add after Base in build.sh.
- Complex numbers are records `Cx R` over any `Lean.Grind.CommRing` (Int = Gaussian integers); cube roots of unity are
  Eisenstein records `Eis R` with ω² = −1 − ω; Möbius maps are 2×2 records `Mob R`.

Sections: spin and stretch about the wall (ma_lift, spin_about_wall, step_is_turned_w, wall_distance_mul, circle_closed,
circle_equation); affine maps fixing −1 (affine_fixing_wall); circle roots of 0 (mroot_of_unity, gauss_root_four,
eis_root_three, totals); Möbius landscape (fixed_points_are_walls, trans_compose, trans_det, trans_type, two_walls_lift,
cayley_tan, double_wall_adds); products (fold_wall_iff_int, corner_bound, fold_bound); Tate gluing (tate_glue).
Not in Lean (no reals or exp in core): log branch and winding, e^z − 1 petals and the ripple a1, sine and Gamma products,
argument principle. Those are checked in checks.py.
