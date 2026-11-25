# DiamondLattice.lean

- 48 theorems (2 of them small division helpers: quot_helper, slope_helper), core Lean 4.37.0-rc1, no mathlib, no `sorry`.
- Compiles clean (no errors, no warnings) with `LEAN_PATH=/tmp/lib/build lean lean/DiamondLattice.lean`.
- `#print axioms` on lift_mul, slope_at_fixed, dmap_fixed_param, q_two_inversion, mirror_only_three, fold_identity:
  only propext, Classical.choice, Quot.sound.
- Imports FusedAlgebra.Base and FusedAlgebra.Topics.Statistics (for `StatFused.lens`).
  Intended module path: FusedAlgebra/Topics/DiamondLattice.lean (builds after Topics/Statistics).

Main theorems: dmap_formula, dmap_even, quotient_rule, wall_one, wall_colour, colour_wall_doubles, lift_mul, lift_weight,
lift_dmap (one step = square in t, then square in the lift), disc_square (p^2 - 4q = Q^2), wall_gap, bend_recip, fixed_factor,
param_cubic, param_fixed, dmap_fixed_param (t = 1/y fixed when Q = (y-1)^2(y+1)), crs_cubic, dden_at_fixed,
slope_at_fixed (slope 4/(1+t)^2 at every nontrivial fixed point), slope_in_y, three/sixteen/fortyfive fixed and slope,
golden_slope, ising_slope, fold_identity (32/27 - Q(y) = (y+1/3)^2(5/3-y)), fold_ferro_fixed/slope, fold_neutral,
attract_identity, q_zero_map/double/parabolic, q_one_poly, q_two_inversion, q_three_mirror, mirror_only_three, parallel_weights.
