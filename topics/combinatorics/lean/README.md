# Combinatorics.lean

- 35 theorems, core Lean 4.37.0-rc1, no mathlib, no `sorry`.
- Compiles clean (no errors, no warnings) with
  `LEAN_PATH=/tmp/lib/build lean lean/Combinatorics.lean`.
- `#print axioms` on graphs_from_connected, flow_order4, wall_entry_absorbs, compositions_cube,
  quarter_rung_negative, painted_subsets: only propext, Classical.choice, Quot.sound (some none).
- Imports FusedAlgebra.Base and FusedAlgebra.Cube. Intended module path: FusedAlgebra/Topics/Combinatorics.lean.

Main theorems: groups_length, painted_subsets, nonempty_subsets, alternating_corners, parity_corners,
wfold_lift, wall_entry_absorbs, ffold_fact (n! − 1 = 1 ∘₁ 2 ∘₁ ⋯), kfold_lift, sign_balance,
cube_is_stirling, derangements_incl_excl, bell_is_E_twice, rung_three, graphs_from_connected,
arrangements, flow_order2/3/4 (group law of Eᵗ to order 4), flow_at_one/two/three/minus_one,
half_rung_x4, quarter_rung_negative, S2_rows_bell, dial_blocks_minus_one, compositions_cube.
