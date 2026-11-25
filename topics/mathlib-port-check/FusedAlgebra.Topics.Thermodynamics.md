# FusedAlgebra.Topics.Thermodynamics

Theorems (68): GENERAL 65, MATHLIB 2, DROPPED-OK 1. No WEAKER / DIFFERENT / MISSING / VACUOUS. No low-confidence verdicts.

## The removed-hypothesis claims, checked one by one (all hold)
- `carnot_series`: `h ≠ 0` was dropped. At `h = 0`, Lean's `x / 0 = 0` makes both sides equal to 1, and `m ≠ 0` is kept, so the stronger statement is true.
- `reverse_is_inverse`: `1 - eta ≠ 0` was dropped. `multAddInv (-eta) = eta / (1 + -eta)` is literally the same expression, so it holds with no hypothesis.
- `shields`: `1 - r ≠ 0` was dropped. At `r = 1`, `h2` forces `n ≠ 0`, and then both sides are 0.
- `wick_num`: `i * i = -1` was dropped. Both sides equal `ia + ib + k i² ab` for every `i`.
- `wick_den` correctly keeps its `i * i = -1`.

## Other checks
- Added class assumptions are automatic over ℚ: `pauling_ice` needs `CharZero` (and is genuinely false in characteristic 2), `lens_involution` needs `NeZero 2`, `tsallis_composes` needs `IsDomain`.
- `outer` is rebuilt with `flatMap`. It gives the same list in the same order, and `outer_nil`/`outer_cons` hold by `rfl`.
- `fermion_cube`: the original's `groups` has exactly the recursion of `List.sublists'` (`sublists'_cons`), so the lists are identical. `prodL`/`sumL` are `List.prod`/`List.sum`.
- `foldMA` → `multAddList`: same value.
- The multiplicativity hypothesis on `f : ℚ → ℚ` became `MulHomClass`. The original follows by taking the bundled map `⟨f, hf⟩`.
- The `boson_mode` hypothesis `1 - z ≠ 0` is equivalent to the original `1 + -z ≠ 0`.
- `map_mul_map` (helper) is dropped; its content is a `simp` step inside `map_outer`.

## Route changes
- `unattainable`: his list induction is replaced by `one_add_multAddList` plus mathlib's `List.prod_pos`.
- `fermion_lift` and `fermion_cube` now rest on Basic's monoid-hom proofs (`map_list_prod`), not on Base's inductions.
- `two_spin_no_real_zero`, `chain3_*`, `second_law_composes`: hand sign arguments replaced by `nlinarith`/`linarith`.
- Throughout, `grind` became `ring`/`field_simp`. That is a change of tactic, not of argument.

New port theorems with no original: 2 (`outer_nil`, `outer_cons`).
