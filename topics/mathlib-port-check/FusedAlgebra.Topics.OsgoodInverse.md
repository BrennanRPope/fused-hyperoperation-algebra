# FusedAlgebra.Topics.OsgoodInverse

Theorems (53): SAME 20, GENERAL 16, MATHLIB 14, DROPPED-OK 2, WEAKER 1. No DIFFERENT / MISSING / VACUOUS. No low-confidence verdicts.

## Findings
- **WEAKER: `const_block` → `blk_const`** (the narrowing that PORT-NOTES admits; confirmed).
  Original: `∀ s, fsum (fun _ => c) s (2 ^ j) = (2 : Rat) ^ j * c`, so any start `s`.
  Port: `blk (fun _ ↦ c) j = 2 ^ j * c`, so only start `2 ^ j`, though over any field.
  The original used it only at `s = 2 ^ j` (in `blk_le` and `le_blk`), and the general form is a one-line mathlib fact (`Finset.sum_const`, `Nat.card_Ico`). So nothing downstream is lost, but the library no longer states the general lemma.

## Checked and fine
- Definitions were rebuilt as Finset sums and the index sets match exactly: `dsum` = `Ico 1 (2^n)` (same as `fsum a 1 (2^n - 1)`, empty at n = 0), `blk` = `Ico (2^j) (2^(j+1))`, `csum` = `range n`, `fsum a 1 N` = `Ico 1 (N+1)`. `Antitone1` ⇔ `AntitoneOn a (Set.Ici 1)`. `steps h n` = `n • h`. `nr` = `Nat.cast`.
- `mirror_lift` is covered by `MultAdd.inv_one_add` (same equation, sides swapped). `steps_quarter` (`4 • (1/4 : ℚ) = 1`) is now only a `norm_num` step inside `riccati_quarter_blows_up`, with the same numbers.
- `pow2_succ` → `Nat.two_pow_succ` (medium): the port never uses this name, so the name itself is unchecked.

## Route changes
- The squeezes `dsum_le_csum`, `csum_le_dsum` and `condensation` keep his own block-by-block proofs. Mathlib's condensation lemmas are only cited in doc comments.
- `fsum_le_dsum` (now `sum_Ico_le_dsum`): his split with `fsum_append` is replaced by mathlib's `sum_le_sum_of_subset_of_nonneg`.
- `harm_anti` and `isq_anti`: his `PiMultAdd.div_le_div_of` is replaced by mathlib's `one_div_le_one_div_of_le`.
- `time_left_pos`, `harm_nonneg`: hand proofs replaced by `positivity`. `csum_harm`: induction replaced by `simp`.
- `block_start`: the NSContinuum induction is replaced by `Cube.mpow_one_add_one`, which goes through the monoid-hom proof of `one_add_mpow`.
- `rung1_never_wall`: his own induction is replaced by `rung1_lift` plus `mul_pos`/`pow_pos`. `blk_const` is proved by `sum_const`/`card_Ico` rather than by doubling induction.

New port theorems with no original: 1 (`csum_succ`).
