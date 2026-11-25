module

public import FusedAlgebra.Basic

/-!
# The layers of the mult-add cube

`1 + mpow a n = (1 + a) ^ n`, and the `n`-cube behind it has `n.choose k` corners with `k`
switches on (Pascal's triangle).

This is the mathlib form of `FusedAlgebra.Cube`. The core-Lean library built the layer counts
`layer n k` and the factorial `fact` by hand; in mathlib they are `Nat.choose` and
`Nat.factorial`, so those definitions and their lemmas are dropped (see `PORT-NOTES.md` for the
name-by-name mapping). What remains here are the statements that tie the layers to mult-add.
-/

@[expose] public section

namespace MultAdd.Cube

open Finset

variable {R : Type*}

/-- The cube by layers: `1 + mpow a n` is the sum over grades `k` of `n.choose k` corners, each of
measure `a ^ k`. -/
theorem one_add_mpow_eq_sum_choose [CommSemiring R] (a : R) (n : ℕ) :
    1 + mpow a n = ∑ k ∈ range (n + 1), (n.choose k : R) * a ^ k := by
  rw [one_add_mpow, add_comm, add_pow]
  exact sum_congr rfl fun k _ ↦ by rw [one_pow, mul_one, mul_comm]

/-- The cube by layers, solved for the mult-add power: the grade-`0` corner is the unit that the
lift adds. -/
theorem mpow_eq_sum_choose [CommRing R] (a : R) (n : ℕ) :
    mpow a n = ∑ k ∈ range n, (n.choose (k + 1) : R) * a ^ (k + 1) := by
  have h := one_add_mpow_eq_sum_choose a n
  rw [sum_range_succ', Nat.choose_zero_right, Nat.cast_one, pow_zero, mul_one, add_comm] at h
  exact add_right_cancel h

/-- There are `n.choose k` corners of grade `k` in the `n`-cube. -/
theorem card_corners (n k : ℕ) : #{s : Finset (Fin n) | #s = k} = n.choose k := by
  have h : ({s | #s = k} : Finset (Finset (Fin n))) = powersetCard k univ := by
    ext s; simp [mem_powersetCard]
  rw [h, card_powersetCard, card_univ, Fintype.card_fin]

/-- The `n`-cube has `2 ^ n` corners: `mpow 1 n + 1 = 2 ^ n`. -/
theorem mpow_one_add_one (n : ℕ) : mpow 1 n + 1 = 2 ^ n := by
  rw [add_comm, one_add_mpow]

end MultAdd.Cube
