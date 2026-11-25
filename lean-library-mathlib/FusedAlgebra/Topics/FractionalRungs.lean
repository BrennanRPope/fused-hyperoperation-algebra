module

public import FusedAlgebra.Rungs

/-!
# The corner dial and the rungs below addition

The corner dial `a + b + k * a * b` is the weighted mult-add `MultAdd.wAdd` of
`FusedAlgebra.Basic`; its commutative, associative, identity, lift, rescale and mirror laws are
proved there (`wAdd_comm`, `wAdd_assoc`, `wAdd_zero`, `one_add_mul_wAdd`, `mul_wAdd`,
`wAdd_neg_neg`).

## Main results

* `dial_three`: three at once; the layers of the cube are weighted `1`, `k`, `k ^ 2`.
* `kdv_weight`: Hirota's two-soliton weight for KdV, stated without division.
* `fuse`: conjugating a binary operation on `ℕ` by `x ↦ x + 1`. Fusing leaves the rungs below
  addition (successor, `max`, `min`) unchanged; addition is the first rung the shift touches
  (`fuse_succ`, `fuse_max`, `fuse_min`, `fuse_add`, `fuse_mul`). Fusing mathlib's
  `hyperoperation n` is the ladder `MultAdd.rung n` (`fuse_hyperoperation`).

The inclusion-exclusion law at `k = -1` (`dial_or` in the core-Lean original) is
`MultAdd.wAdd_neg_one`. The identities were stated over `ℤ` in the original and hold in any
commutative (semi)ring.
-/

@[expose] public section

namespace FractionalRungs

open MultAdd

/-- Three at once: the layers of the cube are weighted `1`, `k`, `k ^ 2` (one factor of `k` per
dimension). -/
theorem dial_three {R : Type*} [CommSemiring R] (k a b c : R) :
    wAdd k (wAdd k a b) c
      = a + b + c + k * (a * b + b * c + c * a) + k * k * (a * b * c) := by
  unfold wAdd; ring

/-- KdV: Hirota's weight `1 - A = 4 * k₁ * k₂ / (k₁ + k₂) ^ 2`, stated without division. -/
theorem kdv_weight {R : Type*} [CommRing R] (k₁ k₂ : R) :
    (k₁ + k₂) * (k₁ + k₂) - (k₁ - k₂) * (k₁ - k₂) = 4 * k₁ * k₂ := by
  ring

/-! ### Fusing (conjugating by `x ↦ x + 1`) leaves the rungs below addition unchanged -/

/-- The fusion of a binary operation on `ℕ`: the same operation seen through `x ↦ x + 1`. -/
def fuse (op : ℕ → ℕ → ℕ) (a b : ℕ) : ℕ := op (a + 1) (b + 1) - 1

/-- Fusing the hyperoperation `n` gives rung `n` of the fused ladder. -/
theorem fuse_hyperoperation (n : ℕ) : fuse (hyperoperation n) = rung n := rfl

theorem fuse_succ (a b : ℕ) : fuse (fun _ y ↦ y + 1) a b = b + 1 := rfl

theorem fuse_max (a b : ℕ) : fuse max a b = max a b := by unfold fuse; lia

theorem fuse_min (a b : ℕ) : fuse min a b = min a b := by unfold fuse; lia

/-- Addition is the first rung the shift touches: fused addition is `sadd`. -/
theorem fuse_add (a b : ℕ) : fuse (· + ·) a b = a + b + 1 := by unfold fuse; lia

/-- Fused multiplication is mult-add. -/
theorem fuse_mul (a b : ℕ) : fuse (· * ·) a b = a ✳ b := by
  rw [fuse, ← multAdd_add_one, Nat.add_sub_cancel]

end FractionalRungs
