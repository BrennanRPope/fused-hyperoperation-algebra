/-
  Fused Hyperoperation Algebra: the corner wAdd and the rungs below addition.
  Core Lean 4 only (no mathlib).

  STATUS: the original file was not compiled where it was written (Nov 25); in the
  shared library it compiles with core Lean 4.37.0-rc1, no `sorry` (three `fuse`
  proofs needed a `dsimp only` before `omega`). Every statement below is also verified as an
  exact algebraic identity in ../symbolic_proofs.py (sympy), which passes 42 of 42.
  To check: `lean FractionalRungs.lean` with core Lean 4.37 (grind's ring solver handles Int).

  LIBRARY NOTE: copied from the finished fractional-rungs file (Nov 25).
-/
import FusedAlgebra.Base
import FusedAlgebra.Rungs

namespace FractionalRungs

open MultAdd

/-! The corner dial a + b + k·a·b is the weighted mult-add `wAdd` of FusedAlgebra.Base;
  its commutative, associative, identity, lift, rescale and mirror laws are
  wAdd_comm, wAdd_assoc, wAdd_zero, wAdd_lift, wAdd_scale and wAdd_flip there. -/

/-- k = -1 is inclusion–exclusion: 1 - (1-a)(1-b) -/
theorem dial_or (a b : Int) : wAdd (-1) a b = 1 - (1 - a) * (1 - b) := by unfold wAdd; grind
/-- three at once: layers of the cube weighted 1, k, k² (one per dimension) -/
theorem dial_three (k a b c : Int) :
    wAdd k (wAdd k a b) c = a + b + c + k * (a * b + b * c + c * a) + k * k * (a * b * c) := by
  unfold wAdd; grind
/-- KdV: Hirota's weight 1 - A = 4 k1 k2 / (k1+k2)², stated without division -/
theorem kdv_weight (k1 k2 : Int) :
    (k1 + k2) * (k1 + k2) - (k1 - k2) * (k1 - k2) = 4 * k1 * k2 := by grind

/-! ### Fusing (conjugating by x ↦ 1 + x) leaves the rungs below addition unchanged -/
def fuse (op : Nat → Nat → Nat) (a b : Nat) : Nat := op (a + 1) (b + 1) - 1

theorem fuse_succ (a b : Nat) : fuse (fun _ y => y + 1) a b = b + 1 := by unfold fuse; dsimp only; omega
theorem fuse_max (a b : Nat) : fuse max a b = max a b := by unfold fuse; omega
theorem fuse_min (a b : Nat) : fuse min a b = min a b := by unfold fuse; omega
/-- addition is the first rung the shift touches -/
theorem fuse_add (a b : Nat) : fuse (· + ·) a b = a + b + 1 := by unfold fuse; dsimp only; omega
theorem fuse_mul (a b : Nat) : fuse (· * ·) a b = a + b + a * b := by
  unfold fuse; dsimp only
  have : (a + 1) * (b + 1) = a + b + a * b + 1 := by grind
  omega

end FractionalRungs
