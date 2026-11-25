/-
  Fused Hyperoperation Algebra: the corner dial and the rungs below addition.
  Core Lean 4 only (no mathlib).

  STATUS: NOT COMPILED. Unpacking the Lean toolchain was refused in this workspace
  (Nov 25), so these proofs are unchecked. Every statement below is also verified as an
  exact algebraic identity in ../symbolic_proofs.py (sympy), which passes 42 of 42.
  To check: `lean FractionalRungs.lean` with core Lean 4.37 (grind's ring solver handles Int).
-/
namespace FractionalRungs

/-- The corner dial: a + b + k·a·b.  k = 0 addition, k = 1 the circle operation (called mult-add before 1.0), k = -1 "or". -/
def dial (k a b : Int) : Int := a + b + k * a * b

theorem dial_comm (k a b : Int) : dial k a b = dial k b a := by unfold dial; grind
theorem dial_assoc (k a b c : Int) : dial k (dial k a b) c = dial k a (dial k b c) := by
  unfold dial; grind
theorem dial_zero (k a : Int) : dial k a 0 = a := by unfold dial; grind
/-- lift: the dial is multiplication seen through x ↦ 1 + k x -/
theorem dial_lift (k a b : Int) : 1 + k * dial k a b = (1 + k * a) * (1 + k * b) := by
  unfold dial; grind
/-- the dial is the circle operation at another scale -/
theorem dial_rescale (k a b : Int) : k * dial k a b = dial 1 (k * a) (k * b) := by
  unfold dial; grind
/-- mirror: flipping signs flips the corner -/
theorem dial_mirror (k a b : Int) : -(dial k (-a) (-b)) = dial (-k) a b := by
  unfold dial; grind
/-- k = -1 is inclusion–exclusion: 1 - (1-a)(1-b) -/
theorem dial_or (a b : Int) : dial (-1) a b = 1 - (1 - a) * (1 - b) := by unfold dial; grind
/-- three at once: layers of the cube weighted 1, k, k² (one per dimension) -/
theorem dial_three (k a b c : Int) :
    dial k (dial k a b) c = a + b + c + k * (a * b + b * c + c * a) + k * k * (a * b * c) := by
  unfold dial; grind
/-- KdV: Hirota's weight 1 - A = 4 k1 k2 / (k1+k2)², stated without division -/
theorem kdv_weight (k1 k2 : Int) :
    (k1 + k2) * (k1 + k2) - (k1 - k2) * (k1 - k2) = 4 * k1 * k2 := by grind

/-! ### Fusing (conjugating by x ↦ 1 + x) leaves the rungs below addition unchanged -/
def fuse (op : Nat → Nat → Nat) (a b : Nat) : Nat := op (a + 1) (b + 1) - 1

theorem fuse_succ (a b : Nat) : fuse (fun _ y => y + 1) a b = b + 1 := by unfold fuse; omega
theorem fuse_max (a b : Nat) : fuse max a b = max a b := by unfold fuse; omega
theorem fuse_min (a b : Nat) : fuse min a b = min a b := by unfold fuse; omega
/-- addition is the first rung the shift touches -/
theorem fuse_add (a b : Nat) : fuse (· + ·) a b = a + b + 1 := by unfold fuse; omega
theorem fuse_mul (a b : Nat) : fuse (· * ·) a b = a + b + a * b := by
  unfold fuse
  have : (a + 1) * (b + 1) = a + b + a * b + 1 := by grind
  omega

end FractionalRungs
