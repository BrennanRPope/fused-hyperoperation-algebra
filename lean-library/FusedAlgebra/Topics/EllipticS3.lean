/-
  Fused Hyperoperation Algebra: why the primes' turning angles live on the 3-sphere.
  Core Lean 4 only (no mathlib), no `sorry`. NOT YET COMPILED: Lean is not installed in the
  workspace where this was written. Style follows the shared library (`grind`, `decide`).
  Every identity was also checked with sympy (../s3_checks.py).
  To compile: copy to lean-library/FusedAlgebra/Topics/EllipticS3.lean, add it to build.sh and
  FusedAlgebra.lean, run ./build.sh.

  LIBRARY NOTE: copied from elliptic/lean/EllipticS3.lean (Nov 25) as Topics/EllipticS3 and
  compiled here with core Lean 4.37.0-rc1. circ_is_1x1 is now a one-line reading of the
  base's lift (circ_lift).
-/
import FusedAlgebra.Base

namespace EllipticS3
open CircleOp

/-! ### 1. Every law on the map is the circle operation's lift in 2×2 matrices
    The number 1 + a becomes the matrix I + a·G with G = [[k, 1], [−q, 0]] (trace k, det q).
    Entries of (I + aG)(I + bG) against (1 − q ab)·(I + F G), F = (a + b + k ab)/(1 − q ab),
    written with N = a + b + k ab and D = 1 − q ab so nothing is divided. -/

/-- top-left entry: (1 + ka)(1 + kb) − q a b = D + k N -/
theorem lift_11 (k q a b : Int) :
    (1 + k * a) * (1 + k * b) + a * (-(q * b)) = (1 - q * a * b) + k * (a + b + k * a * b) := by grind

/-- top-right entry: (1 + ka) b + a = N -/
theorem lift_12 (k q a b : Int) : (1 + k * a) * b + a * 1 = a + b + k * a * b := by grind

/-- bottom-left entry: −q a (1 + kb) − q b = −q N -/
theorem lift_21 (k q a b : Int) :
    (-(q * a)) * (1 + k * b) + 1 * (-(q * b)) = -q * (a + b + k * a * b) := by grind

/-- bottom-right entry: −q a b + 1 = D -/
theorem lift_22 (k q a b : Int) : (-(q * a)) * b + 1 * 1 = 1 - q * a * b := by grind

/-- Cayley–Hamilton for G: G² = kG − qI, entry by entry -/
theorem cayley_hamilton (k q : Int) :
    k * k + 1 * (-q) = k * k - q ∧ k * 1 + 1 * 0 = k * 1 - q * 0 ∧
    (-q) * k + 0 * (-q) = k * (-q) - q * 0 ∧ (-q) * 1 + 0 * 0 = k * 0 - q * 1 := by
  refine ⟨?_, ?_, ?_, ?_⟩ <;> grind

/-- the circle operation (called mult-add before 1.0) is the 1×1 case G = [1]: (1 + a)(1 + b) = 1 + (a ◯ b) -/
theorem circ_is_1x1 (a b : Int) : (1 + a) * (1 + b) = 1 + circ a b := (circ_lift a b).symm

/-! ### 2. The turning country has no real wall: no rational eigenvalue when k² < 4q -/

theorem sq_nonneg_rat (x : Rat) : 0 ≤ x * x := by
  rcases (Rat.le_total : (0 : Rat) ≤ x ∨ x ≤ 0) with h | h
  · exact Rat.mul_nonneg h h
  · have h' : 0 ≤ -x := by grind
    have := Rat.mul_nonneg h' h'
    grind

theorem turning_no_wall (k q l : Rat) (h : k * k < 4 * q) : l * l - k * l + q ≠ 0 := by
  have hs := sq_nonneg_rat (2 * l - k)
  intro h0
  have : 4 * (l * l - k * l + q) = (2 * l - k) * (2 * l - k) + (4 * q - k * k) := by grind
  grind

/-! ### 3. Bad primes: Frobenius on the lift is the p-th power -/

/-- turning (non-split) torus: u^(p+1) = 1, so v = u^p satisfies v·u = 1 and the circle operation of
    v − 1 and u − 1 is 0: Frobenius sends the point to its circle inverse (a_p = −1) -/
theorem nonsplit_frob_is_inverse (u v : Int) (h : v * u = 1) : circ (v - 1) (u - 1) = 0 := by
  unfold circ; grind

/-- split torus: u^p = u, Frobenius fixes the lift (a_p = +1) -/
theorem split_frob_fixes (u v : Int) (h : v = u) : circ (v - 1) 0 = u - 1 := by
  unfold circ; grind

/-! ### 4. Points are fixed points: the wall polynomial at 1
    Lefschetz on the torus's cube of holes (layers 1, 2, 1): 1 − (α + β) + αβ = (1 − α)(1 − β) -/
theorem lefschetz_count (α β : Int) : 1 - (α + β) + α * β = (1 - α) * (1 - β) := by grind

/-! ### 5. The 3-sphere: unit quaternions as 2×2 matrices
    U = [[w + ix, y + iz], [−y + iz, w − ix]]; complex numbers written as (re, im) pairs. -/

/-- det U: real part w² + x² + y² + z², imaginary part 0 -/
theorem quat_det_re (w x y z : Int) :
    (w * w - x * (-x)) - ((y * (-y)) - (z * z)) = w * w + x * x + y * y + z * z := by grind
theorem quat_det_im (w x y z : Int) :
    (w * (-x) + x * w) - (y * z + z * (-y)) = 0 := by grind

/-- trace U = 2w; its wall discriminant is −4 times the squared radius of its latitude 2-sphere:
    how deep a prime sits in the turning country is how big the sphere of rotations with its angle is -/
theorem quat_discriminant (w x y z : Int) :
    (2 * w) * (2 * w) - 4 * (w * w + x * x + y * y + z * z) = -4 * (x * x + y * y + z * z) := by grind

/-! ### 6. Moments through the cube
    Corners of the m-cube read as ±1 walks. `walks m h` counts walks of m steps from height h that
    end at 0; `dyck m h` those that never dip below 0. The circle keeps every balanced corner
    (the middle Pascal layer), the 3-sphere keeps only the Dyck ones (Catalan numbers). -/
def walks : Nat → Int → Nat
  | 0, h => if h = 0 then 1 else 0
  | m + 1, h => walks m (h + 1) + walks m (h - 1)

def dyck : Nat → Nat → Nat
  | 0, h => if h = 0 then 1 else 0
  | m + 1, h => dyck m (h + 1) + (if h = 0 then 0 else dyck m (h - 1))

theorem circle_moments : walks 2 0 = 2 ∧ walks 4 0 = 6 ∧ walks 6 0 = 20 ∧ walks 8 0 = 70 := by decide

theorem sphere_moments : dyck 2 0 = 1 ∧ dyck 4 0 = 2 ∧ dyck 6 0 = 5 ∧ dyck 8 0 = 14 := by decide

/-- middle layer minus its neighbour: 20 − 15 = 5 and 70 − 56 = 14 -/
theorem middle_minus_neighbour : (20 : Nat) - 15 = dyck 6 0 ∧ (70 : Nat) - 56 = dyck 8 0 := by decide

end EllipticS3
