module

public import FusedAlgebra.Basic

/-!
# Gravity Probe B through "twist = curvature × area"

A gyroscope in a circular orbit turns, per orbit, by two twists: the space twist, a fraction
`p_s = 1 - √(1 - 2 * M / r)` of a turn, and the velocity twist, a fraction `p_v = 1 - 1 / γ` with
`γ ^ 2 = (r - 2 * M) / (r - 3 * M)`. They combine by multiplying the factors,
`1 - p = (1 - p_s) * (1 - p_v)`, so `-p = -p_s ✳ -p_v`: **the precession fractions mult-add**.
The total is the known exact geodetic precession `1 - √(1 - 3 * M / r)` per orbit.

Here `x = M / r` (in units `G = c = 1`), and the file works with the squares of the factors, so
no square roots are needed.

## Main results

* `gamma_sq`, `factors_multiply`: `(1 - 2 * x) * (1 / γ ^ 2) = 1 - 3 * x`.
* `fractions_fuse`, `fuse_factors`: the precession fractions mult-add.
* `weak_split`, `thirds`: in the weak field the total is `3 / 2 * x`, one third velocity and two
  thirds space.
* `geodetic_within_error`, `parts_alone_ruled_out`: the Gravity Probe B numbers.

The physics (the geodetic effect, its split into a Thomas part and a space-curvature part, the
measured values) is standard; the algebra here is elementary. The core-Lean original worked over
`ℚ`; the identities are stated in any commutative ring or field (of characteristic zero where a
numeral is divided by), and the measured numbers stay in `ℚ`.
-/

@[expose] public section

namespace GravityProbeB

open MultAdd

variable {R : Type*} [CommRing R] {K : Type*} [Field K]

/-! ### The two factors and their product -/

/-- Local orbital speed squared, relative to static observers: `v ^ 2 = x / (1 - 2 * x)`. -/
def v2 (x : K) : K := x / (1 - 2 * x)

/-- `γ ^ 2 = 1 / (1 - v ^ 2) = (1 - 2 * x) / (1 - 3 * x)`. -/
theorem gamma_sq (x : K) (h2 : 1 - 2 * x ≠ 0) : 1 / (1 - v2 x) = (1 - 2 * x) / (1 - 3 * x) := by
  have e : 1 - x / (1 - 2 * x) = (1 - 3 * x) / (1 - 2 * x) := by
    rw [eq_div_iff h2, sub_mul, div_mul_cancel₀ _ h2]; ring
  rw [v2, e, one_div_div]

/-- **The factors multiply to the known total** (squared):
`(1 - 2 * x) * (1 / γ ^ 2) = 1 - 3 * x`. -/
theorem factors_multiply (x : K) (h2 : 1 - 2 * x ≠ 0) : (1 - 2 * x) * (1 - v2 x) = 1 - 3 * x := by
  unfold v2; field_simp; ring

/-- **The precession fractions mult-add**: if `1 - p = (1 - a) * (1 - b)` then
`-p = -a ✳ -b`. -/
theorem fractions_fuse (p a b : R) (h : 1 - p = (1 - a) * (1 - b)) : -p = -a ✳ -b := by
  unfold multAdd; linear_combination h

/-- And back: the mult-add form gives the product of factors. -/
theorem fuse_factors (a b : R) : 1 + (-a ✳ -b) = (1 - a) * (1 - b) := by
  unfold multAdd; ring

/-- The same law is the dial's weight `-1`: `p = a + b - a * b` ("either twist"). See also
`MultAdd.wAdd_neg_one_eq_neg_multAdd`. -/
theorem fuse_weight_minus_one (a b : R) : -(-a ✳ -b) = a + b + (-1) * a * b := by
  unfold multAdd; ring

/-! ### Weak field: one part velocity, two parts space -/

/-- First-order parts: space `1 * x`, velocity `1 / 2 * x`, total `3 / 2 * x`. -/
theorem weak_split [CharZero K] (x : K) : x + x / 2 = 3 / 2 * x := by ring

/-- The fused total differs from the plain sum only by the corner `a * b` (second order). -/
theorem corner_is_second_order (a b : R) : (a + b) - (-(-a ✳ -b)) = a * b := by
  unfold multAdd; ring

/-- In PPN form the geodetic coefficient is `γ_PPN + 1 / 2`; general relativity's `γ_PPN = 1`
gives `3 / 2`. -/
theorem ppn_geodetic [CharZero K] : (1 : K) + 1 / 2 = 3 / 2 := by norm_num

/-- The velocity part alone is one third, the space part alone two thirds. -/
theorem thirds [CharZero K] : (1 / 2 : K) / (3 / 2) = 1 / 3 ∧ (1 : K) / (3 / 2) = 2 / 3 := by
  norm_num

/-! ### The Gravity Probe B numbers (milliarcseconds per year) -/

/-- Twist bookkeeping at `r = 7027.4` km gives `6603.9`; GP-B measured `6601.8 ± 18.3`. -/
theorem geodetic_within_error :
    (66039 / 10 : ℚ) - 66018 / 10 < 183 / 10 ∧ (66018 / 10 : ℚ) - 66039 / 10 < 183 / 10 := by
  norm_num

/-- Space part `4402.6` plus velocity part `2201.3` is `6603.9`. -/
theorem geodetic_parts : (44026 / 10 : ℚ) + 22013 / 10 = 66039 / 10 := by norm_num

/-- Either part alone misses the measurement by more than `100` standard errors. -/
theorem parts_alone_ruled_out :
    100 * (183 / 10 : ℚ) < 66018 / 10 - 22013 / 10 ∧
      100 * (183 / 10 : ℚ) < 66018 / 10 - 44026 / 10 := by
  norm_num

/-- Frame dragging: about `40.8` against the measured `37.2 ± 7.2`. -/
theorem frame_dragging_within_error : (408 / 10 : ℚ) - 372 / 10 < 72 / 10 := by norm_num

/-! ### Frame dragging over a polar orbit -/

/-- `3 * sin ^ 2 - 1` averaged at the quarter points (`sin ^ 2 = 0, 1, 0, 1`) is `1 / 2`. -/
theorem polar_average_four [CharZero K] :
    ((3 * 0 - 1) + (3 * 1 - 1) + (3 * 0 - 1) + (3 * 1 - 1) : K) / 4 = 1 / 2 := by
  norm_num

/-- At the eighth points (`sin ^ 2 = 0, 1 / 2, 1, 1 / 2, …`) the average is `1 / 2` too. -/
theorem polar_average_eight [CharZero K] :
    ((3 * 0 - 1) + (3 * (1 / 2) - 1) + (3 * 1 - 1) + (3 * (1 / 2) - 1)
      + (3 * 0 - 1) + (3 * (1 / 2) - 1) + (3 * 1 - 1) + (3 * (1 / 2) - 1) : K) / 8 = 1 / 2 := by
  norm_num

/-- Any orbit point and its quarter-turn partner have `sin ^ 2 + cos ^ 2 = 1`, so a pair averages
to `1 / 2` whatever `s` is: the average is `1 / 2` for every set of points closed under quarter
turns. -/
theorem polar_average_pair [CharZero K] (s : K) :
    ((3 * s - 1) + (3 * (1 - s) - 1)) / 2 = 1 / 2 := by
  ring

end GravityProbeB
