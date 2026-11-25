module

public import FusedAlgebra.Cube

/-!
# When the Osgood test can be turned around

The Osgood integral is a clock: in the coordinate `G y = ∫ dy / f y` the solution of `y' = f y`
moves at unit speed, so the integral is the time to reach the end. This file holds the exact
algebraic and discrete forms of that picture; no integral or limit is taken.

1. The Riccati climb `y' = (1 + y) ^ 2`: its solution `ric t = t / (1 - t)` obeys the exact
   difference law `ric_step`, `ric_step_multAdd`, and passes every height before `t = 1`
   (`ric_near_one`).
2. The discrete climb with that law hits the wall: its clock `(1 + y)⁻¹` falls by exactly `h` per
   step (`riccati_clock`), so it cannot take `m` steps when the clock starts at `m • h`
   (`riccati_hits_wall`). The rung-1 climb `y ↦ y ✳ h` never does (`rung1_never_wall`). Through
   the mult-add inverse the blow-up is a walk at constant speed into the wall `-1`
   (`mirror_step`).
3. Non-uniqueness: the trapezoid scheme of `x' = 2 * √x` has a whole fan of solutions out of `0`
   (`fan_step`, `fan_distinct`); the same scheme for the Lipschitz `x' = x` has only the zero
   solution (`lipschitz_unique`).
4. The test moves down one rung exactly: Cauchy condensation (`dsum_le_csum`, `csum_le_dsum`,
   `condensation`). The blocks `[2 ^ j, 2 ^ (j + 1))` start at `1 + mpow 1 j` (`block_start`).
   Instances: `1 / k` grows at least like `n / 2` (`harmonic_grows`); `1 / k ^ 2` stays below `2`
   (`inverse_squares_bounded`).
5. Direction decides: the rotation field and the radial field have the same speed `r ^ 3`, but
   the first is perpendicular to `x` and the second pushes straight out (`rot_perp`, `rad_out`).

This is the mathlib form of `FusedAlgebra.Topics.OsgoodInverse`. Cauchy condensation is known
mathematics and is in mathlib: the two squeezes are `Finset.le_sum_condensed'` and (in a slightly
different normalisation) `Finset.sum_condensed_le`, and the summability statement is
`summable_condensed_iff_of_nonneg`. They are re-derived here from the dyadic blocks, to keep the
block picture. The hand-made sums are `Finset.sum`s over `Finset.Ico` and `Finset.range`, the
hand-made cast `nr` is `Nat.cast`, and "shrinking from index `1` on" is `AntitoneOn a (Set.Ici 1)`.
-/

@[expose] public section

namespace MultAdd.OsgoodInverse

open MultAdd Finset

/-! ### The Riccati climb `y' = (1 + y) ^ 2` and the wall it hits -/

section Field

variable {K : Type*} [Field K]

/-- The solution of `y' = (1 + y) ^ 2`, `y 0 = 0`. -/
def ric (t : K) : K := t / (1 - t)

/-- Under the lift `ric` is a reciprocal: `1 + y = (1 - t)⁻¹`. The time left, `1 - t`, is the
clock. -/
theorem one_add_ric (t : K) (h : 1 - t ≠ 0) : 1 + ric t = (1 - t)⁻¹ := by
  unfold ric; field_simp; ring

/-- The exact difference law: `y (t + h) - y t = h * (1 + y t) * (1 + y (t + h))`, the discrete
form of `y' = (1 + y) ^ 2` with no error term. -/
theorem ric_step (t h : K) (h1 : 1 - t ≠ 0) (h2 : 1 - (t + h) ≠ 0) :
    ric (t + h) - ric t = h * (1 + ric t) * (1 + ric (t + h)) := by
  rw [one_add_ric _ h1, one_add_ric _ h2, ric, ric]; field_simp; ring

/-- The same law in mult-add: the step is `h` times the lift of `y t ✳ y (t + h)`. -/
theorem ric_step_multAdd (t h : K) (h1 : 1 - t ≠ 0) (h2 : 1 - (t + h) ≠ 0) :
    ric (t + h) - ric t = h * (1 + (ric t ✳ ric (t + h))) := by
  rw [ric_step t h h1 h2, one_add_multAdd, mul_assoc]

/-- Read in the clock coordinate, the climb is uniform: one step of the law lowers `(1 + y)⁻¹` by
exactly `h`. -/
theorem riccati_clock_step (a b h : K) (ha : 1 + a ≠ 0) (hb : 1 + b ≠ 0)
    (hs : b - a = h * (1 + a) * (1 + b)) : (1 + b)⁻¹ = (1 + a)⁻¹ - h := by
  have key : (1 + a)⁻¹ - (1 + b)⁻¹ = h := by
    rw [inv_sub_inv ha hb, div_eq_iff (mul_ne_zero ha hb)]
    linear_combination hs
  linear_combination -key

/-- The discrete Riccati climb runs on a straight clock: as long as the lift stays off the wall,
`(1 + y n)⁻¹ = (1 + y 0)⁻¹ - n • h`. -/
theorem riccati_clock (y : ℕ → K) (h : K)
    (hs : ∀ n, y (n + 1) - y n = h * (1 + y n) * (1 + y (n + 1))) (n : ℕ)
    (hne : ∀ k ≤ n, 1 + y k ≠ 0) : (1 + y n)⁻¹ = (1 + y 0)⁻¹ - n • h := by
  induction n with
  | zero => simp
  | succ n ih =>
    rw [riccati_clock_step _ _ h (hne n n.le_succ) (hne (n + 1) le_rfl) (hs n),
      ih fun k hk ↦ hne k (hk.trans n.le_succ), succ_nsmul]
    ring

/-- Blow-up: if the clock starts at `m • h`, the climb cannot take `m` steps with finite values.
At step `m` the clock reads `0`, and no finite `y` has `(1 + y)⁻¹ = 0`. (Continuous form: the
blow-up time is `∫ ds / (1 + s) ^ 2 = 1`.) -/
theorem riccati_hits_wall (y : ℕ → K) (h : K) (m : ℕ)
    (hs : ∀ n, y (n + 1) - y n = h * (1 + y n) * (1 + y (n + 1)))
    (hstart : (1 + y 0)⁻¹ = m • h) : ¬ ∀ k ≤ m, 1 + y k ≠ 0 := by
  intro hne
  have hc := riccati_clock y h hs m hne
  rw [hstart, sub_self] at hc
  exact hne m le_rfl (inv_eq_zero.1 hc)

/-- The mirror: through the mult-add inverse `v = -u / (1 + u)`, a Riccati step is a step of
exactly `-h`. The climb to infinity is a walk at constant speed into the wall `-1` (which the
mirrored walk reaches exactly when the climb blows up, by `MultAdd.inv_one_add`). -/
theorem mirror_step (a b h : K) (ha : 1 + a ≠ 0) (hb : 1 + b ≠ 0)
    (hs : b - a = h * (1 + a) * (1 + b)) : multAddInv b - multAddInv a = -h := by
  have hc := riccati_clock_step a b h ha hb hs
  rw [inv_one_add b hb, inv_one_add a ha] at hc
  linear_combination hc

/-- Lipschitz means unique: the trapezoid scheme of `x' = x`,
`x (n + 1) - x n = h / 2 * (x n + x (n + 1))`, started at `0` stays at `0` (for `h ≠ 2`). -/
theorem lipschitz_unique (x : ℕ → K) (h : K) (hh : 2 - h ≠ 0) (h0 : x 0 = 0)
    (hs : ∀ n, x (n + 1) - x n = h / 2 * (x n + x (n + 1))) (n : ℕ) : x n = 0 := by
  induction n with
  | zero => exact h0
  | succ n ih =>
    have e := hs n
    rw [ih, sub_zero, zero_add] at e
    by_cases h2 : (2 : K) = 0
    · rwa [h2, div_zero, zero_mul] at e
    · have hp : x (n + 1) * (2 - h) = 0 := by
        linear_combination 2 * e + x (n + 1) * mul_div_cancel₀ h h2
      exact (mul_eq_zero.1 hp).resolve_right hh

end Field

section Semiring

variable {R : Type*} [Semiring R]

/-- The rung-1 climb: if `y (n + 1) = y n ✳ h` then `1 + y n = (1 + y 0) * (1 + h) ^ n`. The
clock `log (1 + y)` advances by `log (1 + h)` per step forever. -/
theorem rung1_lift (y : ℕ → R) (h : R) (hs : ∀ n, y (n + 1) = y n ✳ h) (n : ℕ) :
    1 + y n = (1 + y 0) * (1 + h) ^ n := by
  induction n with
  | zero => simp
  | succ n ih => rw [hs n, one_add_multAdd, ih, pow_succ, mul_assoc]

/-- In an ordered semiring the rung-1 climb stays on the safe side of the wall forever. -/
theorem rung1_never_wall [PartialOrder R] [IsStrictOrderedRing R] (y : ℕ → R) (h : R)
    (hh : 0 < 1 + h) (h0 : 0 < 1 + y 0) (hs : ∀ n, y (n + 1) = y n ✳ h) (n : ℕ) :
    0 < 1 + y n := by
  rw [rung1_lift y h hs n]; exact mul_pos h0 (pow_pos hh n)

end Semiring

/-- With `h = 1 / 4` the Riccati clock started at `1` reads `0` after `4` steps: the climb from
`y 0 = 0` blows up at step `4 = 1 / h`. -/
theorem riccati_quarter_blows_up (y : ℕ → ℚ) (h0 : y 0 = 0)
    (hs : ∀ n, y (n + 1) - y n = (1 / 4) * (1 + y n) * (1 + y (n + 1))) :
    ¬ ∀ k ≤ 4, 1 + y k ≠ 0 :=
  riccati_hits_wall y (1 / 4) 4 hs (by rw [h0]; norm_num)

section LinearOrderedField

variable {K : Type*} [Field K] [LinearOrder K] [IsStrictOrderedRing K]

/-- The continuous climb passes every height before `t = 1`: `y (1 - 1 / (M + 2)) = M + 1`. -/
theorem ric_near_one (M : ℕ) : ric (1 - 1 / ((M : K) + 2)) = M + 1 := by
  have hM : (M : K) + 2 ≠ 0 := by positivity
  unfold ric; field_simp; ring

/-- That time is still before `1`: the time left `1 / (M + 2)` is positive. -/
theorem time_left_pos (M : ℕ) : 0 < 1 / ((M : K) + 2) := by positivity

end LinearOrderedField

/-! ### Non-uniqueness out of `0`: `x' = 2 * √x` versus `x' = x`

The trapezoid scheme of `x' = 2 * √x` on the integer grid is
`x (n + 1) = x n + √(x n) + √(x (n + 1))`. For every start time `c` the curve
`fan c n = (n - c) ^ 2` (zero until `c`, truncated subtraction) solves it with `√x = n - c`: the
solution may sit at `0` as long as it likes and then leave. -/

/-- The fan member that leaves `0` at time `c`. -/
def fan (c n : ℕ) : ℕ := (n - c) ^ 2

/-- The square root of the fan member. -/
def root (c n : ℕ) : ℕ := n - c

theorem root_sq (c n : ℕ) : root c n ^ 2 = fan c n := rfl

/-- Every fan member solves the trapezoid scheme of `x' = 2 * √x`. -/
theorem fan_step (c n : ℕ) : fan c (n + 1) = fan c n + root c n + root c (n + 1) := by
  unfold fan root
  rcases le_or_gt c n with h | h
  · obtain ⟨m, rfl⟩ := Nat.exists_eq_add_of_le h
    have e1 : c + m + 1 - c = m + 1 := by lia
    have e2 : c + m - c = m := by lia
    rw [e1, e2]; ring
  · rw [Nat.sub_eq_zero_of_le h, Nat.sub_eq_zero_of_le h.le]; simp

/-- Every fan member starts at `0`. -/
@[simp]
theorem fan_zero (c : ℕ) : fan c 0 = 0 := by simp [fan]

/-- Different start times give different solutions, so the scheme has infinitely many solutions
from `x 0 = 0`. -/
theorem fan_distinct (c d : ℕ) (h : c < d) : fan c ≠ fan d := by
  intro heq
  have hd : (d - c) ^ 2 = (d - d) ^ 2 := congrFun heq d
  rw [Nat.sub_self, zero_pow two_ne_zero, pow_eq_zero_iff two_ne_zero] at hd
  lia

/-! ### The test moves down one rung exactly (Cauchy condensation)

The substitution `1 + s = exp x` turns `∫ ds / ((1 + s) * g s)` into `∫ dx / g (E x)`: the test
for `f = (1 + s) * g` is the test for `g ∘ E`, one rung down. The discrete form is Cauchy
condensation: sample at the points `1 + mpow 1 j = 2 ^ j`, and weight each sample by the block
length `2 ^ j` (the lift factor). -/

/-- The blocks are the mult-add lattice: the `j`-th block starts at `2 ^ j = mpow 1 j + 1`, the
lift of the `j`-fold mult-add of `1`. -/
theorem block_start (j : ℕ) : 2 ^ j = mpow 1 j + 1 :=
  (Cube.mpow_one_add_one j).symm

section Condensation

variable {K : Type*} [Field K]

/-- The sum of `a` over the `j`-th block `[2 ^ j, 2 ^ (j + 1))`. -/
def blk (a : ℕ → K) (j : ℕ) : K := ∑ k ∈ Ico (2 ^ j) (2 ^ (j + 1)), a k

/-- Dyadic partial sums: `a 1 + ⋯ + a (2 ^ n - 1)`. -/
def dsum (a : ℕ → K) (n : ℕ) : K := ∑ k ∈ Ico 1 (2 ^ n), a k

/-- Condensed partial sums: `∑ j < n, 2 ^ j * a (2 ^ j)`. -/
def csum (a : ℕ → K) (n : ℕ) : K := ∑ j ∈ range n, 2 ^ j * a (2 ^ j)

theorem dsum_succ (a : ℕ → K) (n : ℕ) : dsum a (n + 1) = dsum a n + blk a n :=
  (sum_Ico_consecutive a Nat.one_le_two_pow (Nat.pow_le_pow_right two_pos n.le_succ)).symm

theorem csum_succ (a : ℕ → K) (n : ℕ) : csum a (n + 1) = csum a n + 2 ^ n * a (2 ^ n) :=
  sum_range_succ _ _

/-- A constant over the `j`-th block sums to `2 ^ j * c`. -/
theorem blk_const (c : K) (j : ℕ) : blk (fun _ ↦ c) j = 2 ^ j * c := by
  rw [blk, sum_const, Nat.card_Ico, pow_succ, mul_two, Nat.add_sub_cancel, nsmul_eq_mul]
  norm_num

variable [LinearOrder K] [IsStrictOrderedRing K]

theorem blk_le (a : ℕ → K) (ha : AntitoneOn a (Set.Ici 1)) (j : ℕ) :
    blk a j ≤ 2 ^ j * a (2 ^ j) := by
  rw [← blk_const]
  refine sum_le_sum fun k hk ↦ ?_
  have h1 : 1 ≤ 2 ^ j := Nat.one_le_two_pow
  have h2 : 2 ^ j ≤ k := (mem_Ico.1 hk).1
  exact ha h1 (h1.trans h2) h2

theorem le_blk (a : ℕ → K) (ha : AntitoneOn a (Set.Ici 1)) (j : ℕ) :
    2 ^ j * a (2 ^ (j + 1)) ≤ blk a j := by
  rw [← blk_const]
  refine sum_le_sum fun k hk ↦ ?_
  have h1 : 1 ≤ 2 ^ j := Nat.one_le_two_pow
  have h2 := mem_Ico.1 hk
  exact ha (h1.trans h2.1) (h1.trans (h2.1.trans h2.2.le)) h2.2.le

/-- Upper squeeze: `a 1 + ⋯ + a (2 ^ n - 1) ≤ ∑ j < n, 2 ^ j * a (2 ^ j)`. This is mathlib's
`Finset.le_sum_condensed'`. -/
theorem dsum_le_csum (a : ℕ → K) (ha : AntitoneOn a (Set.Ici 1)) (n : ℕ) :
    dsum a n ≤ csum a n := by
  induction n with
  | zero => simp [dsum, csum]
  | succ n ih => rw [dsum_succ, csum_succ]; exact add_le_add ih (blk_le a ha n)

/-- Lower squeeze: `∑ j ≤ n, 2 ^ j * a (2 ^ j) ≤ a 1 + 2 * (a 1 + ⋯ + a (2 ^ n - 1))`. Compare
mathlib's `Finset.sum_condensed_le`. -/
theorem csum_le_dsum (a : ℕ → K) (ha : AntitoneOn a (Set.Ici 1)) (n : ℕ) :
    csum a (n + 1) ≤ a 1 + 2 * dsum a n := by
  induction n with
  | zero => simp [dsum, csum]
  | succ n ih =>
    rw [dsum_succ, csum_succ, pow_succ]
    linarith [le_blk a ha n]

/-- The test is the same one rung down: for shrinking nonnegative terms, the dyadic partial sums
are bounded iff the condensed sums are. The summability form is mathlib's
`summable_condensed_iff_of_nonneg`. -/
theorem condensation (a : ℕ → K) (ha : AntitoneOn a (Set.Ici 1)) (hpos : ∀ k, 0 ≤ a k) :
    (∃ B, ∀ n, csum a n ≤ B) ↔ ∃ B, ∀ n, dsum a n ≤ B := by
  constructor
  · rintro ⟨B, hB⟩
    exact ⟨B, fun n ↦ (dsum_le_csum a ha n).trans (hB n)⟩
  · rintro ⟨B, hB⟩
    refine ⟨a 1 + 2 * B, fun n ↦ ?_⟩
    have hB0 : 0 ≤ B := by simpa [dsum] using hB 0
    cases n with
    | zero => simpa [csum] using add_nonneg (hpos 1) (mul_nonneg zero_le_two hB0)
    | succ n => linarith [csum_le_dsum a ha n, hB n]

/-- Every partial sum `a 1 + ⋯ + a N` sits under a dyadic one, so the dyadic sums decide
boundedness of all partial sums. -/
theorem sum_Ico_le_dsum (a : ℕ → K) (hpos : ∀ k, 0 ≤ a k) (N : ℕ) :
    ∑ k ∈ Ico 1 (N + 1), a k ≤ dsum a N :=
  sum_le_sum_of_subset_of_nonneg (Ico_subset_Ico_right N.lt_two_pow_self) fun k _ _ ↦ hpos k

end Condensation

/-! ### Instances: rung 0 passes, one power past it fails -/

/-- The harmonic terms `1 / k`: the discrete `∫ ds / s`, rung 0. -/
def harm (k : ℕ) : ℚ := 1 / k

theorem harm_nonneg (k : ℕ) : 0 ≤ harm k := by
  unfold harm; positivity

theorem harm_anti : AntitoneOn harm (Set.Ici 1) := fun _ hi _ _ hik ↦
  one_div_le_one_div_of_le (Nat.cast_pos.2 hi) (Nat.cast_le.2 hik)

/-- One rung down, the harmonic terms become all ones: `2 ^ j * (1 / 2 ^ j) = 1`. -/
theorem csum_harm (n : ℕ) : csum harm n = n := by
  simp [csum, harm]

/-- Rung 0 passes: `1 + 1 / 2 + ⋯ + 1 / (2 ^ n - 1) ≥ n / 2`, the clock never stops. -/
theorem harmonic_grows (n : ℕ) : (n : ℚ) ≤ 2 * dsum harm n := by
  have h := csum_le_dsum harm harm_anti n
  rw [csum_harm] at h
  have h1 : harm 1 = 1 := by simp [harm]
  push_cast at h
  linarith

/-- The inverse squares `1 / k ^ 2`: one full power past rung 0. -/
def isq (k : ℕ) : ℚ := 1 / (k * k)

theorem isq_anti : AntitoneOn isq (Set.Ici 1) := fun i hi k _ hik ↦ by
  have hi' : (0 : ℚ) < i := Nat.cast_pos.2 hi
  have hik' : (i : ℚ) ≤ k := Nat.cast_le.2 hik
  exact one_div_le_one_div_of_le (mul_pos hi' hi') (mul_le_mul hik' hik' hi'.le (hi'.le.trans hik'))

/-- One rung down, the inverse squares become a geometric series:
`∑ j < n, 2 ^ j / 4 ^ j = 2 - 2 / 2 ^ n`. -/
theorem csum_isq (n : ℕ) : csum isq n = 2 - 2 * (1 / 2 ^ n) := by
  induction n with
  | zero => simp [csum]
  | succ n ih =>
    rw [csum_succ, ih, isq]
    push_cast
    field_simp
    ring

/-- One power past rung 0 fails the test (the sum converges):
`1 + 1 / 4 + 1 / 9 + ⋯ + 1 / (2 ^ n - 1) ^ 2 ≤ 2` for every `n`. -/
theorem inverse_squares_bounded (n : ℕ) : dsum isq n ≤ 2 := by
  have h := dsum_le_csum isq isq_anti n
  rw [csum_isq] at h
  have hp : (0 : ℚ) < 1 / 2 ^ n := by positivity
  linarith

/-! ### Direction decides: same speed, opposite outcomes -/

section Direction

variable {R : Type*} [CommRing R]

/-- First component of the rotation field `|x| ^ 2 • J x` in the plane. -/
def rot1 (x y : R) : R := (x ^ 2 + y ^ 2) * (-y)

/-- Second component of the rotation field `|x| ^ 2 • J x`. -/
def rot2 (x y : R) : R := (x ^ 2 + y ^ 2) * x

/-- First component of the radial field `|x| ^ 2 • x`. -/
def rad1 (x y : R) : R := (x ^ 2 + y ^ 2) * x

/-- Second component of the radial field `|x| ^ 2 • x`. -/
def rad2 (x y : R) : R := (x ^ 2 + y ^ 2) * y

/-- The rotation field has speed `r ^ 3` (speed squared `r ^ 6`): the size alone fails the Osgood
test. -/
theorem rot_speed (x y : R) : rot1 x y ^ 2 + rot2 x y ^ 2 = (x ^ 2 + y ^ 2) ^ 3 := by
  unfold rot1 rot2; ring

/-- The radial field has the same speed `r ^ 3`. -/
theorem rad_speed (x y : R) : rad1 x y ^ 2 + rad2 x y ^ 2 = (x ^ 2 + y ^ 2) ^ 3 := by
  unfold rad1 rad2; ring

/-- The rotation pushes sideways: `x · F = 0`, so `d/dt |x| ^ 2 = 2 * x · F = 0` and the size
never changes, though the speed fails the Osgood test. -/
theorem rot_perp (x y : R) : x * rot1 x y + y * rot2 x y = 0 := by
  unfold rot1 rot2; ring

/-- The radial field pushes straight out: `x · F = r ^ 4`, so `ρ = r ^ 2` obeys `ρ' = 2 * ρ ^ 2`,
the Riccati law, and blows up at time `1 / (2 * ρ 0)`. -/
theorem rad_out (x y : R) : x * rad1 x y + y * rad2 x y = (x ^ 2 + y ^ 2) ^ 2 := by
  unfold rad1 rad2; ring

end Direction

end MultAdd.OsgoodInverse
