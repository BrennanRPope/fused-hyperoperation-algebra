module

public import FusedAlgebra.Basic

/-!
# Thermodynamics as mult-add

Independent systems multiply their partition functions, so the excess `Z - 1` composes by
mult-add; efficiencies of engines in series compose by the weighted family at `k = -1`; Tsallis
entropies compose at `k = 1 - q`; and several two-walled laws (Fermi functions, Ising bonds in
parallel, boosts) are velocity addition `vel u v = (u + v) / (1 + u * v)`, which the lens
`lens t = (1 - t) / (1 + t)` turns into multiplication.

## Main definitions

* `ThermoFused.outer`: all products `xᵢ * yⱼ` of two lists (joint states of two systems).
* `ThermoFused.vel`: velocity addition, `(u + v) / (1 + u * v)`.
* `ThermoFused.lens`: the lens `(1 - t) / (1 + t)`.
* `ThermoFused.stokes`: Stokes' law for the reflectance of two plates.
* `ThermoFused.chain3`: the Lee–Yang polynomial of a three-spin chain.

## Main results

* `sum_outer`, `partition_excess`: a grid's total is its edges multiplied, so excesses mult-add.
* `unattainable`, `never_the_wall`: finitely many cooling steps short of the wall stay short.
* `engine_cascade`, `carnot_series`, `second_law_composes`: engines in series at weight `-1`.
* `tsallis_composes`, `qlog_wAdd`: Tsallis pseudo-additivity is `wAdd k`.
* `lens_vel`, `lens_involution`: the lens carries `vel` to multiplication and is an involution
  (Kramers–Wannier duality).
* `stokes_additive`, `shields`: `R / (1 - R)` turns Stokes' law into addition.
* `chain3_factor`, `chain3_disc`: the three-spin chain and its zeros.

Everything here is elementary algebra and known physics bookkeeping; what the file records is
that each rule is a statement about mult-add or its weighted family. Real `exp`, `log`, square
roots and real powers enter only as abstract quantities with the one property the proof uses (a
multiplicative map `f`, a number `s` with `s * s = r`). The core-Lean original worked over `ℚ`;
here the identities hold in any commutative ring or field and the inequalities in any linearly
ordered commutative ring or field.
-/

@[expose] public section

namespace ThermoFused

open MultAdd

variable {R : Type*} [CommRing R] {K : Type*} [Field K]

/-! ### Grids of products -/

/-- All products `xᵢ * yⱼ`: the states of two independent systems side by side. -/
def outer {M : Type*} [Mul M] (xs ys : List M) : List M := xs.flatMap fun x ↦ ys.map (x * ·)

@[simp]
theorem outer_nil {M : Type*} [Mul M] (ys : List M) : outer [] ys = [] := rfl

@[simp]
theorem outer_cons {M : Type*} [Mul M] (x : M) (xs ys : List M) :
    outer (x :: xs) ys = ys.map (x * ·) ++ outer xs ys := rfl

/-- A grid's total is its edges multiplied. -/
theorem sum_outer {S : Type*} [NonUnitalNonAssocSemiring S] (xs ys : List S) :
    (outer xs ys).sum = xs.sum * ys.sum := by
  induction xs with
  | nil => simp
  | cons x xs ih =>
    have h : (ys.map (x * ·)).sum = x * ys.sum := by
      clear ih; induction ys <;> simp_all [mul_add]
    simp [ih, h, add_mul]

/-- A multiplicative map carries a grid to the grid of the images. -/
theorem map_outer {M N F : Type*} [Mul M] [Mul N] [FunLike F M N] [MulHomClass F M N] (f : F)
    (xs ys : List M) : (outer xs ys).map f = outer (xs.map f) (ys.map f) := by
  induction xs with
  | nil => simp
  | cons x xs ih => simp [ih, Function.comp_def]

/-! ### Partition functions: the corner is "both excited" -/

/-- Independent systems: `Z_AB = Z_A * Z_B`, so `Z - 1` composes by mult-add. -/
theorem partition_square (ZA ZB : R) : ZA * ZB - 1 = (ZA - 1) ✳ (ZB - 1) := by
  unfold multAdd; ring

/-- The joint partition function of two systems is the product of theirs (the weights of joint
states are products of weights). -/
theorem partition_multiplies (wA wB : List R) : (outer wA wB).sum = wA.sum * wB.sum :=
  sum_outer wA wB

/-- Put each ground state first with weight `1`. Then the joint partition function is
`1 + (zA ✳ zB)`, where `z` is the sum of the excited weights: the corner `zA * zB` is "both
excited". -/
theorem partition_excess (xs ys : List R) :
    (outer (1 :: xs) (1 :: ys)).sum = 1 + (xs.sum ✳ ys.sum) := by
  rw [sum_outer, List.sum_cons, List.sum_cons, one_add_multAdd]

/-- Free energy (and any logarithm) is the linear coordinate: it adds. -/
theorem free_energy_adds {G : Type*} [Add G] (L : R → G) (hL : ∀ x y, L (x * y) = L x + L y)
    (zA zB : R) : L (1 + (zA ✳ zB)) = L (1 + zA) + L (1 + zB) := by
  rw [one_add_multAdd, hL]

/-! ### The third law: the identity, and the wall you cannot reach -/

/-- As `T → 0` every excess goes to `0`, the identity, and the identity composes to the
identity. -/
theorem third_law_identity : (0 : R) ✳ 0 = 0 := multAdd_zero 0

/-- Residual entropy: ground-state degeneracies compose by mult-add of `g - 1`. -/
theorem residual_degeneracy (g1 g2 : R) : g1 * g2 - 1 = (g1 - 1) ✳ (g2 - 1) :=
  partition_square g1 g2

/-- Pauling's ice: `W = (3 / 2) ^ N`, so `W - 1` is the mult-add power of `1 / 2`. -/
theorem pauling_ice [CharZero K] (n : ℕ) : 1 + mpow (1 / 2 : K) n = (3 / 2) ^ n := by
  rw [one_add_mpow]; norm_num

/-- Landauer: `n` bits have `2 ^ n` states, and `2 ^ n - 1` is the mult-add power of `1`. -/
theorem landauer_states (n : ℕ) : 1 + mpow (1 : R) n = 2 ^ n := by
  rw [one_add_mpow, one_add_one_eq_two]

/-- One cooling step that keeps `T` above `0` (ratio `1 + a > 0`), followed by another, still
keeps `T` above `0`. -/
theorem cooling_step_pos [LinearOrder R] [IsStrictOrderedRing R] (a b : R) (ha : 0 < 1 + a)
    (hb : 0 < 1 + b) : 0 < 1 + (a ✳ b) := by
  rw [one_add_multAdd]; exact mul_pos ha hb

/-- Unattainability: any finite list of cooling steps, each short of the wall, stays short of
it. -/
theorem unattainable [LinearOrder R] [IsStrictOrderedRing R] (l : List R)
    (h : ∀ x ∈ l, 0 < 1 + x) : 0 < 1 + multAddList l := by
  rw [one_add_multAddList]
  exact List.prod_pos (by simpa using h)

/-- `T = 0` is the wall `-1`, and mult-add only reaches the wall if some step is already on
it. -/
theorem never_the_wall [LinearOrder R] [IsStrictOrderedRing R] (l : List R)
    (h : ∀ x ∈ l, 0 < 1 + x) : multAddList l ≠ -1 := by
  intro heq
  have := unattainable l h
  rw [heq] at this
  simp at this

/-- Equipartition: `n` quadratic degrees of freedom are `n` half-steps of `T`. If
`s * s = T / T₀` then `1 + mpow (s - 1) (2 * m) = (T / T₀) ^ m`. -/
theorem equipartition (s r : R) (hs : s * s = r) (m : ℕ) :
    1 + mpow (s - 1) (2 * m) = r ^ m := by
  rw [one_add_mpow, add_sub_cancel, pow_mul, sq, hs]

/-! ### Fermions are the cube; bosons are the way back -/

/-- Fermion grand partition function: each mode is empty or full (Pauli), so
`Ξ = ∏ (1 + zᵢ) = 1 + (z₁ ✳ ⋯ ✳ zₙ)`. This is `MultAdd.one_add_multAddList`. -/
theorem fermion_lift (l : List R) : 1 + multAddList l = (l.map (1 + ·)).prod :=
  one_add_multAddList l

/-- Every many-fermion state is one corner of the cube (one subset of occupied modes). This is
`MultAdd.one_add_multAddList_eq_sum_sublists`. -/
theorem fermion_cube (l : List R) : 1 + multAddList l = (l.sublists'.map List.prod).sum :=
  one_add_multAddList_eq_sum_sublists l

/-- A boson mode: `Ξ_B = 1 / (1 - z) = 1 + multAddInv (-z)`. The excess `z / (1 - z)` is the way
back from a negative fermion. -/
theorem boson_mode (z : K) (h : 1 - z ≠ 0) : (1 + multAddInv (-z)) * (1 - z) = 1 := by
  rw [sub_eq_add_neg] at h ⊢
  rw [← inv_one_add _ h, inv_mul_cancel₀ h]

/-- Boson–fermion duality, one mode: `Ξ_B(z) * Ξ_F(-z) = 1`. -/
theorem boson_fermion_dual (z : K) (h : 1 - z ≠ 0) : multAddInv (-z) ✳ -z = 0 :=
  multAddInv_multAdd (-z) (by rwa [← sub_eq_add_neg])

/-- Fermi function: with Boltzmann odds `E` the occupation is `1 / (1 + E)` and
`t = (E - 1) / (E + 1)`. Adding energies multiplies odds, and then `t` composes by velocity
addition: this is the numerator. -/
theorem fermi_velocity_num (E1 E2 : R) :
    (E1 - 1) * (E2 + 1) + (E1 + 1) * (E2 - 1) = 2 * (E1 * E2 - 1) := by ring

/-- The denominator: `(t₁ + t₂) / (1 + t₁ * t₂) = (E₁ * E₂ - 1) / (E₁ * E₂ + 1)`. -/
theorem fermi_velocity_den (E1 E2 : R) :
    (E1 + 1) * (E2 + 1) + (E1 - 1) * (E2 - 1) = 2 * (E1 * E2 + 1) := by ring

/-- Velocity addition (the Fermi function, Ising bonds in parallel, boosts). -/
def vel (u v : K) : K := (u + v) / (1 + u * v)

/-- The wall `1` absorbs (occupation `n = 0`). -/
theorem vel_wall_pos (v : K) (h : 1 + v ≠ 0) : vel 1 v = 1 := by
  rw [vel, one_mul, div_self h]

/-- The wall `-1` absorbs (occupation `n = 1`). -/
theorem vel_wall_neg (v : K) (h : 1 - v ≠ 0) : vel (-1) v = -1 := by
  have h' : 1 + -1 * v ≠ 0 := by rwa [neg_one_mul, ← sub_eq_add_neg]
  rw [vel, div_eq_iff h']; ring

/-- Infinite temperature (`t = 0`, `n = 1 / 2`) is the identity, not a wall. -/
@[simp]
theorem vel_identity (v : K) : vel 0 v = v := by simp [vel]

/-! ### Engines, fridges, heat pumps -/

/-- Two engines in series, the second running on the first's waste heat: total work over heat in
is `wAdd (-1) η₁ η₂` (the failure-union rule). -/
theorem engine_cascade (Qh e1 e2 : R) :
    e1 * Qh + e2 * (Qh - e1 * Qh) = wAdd (-1) e1 e2 * Qh := by
  unfold wAdd; ring

/-- Carnot: `η(Th, Tc) = wAdd (-1) (η(Th, Tm)) (η(Tm, Tc))`; the lift `1 - η` is `Tc / Th`. -/
theorem carnot_series (h m c : K) (hm : m ≠ 0) :
    1 - c / h = wAdd (-1) (1 - m / h) (1 - c / m) := by
  unfold wAdd; field_simp; ring

/-- The `k = -1` wall is the perfect engine: one perfect stage makes all perfect. -/
theorem perfect_stage (e : R) : wAdd (-1) 1 e = 1 := by
  unfold wAdd; ring

/-- Two fridges in cascade (the second pumps the first's rejected heat): work per unit of cold
heat composes by plain mult-add. -/
theorem fridge_cascade (Qc u1 u2 : R) : u1 * Qc + u2 * (Qc + u1 * Qc) = (u1 ✳ u2) * Qc := by
  unfold multAdd; ring

/-- Running an engine backwards is the mult-add inverse: `η / (1 - η) = multAddInv (-η)`. -/
theorem reverse_is_inverse (eta : K) : eta / (1 - eta) = multAddInv (-eta) := by
  rw [multAddInv, neg_neg, sub_eq_add_neg]

/-- So `1 + u` is the reciprocal of the engine's lift `1 - η`. -/
theorem fridge_lift (eta : K) (h : 1 - eta ≠ 0) : (1 + eta / (1 - eta)) * (1 - eta) = 1 := by
  field_simp; ring

/-- Carnot fridge: `1 + u = Th / Tc`. -/
theorem carnot_fridge (h c : K) (hc : c ≠ 0) : 1 + (h - c) / c = h / c := by
  field_simp; ring

/-- Curzon–Ahlborn (efficiency at maximum power): with `s * s = Tc / Th`, `η_CA = 1 - s` and
`η_C = 1 - s * s`. Two CA steps make one Carnot step: CA is the half-step of the `k = -1`
climb. -/
theorem curzon_ahlborn_half (s : R) : wAdd (-1) (1 - s) (1 - s) = 1 - s * s := by
  unfold wAdd; ring

/-- CA engines in series also compose by `wAdd (-1)` (and so do Otto cycles). -/
theorem ratio_series (s1 s2 : R) : 1 - s1 * s2 = wAdd (-1) (1 - s1) (1 - s2) := by
  unfold wAdd; ring

/-- The second law composes: stages below their Carnot bounds (`rᵢ ≤ 1 - ηᵢ`, `0 ≤ rᵢ`) make a
cascade below the overall bound `1 - r₁ * r₂`. -/
theorem second_law_composes [LinearOrder R] [IsStrictOrderedRing R] (e1 e2 r1 r2 : R)
    (h1 : 0 ≤ r1) (h2 : 0 ≤ r2) (b1 : r1 ≤ 1 - e1) (b2 : r2 ≤ 1 - e2) :
    wAdd (-1) e1 e2 ≤ 1 - r1 * r2 := by
  rw [wAdd_neg_one]
  linarith [mul_le_mul b1 b2 h2 (h1.trans b1)]

/-! ### Tsallis: the `q`-sum is weighted mult-add with `k = 1 - q` -/

/-- Power sums of independent sources multiply (`f = (· ^ q)` is multiplicative). -/
theorem power_sum_multiplies {M F : Type*} [Mul M] [FunLike F M R] [MulHomClass F M R] (f : F)
    (xs ys : List M) : ((outer xs ys).map f).sum = (xs.map f).sum * (ys.map f).sum := by
  rw [map_outer, sum_outer]

/-- Tsallis pseudo-additivity: if `1 + k * S` is the power sum for `A`, for `B` and for the pair
(whose power sum is the product), then `S_AB = wAdd k S_A S_B`. -/
theorem tsallis_composes [IsDomain R] (k SA SB SAB PA PB : R) (hk : k ≠ 0)
    (hA : 1 + k * SA = PA) (hB : 1 + k * SB = PB) (hAB : 1 + k * SAB = PA * PB) :
    SAB = wAdd k SA SB := by
  apply mul_left_cancel₀ hk
  apply add_left_cancel (a := 1)
  rw [hAB, one_add_mul_wAdd, hA, hB]

/-- The `q`-logarithm `ln_q u = (u ^ k - 1) / k` obeys `ln_q (u * v) = wAdd k (ln_q u) (ln_q v)`.
Here `f` stands for `u ↦ u ^ k`. -/
theorem qlog_wAdd {M F : Type*} [Mul M] [FunLike F M K] [MulHomClass F M K] (f : F) (k : K)
    (u v : M) (hk : k ≠ 0) :
    (f (u * v) - 1) / k = wAdd k ((f u - 1) / k) ((f v - 1) / k) := by
  rw [map_mul]; unfold wAdd; field_simp; ring

/-- Rényi entropy `log (1 + k * S) / k` adds: the logarithm of the lift. -/
theorem renyi_adds {G : Type*} [Add G] (L : R → G) (hL : ∀ x y, L (x * y) = L x + L y)
    (k SA SB : R) : L (1 + k * wAdd k SA SB) = L (1 + k * SA) + L (1 + k * SB) := by
  rw [one_add_mul_wAdd, hL]

/-! ### The Ising chain: series is parity, parallel is velocity, duality is the lens -/

/-- Summing out the middle spin of two bonds (`X = exp (2 * K)`): the "agree" weight minus the
"disagree" weight is `(X₁ - 1) * (X₂ - 1)`. -/
theorem decimation_num (X1 X2 : R) : (X1 * X2 + 1) - (X1 + X2) = (X1 - 1) * (X2 - 1) := by ring

/-- Their sum is `(X₁ + 1) * (X₂ + 1)`, so `tanh K' = tanh K₁ * tanh K₂`. -/
theorem decimation_den (X1 X2 : R) : (X1 * X2 + 1) + (X1 + X2) = (X1 + 1) * (X2 + 1) := by ring

/-- In the chance `p` that two linked spins disagree (`tanh K = 1 - 2 * p`), bonds in series
compose at `k = -2`, the parity weight (Matsui's piling-up rule). -/
theorem ising_series (p1 p2 : R) :
    1 - 2 * wAdd (-2) p1 p2 = (1 - 2 * p1) * (1 - 2 * p2) := by
  unfold wAdd; ring

/-- The lens of velocity addition, `(1 - t) / (1 + t)`. -/
def lens (t : K) : K := (1 - t) / (1 + t)

/-- The numerator of the lens of a velocity sum. -/
theorem lens_num (u v : R) : (1 + u * v) - (u + v) = (1 - u) * (1 - v) := by ring

/-- The denominator of the lens of a velocity sum. -/
theorem lens_den (u v : R) : (1 + u * v) + (u + v) = (1 + u) * (1 + v) := by ring

/-- The lens turns velocity addition (bonds in parallel) into multiplication (bonds in
series). -/
theorem lens_vel (u v : K) (hu : 1 + u ≠ 0) (hv : 1 + v ≠ 0) (huv : 1 + u * v ≠ 0) :
    lens (vel u v) = lens u * lens v := by
  have e1 : 1 - (u + v) / (1 + u * v) = (1 - u) * (1 - v) / (1 + u * v) := by
    field_simp; ring
  have e2 : 1 + (u + v) / (1 + u * v) = (1 + u) * (1 + v) / (1 + u * v) := by
    field_simp; ring
  rw [lens, vel, e1, e2, lens, lens]
  field_simp

/-- Kramers–Wannier duality `t ↦ lens t` is an involution. -/
theorem lens_involution [NeZero (2 : K)] (t : K) (h : 1 + t ≠ 0) : lens (lens t) = t := by
  have h2 : (2 : K) ≠ 0 := two_ne_zero
  have e1 : 1 - (1 - t) / (1 + t) = 2 * t / (1 + t) := by field_simp; ring
  have e2 : 1 + (1 - t) / (1 + t) = 2 / (1 + t) := by field_simp; ring
  rw [lens, lens, e1, e2]
  field_simp

/-- The self-dual point (the 2D Ising critical point, `tanh K_c = √2 - 1`) is the fixed point of
the lens: any `t` with `t ^ 2 + 2 * t - 1 = 0`. -/
theorem lens_fixed_point (t : K) (h : 1 + t ≠ 0) (ht : t * t + 2 * t - 1 = 0) : lens t = t := by
  rw [lens, div_eq_iff h]; linear_combination -ht

/-- A Doppler factor `D ^ 2 = (1 + β) / (1 - β)` is the reciprocal of the lens, so temperature
ratios of boosted blackbody light multiply as boosts add. -/
theorem doppler_is_lens (bb : K) (h1 : 1 + bb ≠ 0) (h2 : 1 - bb ≠ 0) :
    lens bb * ((1 + bb) / (1 - bb)) = 1 := by
  unfold lens; field_simp

/-! ### Gray plates and Stokes' pile: a law on the border parabola -/

/-- Reflectance of two non-absorbing plates (Stokes 1862), equivalently `1` minus the effective
emissivity of two gray plates. -/
def stokes (a b : K) : K := (a + b - 2 * a * b) / (1 - a * b)

@[simp]
theorem stokes_identity (a : K) : stokes a 0 = a := by simp [stokes]

/-- A perfect mirror (`R = 1`) is the wall. -/
theorem stokes_wall (b : K) (h : 1 - b ≠ 0) : stokes 1 b = 1 := by
  rw [stokes, div_eq_iff (by simpa using h)]; ring

/-- The numerator of `N (stokes a b)`, where `N R = R / (1 - R)`. -/
theorem stokes_num (a b : R) : a * (1 - b) + b * (1 - a) = a + b - 2 * a * b := by ring

/-- The denominator of `N (stokes a b)`. -/
theorem stokes_den (a b : R) : (1 - a * b) - (a + b - 2 * a * b) = (1 - a) * (1 - b) := by ring

/-- The lens `N R = R / (1 - R)` turns Stokes' law into plain addition. -/
theorem stokes_additive (a b : K) (ha : 1 - a ≠ 0) (hb : 1 - b ≠ 0) (hab : 1 - a * b ≠ 0) :
    stokes a b / (1 - stokes a b) = a / (1 - a) + b / (1 - b) := by
  have e : 1 - (a + b - 2 * a * b) / (1 - a * b) = (1 - a) * (1 - b) / (1 - a * b) := by
    field_simp; ring
  rw [stokes, e]
  field_simp; ring

/-- The rate `1 + k * x + q * x ^ 2` at `(k, q) = (-2, 1)` is `(1 - x) ^ 2`: one double wall. -/
theorem stokes_double_wall (x : R) : 1 + (-2) * x + 1 * x * x = (1 - x) * (1 - x) := by ring

/-- The point `(k, q) = (-2, 1)` is on the border parabola `k ^ 2 = 4 * q`. -/
theorem stokes_on_border : (-2 : R) * (-2) - 4 * 1 = 0 := by ring

/-- `n` identical plates: `Rₙ = n * R / (1 + (n - 1) * R)` has `N Rₙ = n * N R`, so `n` radiation
shields cut the heat by `1 / (n + 1)`, not exponentially. -/
theorem shields (n r : K) (h2 : 1 + (n - 1) * r ≠ 0) :
    (n * r / (1 + (n - 1) * r)) / (1 - n * r / (1 + (n - 1) * r)) = n * (r / (1 - r)) := by
  have e : 1 - n * r / (1 + (n - 1) * r) = (1 - r) / (1 + (n - 1) * r) := by
    rw [eq_div_iff h2, sub_mul, div_mul_cancel₀ _ h2]; ring
  rw [e, div_div_div_cancel_right₀ h2, mul_div_assoc]

/-! ### Lee–Yang: bonds are weighted mult-add squares -/

/-- At infinite temperature (`a = 1`) a bond's Lee–Yang square is the plain mult-add square. -/
theorem bond_infinite_T (x y : R) : 1 + 1 * x + 1 * y + x * y = 1 + (x ✳ y) := by
  unfold multAdd; ring

/-- At zero temperature (`a = 0`) only the diagonal corners survive. -/
theorem bond_zero_T (x y : R) : 1 + 0 * x + 0 * y + x * y = 1 + x * y := by ring

/-- Two spins, ferromagnetic (`0 ≤ a < 1`): `Z = 1 + 2 * a * z + z ^ 2` has no real zero. -/
theorem two_spin_no_real_zero [LinearOrder R] [IsStrictOrderedRing R] (a z : R) (h0 : 0 ≤ a)
    (h1 : a < 1) : 0 < 1 + 2 * a * z + z * z := by
  nlinarith [mul_self_nonneg (z + a), mul_pos (sub_pos.2 h1) (add_pos_of_pos_of_nonneg one_pos h0)]

/-- Three spins in a chain with bond weights `a`, `b`, summed state by state (`z` counts
up-spins): no spin up; `{1}`, `{2}`, `{3}`; `{1, 2}`, `{2, 3}`, `{1, 3}`; all up. -/
def chain3 (a b z : R) : R :=
  1 + (a + a * b + b) * z + (b + a + a * b) * (z * z) + z * z * z

/-- The middle coefficients of the three-spin chain are `a ✳ b`. -/
theorem chain3_multAdd (a b z : R) :
    chain3 a b z = 1 + (a ✳ b) * z + (a ✳ b) * (z * z) + z * z * z := by
  unfold chain3 multAdd; ring

/-- The chain factors as `(1 + z) * (1 + (m - 1) * z + z ^ 2)` with `m = a ✳ b`. -/
theorem chain3_factor (m z : R) :
    1 + m * z + m * (z * z) + z * z * z = (1 + z) * (1 + (m - 1) * z + z * z) := by ring

/-- Ferromagnetic bonds (`0 ≤ a, b ≤ 1`) give `0 ≤ a ✳ b ≤ 3`. -/
theorem chain3_ferro_bound [LinearOrder R] [IsStrictOrderedRing R] (a b : R) (ha0 : 0 ≤ a)
    (ha1 : a ≤ 1) (hb0 : 0 ≤ b) (hb1 : b ≤ 1) : 0 ≤ a ✳ b ∧ a ✳ b ≤ 3 := by
  unfold multAdd
  constructor
  · positivity
  · nlinarith [mul_le_mul ha1 hb1 hb0 zero_le_one]

/-- So the quadratic's discriminant `(m - 1) ^ 2 - 4 = (m - 3) * (m + 1)` is at most `0`: its
zeros lie on the unit circle, as Lee–Yang says. -/
theorem chain3_disc [LinearOrder R] [IsStrictOrderedRing R] (a b : R) (ha0 : 0 ≤ a)
    (ha1 : a ≤ 1) (hb0 : 0 ≤ b) (hb1 : b ≤ 1) : (a ✳ b - 1) * (a ✳ b - 1) - 4 ≤ 0 := by
  obtain ⟨l, u⟩ := chain3_ferro_bound a b ha0 ha1 hb0 hb1
  nlinarith

/-! ### Wick rotation is a change of units by `i` -/

/-- The law at `(k, q)` read in units `i * a` is the law at `(i * k, i ^ 2 * q)`: this is the
numerator, which holds for every `i`. -/
theorem wick_num (i k a b : R) :
    i * a + i * b + k * (i * a) * (i * b) = i * (a + b + (i * k) * a * b) := by ring

/-- The denominator, where `i * i = -1` turns `q` into `-q`. On the `k = 0` axis thermal `tanh`
addition (two walls) becomes quantum `tan` addition (a turn). -/
theorem wick_den (i q a b : R) (hi : i * i = -1) :
    1 - q * (i * a) * (i * b) = 1 - (-q) * a * b := by
  linear_combination -(q * a * b) * hi

/-! ### Fluctuation theorems and chemistry -/

/-- Jarzynski: for independent processes the averages `⟨exp (-β * W)⟩` multiply (weights times
probabilities form a grid), so `⟨exp (-β * W)⟩ - 1` mult-adds. -/
theorem jarzynski_composes (xs ys : List R) :
    (outer xs ys).sum - 1 = (xs.sum - 1) ✳ (ys.sum - 1) := by
  rw [sum_outer, partition_square]

/-- Hess's law: coupled reactions multiply equilibrium constants. -/
theorem hess (K1 K2 : R) : K1 * K2 - 1 = (K1 - 1) ✳ (K2 - 1) := partition_square K1 K2

end ThermoFused
