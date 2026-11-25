/-
  Thermodynamics in Fused Hyperoperation Algebra
  ==============================================

  Companion page: ../index.html.  Numbers: ../checks.py (55 checks).

  STATUS: core Lean 4 only (no mathlib), no `sorry`, written for the shared
  library (`import FusedAlgebra.Base`, Lean 4.37.0-rc1) in the style of its
  compiled topics (`Rat`, `grind`, `Rat.mul_nonneg`, `Rat.mul_pos`).
  NOT YET COMPILED where it was written: Lean is not installed in that
  workspace.  Every ring identity here was also checked with sympy.
  To compile: copy to lean-library/FusedAlgebra/Topics/Thermodynamics.lean,
  add it to build.sh and FusedAlgebra.lean, run ./build.sh.

  The circle operation (called mult-add before 1.0) a ◯ b = a + b + ab, its lift 1 + (a ◯ b) = (1 + a)(1 + b), the
  weighted family circK k a b = a + b + k·ab with 1 + k(a ◯[k] b) = (1+ka)(1+kb),
  the inverse −u/(1+u), list folds and the cube's corners all come from Base.

  What Lean does not check here: anything needing real exp, log, sqrt or
  real powers (Boltzmann factors, tanh, x^q, √(Tc/Th)).  Those enter as
  abstract quantities with the one property the proof needs (for example
  "s·s = Tc/Th", or "f(xy) = f(x)f(y)").  The physics models themselves
  (that a cascade is energy bookkeeping, that a spin sum is a partition
  function) are checked numerically in ../checks.py.
-/
import FusedAlgebra.Base

namespace ThermoFused

open CircleOp

/-! ## 0. Two small helpers -/

theorem sq_nonneg' (s : Rat) : 0 ≤ s * s := by
  rcases (@Rat.le_total 0 s) with h | h
  · exact Rat.mul_nonneg h h
  · have h1 : 0 ≤ -s := by grind
    have := Rat.mul_nonneg h1 h1
    grind

/-- All products xᵢ·yⱼ: the states of two independent systems side by side. -/
def outer : List Rat → List Rat → List Rat
  | [], _ => []
  | x :: xs, ys => ys.map (fun y => x * y) ++ outer xs ys

theorem sumL_map_mul (x : Rat) (ys : List Rat) :
    sumL (ys.map (fun y => x * y)) = x * sumL ys := by
  induction ys with
  | nil => simp only [List.map_nil, sumL]; grind
  | cons y ys ih => simp only [List.map_cons, sumL]; rw [ih]; grind

/-- A grid's total is its edges multiplied. -/
theorem outer_sum (xs ys : List Rat) : sumL (outer xs ys) = sumL xs * sumL ys := by
  induction xs with
  | nil => simp only [outer, sumL]; grind
  | cons x xs ih => simp only [outer, sumL]; rw [sumL_append, sumL_map_mul, ih]; grind

/-! ## 1. Partition functions: the corner is "both excited" -/

/-- Independent systems: Z_AB = Z_A · Z_B, so Z − 1 composes by the circle operation. -/
theorem partition_square (ZA ZB : Rat) : ZA * ZB - 1 = circ (ZA - 1) (ZB - 1) := by
  unfold circ; grind

/-- The joint partition function of two systems is the product of theirs
    (the weights of joint states are products of weights). -/
theorem partition_multiplies (wA wB : List Rat) :
    sumL (outer wA wB) = sumL wA * sumL wB := outer_sum wA wB

/-- Put each ground state first with weight 1 (energies measured from it).
    Then the joint partition function is 1 + (zA ◯ zB), where z is the sum of
    the excited weights: zA is "A excited, B in its ground state", zB the
    reverse, and the corner zA·zB is "both excited". -/
theorem partition_excess (xs ys : List Rat) :
    sumL (outer (1 :: xs) (1 :: ys)) = 1 + circ (sumL xs) (sumL ys) := by
  rw [outer_sum]; simp only [sumL]; rw [circ_lift]

/-- Free energy (and any log) is the linear coordinate: it adds. -/
theorem free_energy_adds {G : Type} [Add G] (L : Rat → G)
    (hL : ∀ x y, L (x * y) = L x + L y) (zA zB : Rat) :
    L (1 + circ zA zB) = L (1 + zA) + L (1 + zB) := by
  rw [circ_lift, hL]

/-! ## 2. The third law: the identity, and the wall you can't reach -/

/-- As T → 0 every excess goes to 0, the identity, and the identity composes
    to the identity: a composite's entropy also goes to 0. -/
theorem third_law_identity : circ (0 : Rat) 0 = 0 := circ_zero 0

/-- Residual entropy: ground-state degeneracies compose by the circle operation of g − 1. -/
theorem residual_degeneracy (g1 g2 : Rat) : g1 * g2 - 1 = circ (g1 - 1) (g2 - 1) :=
  partition_square g1 g2

/-- Pauling's ice: W = (3/2)^N, so W − 1 is the circle power (1/2)^(◯N). -/
theorem pauling_ice (n : Nat) : 1 + mpow (1/2 : Rat) n = (3/2) ^ n := by
  rw [mpow_lift]
  have e : (1 + 1/2 : Rat) = 3/2 := by grind
  rw [e]

/-- Landauer: n bits have 2^n states, 2^n − 1 = 1^(◯n). -/
theorem landauer_states (n : Nat) : 1 + mpow (1 : Rat) n = 2 ^ n := by
  rw [mpow_lift]
  have e : (1 + 1 : Rat) = 2 := by grind
  rw [e]

/-- One cooling step that keeps T above 0 (ratio 1 + a > 0), followed by
    another, still keeps T above 0. -/
theorem cooling_step_pos (a b : Rat) (ha : 0 < 1 + a) (hb : 0 < 1 + b) :
    0 < 1 + circ a b := by
  rw [circ_lift]; exact Rat.mul_pos ha hb

/-- Unattainability: any finite list of cooling steps, each short of the wall,
    stays short of it. T = 0 is the wall −1, and the circle operation only reaches the
    wall if some step is already on it. -/
theorem unattainable (l : List Rat) (h : ∀ x ∈ l, 0 < 1 + x) : 0 < 1 + foldCirc l := by
  induction l with
  | nil => simp only [foldCirc]; grind
  | cons a l ih =>
    simp only [foldCirc]
    rw [circ_lift]
    exact Rat.mul_pos (h a (by simp)) (ih (fun x hx => h x (by simp [hx])))

theorem never_the_wall (l : List Rat) (h : ∀ x ∈ l, 0 < 1 + x) : foldCirc l ≠ -1 := by
  intro heq
  have := unattainable l h
  rw [heq] at this
  grind

/-- Equipartition: n quadratic degrees of freedom are n half-steps of T.
    If s·s = T/T₀ then 1 + (s − 1)^(◯2m) = (T/T₀)^m. -/
theorem equipartition (s r : Rat) (hs : s * s = r) (m : Nat) :
    1 + mpow (s - 1) (2 * m) = r ^ m := by
  rw [mpow_lift]
  have e : (1 + (s - 1)) = s := by grind
  rw [e, ← pw_eq_pow, ← pw_eq_pow, pw_mul]
  have e2 : pw s 2 = r := by
    show pw s 0 * s * s = r
    simp only [pw]; grind
  rw [e2]

/-! ## 3. Fermions are the cube; bosons are the way back -/

/-- Fermion grand partition function: each mode empty or full (Pauli), so
    Ξ = Π(1 + zᵢ) = 1 + (z₁ ◯ … ◯ zₙ). -/
theorem fermion_lift (l : List Rat) : 1 + foldCirc l = prodL (l.map (fun x => 1 + x)) :=
  foldCirc_lift l

/-- …and that is the sum over every corner of the cube: every many-fermion
    state is one corner (one subset of occupied modes). -/
theorem fermion_cube (l : List Rat) : 1 + foldCirc l = sumL ((groups l).map prodL) :=
  one_add_foldCirc_eq_corners l

/-- A boson mode: Ξ_B = 1/(1 − z) = 1 + circInv(−z). The excess z/(1 − z)
    is the way back from a negative fermion. -/
theorem boson_mode (z : Rat) (h : 1 + -z ≠ 0) : (1 + circInv (-z)) * (1 - z) = 1 := by
  unfold circInv; grind

/-- Boson–fermion duality, one mode: Ξ_B(z) · Ξ_F(−z) = 1. -/
theorem boson_fermion_dual (z : Rat) (h : 1 + -z ≠ 0) :
    circ (circInv (-z)) (-z) = 0 := circInv_circ (-z) h

/-- Fermi function: with Boltzmann odds E = e^x the occupation is 1/(1 + E)
    and t = 1 − 2n = (E − 1)/(E + 1).  Adding energies multiplies odds, and
    then t composes by velocity addition: numerator… -/
theorem fermi_velocity_num (E1 E2 : Rat) :
    (E1 - 1) * (E2 + 1) + (E1 + 1) * (E2 - 1) = 2 * (E1 * E2 - 1) := by grind

/-- …and denominator: (t₁ + t₂)/(1 + t₁t₂) = (E₁E₂ − 1)/(E₁E₂ + 1). -/
theorem fermi_velocity_den (E1 E2 : Rat) :
    (E1 + 1) * (E2 + 1) + (E1 - 1) * (E2 - 1) = 2 * (E1 * E2 + 1) := by grind

/-- Velocity addition (the Fermi function, Ising bonds in parallel, boosts). -/
def vel (u v : Rat) : Rat := (u + v) / (1 + u * v)

/-- The two walls ±1 (n = 0 and n = 1, i.e. T = 0⁺ and T = 0⁻) absorb. -/
theorem vel_wall_pos (v : Rat) (h : 1 + v ≠ 0) : vel 1 v = 1 := by
  unfold vel; grind

theorem vel_wall_neg (v : Rat) (h : 1 - v ≠ 0) : vel (-1) v = -1 := by
  unfold vel; grind

/-- Infinite temperature (t = 0, n = 1/2) is the identity, not a wall. -/
theorem vel_identity (v : Rat) : vel 0 v = v := by
  unfold vel; grind

/-! ## 4. Engines, fridges, heat pumps -/

/-- Two engines in series, the second running on the first's waste heat:
    total work / heat in = η₁ ◯[−1] η₂ (the failure-union rule). -/
theorem engine_cascade (Qh e1 e2 : Rat) :
    e1 * Qh + e2 * (Qh - e1 * Qh) = circK (-1) e1 e2 * Qh := by
  unfold circK; grind

/-- Carnot: η(Th,Tc) = η(Th,Tm) ◯[−1] η(Tm,Tc); the lift 1 − η is Tc/Th. -/
theorem carnot_series (h m c : Rat) (hh : h ≠ 0) (hm : m ≠ 0) :
    1 - c / h = circK (-1) (1 - m / h) (1 - c / m) := by
  unfold circK; grind

/-- The k = −1 wall is the perfect engine: one perfect stage makes all perfect. -/
theorem perfect_stage (e : Rat) : circK (-1) 1 e = 1 := by
  unfold circK; grind

/-- Two fridges in cascade (the second pumps the first's rejected heat):
    work per unit of cold heat u = W/Qc composes by the plain circle operation. -/
theorem fridge_cascade (Qc u1 u2 : Rat) :
    u1 * Qc + u2 * (Qc + u1 * Qc) = circ u1 u2 * Qc := by
  unfold circ; grind

/-- Running an engine backwards is the circle inverse:
    u = η/(1 − η) = circInv(−η). -/
theorem reverse_is_inverse (eta : Rat) (h : 1 - eta ≠ 0) :
    eta / (1 - eta) = circInv (-eta) := by
  unfold circInv; grind

/-- …so 1 + u is the reciprocal of the engine's lift 1 − η. -/
theorem fridge_lift (eta : Rat) (h : 1 - eta ≠ 0) :
    (1 + eta / (1 - eta)) * (1 - eta) = 1 := by grind

/-- Carnot fridge: 1 + u = Th/Tc. -/
theorem carnot_fridge (h c : Rat) (hc : c ≠ 0) : 1 + (h - c) / c = h / c := by grind

/-- Curzon–Ahlborn (efficiency at maximum power): with s·s = Tc/Th,
    η_CA = 1 − s and η_C = 1 − s·s.  Two CA steps make one Carnot step:
    η_C = η_CA ◯[−1] η_CA.  CA is the half-step of the k = −1 climb. -/
theorem curzon_ahlborn_half (s : Rat) : circK (-1) (1 - s) (1 - s) = 1 - s * s := by
  unfold circK; grind

/-- CA engines in series also compose by ◯[−1] (and so do Otto cycles, whose
    1 − η = r^(1−γ) multiplies with the compression ratios). -/
theorem ratio_series (s1 s2 : Rat) : 1 - s1 * s2 = circK (-1) (1 - s1) (1 - s2) := by
  unfold circK; grind

/-- The second law composes: stages below their Carnot bounds (1 − ηᵢ ≥ rᵢ,
    rᵢ = T_out/T_in ≥ 0) make a cascade below the overall bound 1 − r₁r₂. -/
theorem second_law_composes (e1 e2 r1 r2 : Rat) (h1 : 0 ≤ r1) (h2 : 0 ≤ r2)
    (b1 : r1 ≤ 1 - e1) (b2 : r2 ≤ 1 - e2) : circK (-1) e1 e2 ≤ 1 - r1 * r2 := by
  have hA := Rat.mul_nonneg (show (0 : Rat) ≤ 1 - e1 - r1 by grind) (show (0 : Rat) ≤ 1 - e2 by grind)
  have hB := Rat.mul_nonneg h1 (show (0 : Rat) ≤ 1 - e2 - r2 by grind)
  unfold circK; grind

/-! ## 5. Tsallis: the q-sum is the circle operation with k = 1 − q -/

theorem map_mul_map (f : Rat → Rat) (hf : ∀ u v, f (u * v) = f u * f v) (x : Rat)
    (ys : List Rat) :
    (ys.map (fun y => x * y)).map f = (ys.map f).map (fun y => f x * y) := by
  induction ys with
  | nil => rfl
  | cons y ys ih => simp only [List.map_cons]; rw [hf, ih]

theorem outer_map (f : Rat → Rat) (hf : ∀ u v, f (u * v) = f u * f v) (xs ys : List Rat) :
    (outer xs ys).map f = outer (xs.map f) (ys.map f) := by
  induction xs with
  | nil => rfl
  | cons x xs ih => simp only [outer, List.map_cons, List.map_append]; rw [map_mul_map f hf, ih]

/-- Power sums of independent sources multiply (f = x ↦ x^q is multiplicative). -/
theorem power_sum_multiplies (f : Rat → Rat) (hf : ∀ u v, f (u * v) = f u * f v)
    (xs ys : List Rat) :
    sumL ((outer xs ys).map f) = sumL (xs.map f) * sumL (ys.map f) := by
  rw [outer_map f hf, outer_sum]

/-- Tsallis pseudo-additivity: if 1 + k·S is the power sum for A, for B and for
    the pair (whose power sum is the product), then S_AB = S_A ◯[k] S_B. -/
theorem tsallis_composes (k SA SB SAB PA PB : Rat) (hk : k ≠ 0)
    (hA : 1 + k * SA = PA) (hB : 1 + k * SB = PB) (hAB : 1 + k * SAB = PA * PB) :
    SAB = circK k SA SB := by
  have := circK_lift k SA SB
  have e : k * SAB = k * circK k SA SB := by grind
  have e2 : SAB = (k * SAB) / k := by grind
  rw [e2, e]; grind

/-- The q-logarithm ln_q u = (u^k − 1)/k obeys ln_q(uv) = ln_q u ◯[k] ln_q v. -/
theorem qlog_qsum (f : Rat → Rat) (hf : ∀ u v, f (u * v) = f u * f v) (k u v : Rat)
    (hk : k ≠ 0) :
    (f (u * v) - 1) / k = circK k ((f u - 1) / k) ((f v - 1) / k) := by
  rw [hf]; unfold circK; grind

/-- Rényi entropy log(1 + k·S)/k adds: the log of the lift. -/
theorem renyi_adds {G : Type} [Add G] (L : Rat → G) (hL : ∀ x y, L (x * y) = L x + L y)
    (k SA SB : Rat) : L (1 + k * circK k SA SB) = L (1 + k * SA) + L (1 + k * SB) := by
  rw [circK_lift, hL]

/-! ## 6. The Ising chain: series is parity, parallel is velocity, duality is the lens -/

/-- Summing out the middle spin of two bonds (X = e^{2K}): the "agree" weight
    minus the "disagree" weight is (X₁ − 1)(X₂ − 1)… -/
theorem decimation_num (X1 X2 : Rat) : (X1 * X2 + 1) - (X1 + X2) = (X1 - 1) * (X2 - 1) := by
  grind

/-- …and their sum is (X₁ + 1)(X₂ + 1), so tanh K' = tanh K₁ · tanh K₂. -/
theorem decimation_den (X1 X2 : Rat) : (X1 * X2 + 1) + (X1 + X2) = (X1 + 1) * (X2 + 1) := by
  grind

/-- In the chance p that two linked spins disagree (tanh K = 1 − 2p), bonds in
    series compose at k = −2, the parity weight (Matsui's piling-up rule). -/
theorem ising_series (p1 p2 : Rat) : 1 - 2 * circK (-2) p1 p2 = (1 - 2 * p1) * (1 - 2 * p2) := by
  unfold circK; grind

/-- The lens of velocity addition, M(t) = (1 − t)/(1 + t). -/
def lens (t : Rat) : Rat := (1 - t) / (1 + t)

/-- The lens turns velocity addition (bonds in parallel) into multiplication
    (bonds in series): polynomial form. -/
theorem lens_num (u v : Rat) : (1 + u * v) - (u + v) = (1 - u) * (1 - v) := by grind
theorem lens_den (u v : Rat) : (1 + u * v) + (u + v) = (1 + u) * (1 + v) := by grind

theorem lens_vel (u v : Rat) (hu : 1 + u ≠ 0) (hv : 1 + v ≠ 0) (huv : 1 + u * v ≠ 0) :
    lens (vel u v) = lens u * lens v := by
  unfold lens vel
  have e1 : 1 - (u + v) / (1 + u * v) = (1 - u) * (1 - v) / (1 + u * v) := by grind
  have e2 : 1 + (u + v) / (1 + u * v) = (1 + u) * (1 + v) / (1 + u * v) := by grind
  rw [e1, e2]; grind

/-- Kramers–Wannier duality t ↦ e^{−2K} = M(tanh K) is an involution. -/
theorem lens_involution (t : Rat) (h : 1 + t ≠ 0) : lens (lens t) = t := by
  unfold lens
  have e1 : 1 - (1 - t) / (1 + t) = 2 * t / (1 + t) := by grind
  have e2 : 1 + (1 - t) / (1 + t) = 2 / (1 + t) := by grind
  rw [e1, e2]; grind

/-- The self-dual point (2D Ising critical point, tanh K_c = √2 − 1) is the
    fixed point of the lens: any t with t² + 2t − 1 = 0. -/
theorem lens_fixed_point (t : Rat) (h : 1 + t ≠ 0) (ht : t * t + 2 * t - 1 = 0) :
    lens t = t := by
  unfold lens; grind

/-- A Doppler factor D² = (1 + β)/(1 − β) is the reciprocal of the lens, so
    temperature ratios of boosted blackbody light multiply as boosts add. -/
theorem doppler_is_lens (bb : Rat) (h1 : 1 + bb ≠ 0) (h2 : 1 - bb ≠ 0) :
    lens bb * ((1 + bb) / (1 - bb)) = 1 := by
  unfold lens; grind

/-! ## 7. Gray plates and Stokes' pile: a law on the border parabola -/

/-- Reflectance of two non-absorbing plates (Stokes 1862), equivalently
    1 − effective emissivity of two gray plates. -/
def stokes (a b : Rat) : Rat := (a + b - 2 * a * b) / (1 - a * b)

theorem stokes_identity (a : Rat) : stokes a 0 = a := by
  unfold stokes; grind

/-- A perfect mirror (R = 1) is the wall. -/
theorem stokes_wall (b : Rat) (h : 1 - b ≠ 0) : stokes 1 b = 1 := by
  unfold stokes; grind

/-- The lens N(R) = R/(1 − R) = 1/ε − 1: polynomial pieces… -/
theorem stokes_num (a b : Rat) : a * (1 - b) + b * (1 - a) = a + b - 2 * a * b := by grind
theorem stokes_den (a b : Rat) : (1 - a * b) - (a + b - 2 * a * b) = (1 - a) * (1 - b) := by grind

/-- …so N turns the law into plain addition (addition in disguise). -/
theorem stokes_additive (a b : Rat) (ha : 1 - a ≠ 0) (hb : 1 - b ≠ 0) (hab : 1 - a * b ≠ 0) :
    stokes a b / (1 - stokes a b) = a / (1 - a) + b / (1 - b) := by
  unfold stokes
  have e : 1 - (a + b - 2 * a * b) / (1 - a * b) = (1 - a) * (1 - b) / (1 - a * b) := by grind
  rw [e]; grind

/-- The rate 1 + k·x + q·x² at (k, q) = (−2, 1) is (1 − x)²: one double wall. -/
theorem stokes_double_wall (x : Rat) : 1 + (-2) * x + 1 * x * x = (1 - x) * (1 - x) := by grind

theorem stokes_on_border : (-2 : Rat) * (-2) - 4 * 1 = 0 := by grind

/-- n identical plates: Rₙ = nR/(1 + (n − 1)R) has N(Rₙ) = n·N(R), so n
    radiation shields cut the heat by 1/(n + 1), not exponentially. -/
theorem shields (n r : Rat) (h1 : 1 - r ≠ 0) (h2 : 1 + (n - 1) * r ≠ 0) :
    (n * r / (1 + (n - 1) * r)) / (1 - n * r / (1 + (n - 1) * r)) = n * (r / (1 - r)) := by
  have e : 1 - n * r / (1 + (n - 1) * r) = (1 - r) / (1 + (n - 1) * r) := by grind
  rw [e]; grind

/-! ## 8. Lee–Yang: bonds are circle-operation squares with dial k -/

/-- At infinite temperature (a = 1) a bond's Lee–Yang square is the plain
    circle-operation square. -/
theorem bond_infinite_T (x y : Rat) : 1 + 1 * x + 1 * y + x * y = 1 + circ x y := by
  unfold circ; grind

/-- At zero temperature (a = 0) only the diagonal corners survive. -/
theorem bond_zero_T (x y : Rat) : 1 + 0 * x + 0 * y + x * y = 1 + x * y := by grind

/-- Two spins, ferromagnetic (0 ≤ a < 1): Z = 1 + 2az + z² has no real zero. -/
theorem two_spin_no_real_zero (a z : Rat) (h0 : 0 ≤ a) (h1 : a < 1) :
    0 < 1 + 2 * a * z + z * z := by
  have hs := sq_nonneg' (z + a)
  have hp := Rat.mul_pos (show (0 : Rat) < 1 - a by grind) (show (0 : Rat) < 1 + a by grind)
  grind

/-- Three spins in a chain with bond weights a, b (a = e^{−2βJ}), summed state
    by state (z counts up-spins).  The middle coefficients are a ◯ b. -/
def chain3 (a b z : Rat) : Rat :=
  1                                       -- no spin up
  + (a + a * b + b) * z                   -- {1}, {2}, {3}
  + (b + a + a * b) * (z * z)             -- {1,2}, {2,3}, {1,3}
  + z * z * z                             -- all up

theorem chain3_circ (a b z : Rat) :
    chain3 a b z = 1 + circ a b * z + circ a b * (z * z) + z * z * z := by
  unfold chain3 circ; grind

/-- It factors as (1 + z)(1 + (m − 1)z + z²) with m = a ◯ b. -/
theorem chain3_factor (m z : Rat) :
    1 + m * z + m * (z * z) + z * z * z = (1 + z) * (1 + (m - 1) * z + z * z) := by grind

/-- Ferromagnetic bonds (0 ≤ a, b ≤ 1) give 0 ≤ a ◯ b ≤ 3, so the quadratic's
    discriminant (m − 1)² − 4 = (m − 3)(m + 1) ≤ 0: its zeros lie on the unit
    circle (the turning country), as Lee–Yang says. -/
theorem chain3_ferro_bound (a b : Rat) (ha0 : 0 ≤ a) (ha1 : a ≤ 1) (hb0 : 0 ≤ b) (hb1 : b ≤ 1) :
    0 ≤ circ a b ∧ circ a b ≤ 3 := by
  have h1 := Rat.mul_nonneg ha0 hb0
  have h2 := Rat.mul_nonneg (show (0 : Rat) ≤ 1 - a by grind) (show (0 : Rat) ≤ 1 + b by grind)
  unfold circ
  constructor <;> grind

theorem chain3_disc (a b : Rat) (ha0 : 0 ≤ a) (ha1 : a ≤ 1) (hb0 : 0 ≤ b) (hb1 : b ≤ 1) :
    (circ a b - 1) * (circ a b - 1) - 4 ≤ 0 := by
  obtain ⟨l, u⟩ := chain3_ferro_bound a b ha0 ha1 hb0 hb1
  have h := Rat.mul_nonneg (show (0 : Rat) ≤ 3 - circ a b by grind)
    (show (0 : Rat) ≤ circ a b + 1 by grind)
  grind

/-! ## 9. Wick rotation is a change of units by i -/

/-- With i·i = −1, the law at (k, q) read in units i·a is the law at (i·k, −q):
    numerator… -/
theorem wick_num (i k a b : Rat) (hi : i * i = -1) :
    i * a + i * b + k * (i * a) * (i * b) = i * (a + b + (i * k) * a * b) := by grind

/-- …and denominator.  On the k = 0 axis thermal tanh addition (two walls)
    becomes quantum tan addition (a turn). -/
theorem wick_den (i q a b : Rat) (hi : i * i = -1) :
    1 - q * (i * a) * (i * b) = 1 - (-q) * a * b := by grind

/-! ## 10. Fluctuation theorems and chemistry -/

/-- Jarzynski: for independent processes the averages ⟨e^{−βW}⟩ multiply
    (weights times probabilities form a grid), so ⟨e^{−βW}⟩ − 1 combines by the circle operation. -/
theorem jarzynski_composes (xs ys : List Rat) :
    sumL (outer xs ys) - 1 = circ (sumL xs - 1) (sumL ys - 1) := by
  rw [outer_sum]; unfold circ; grind

/-- Hess's law: coupled reactions multiply equilibrium constants. -/
theorem hess (K1 K2 : Rat) : K1 * K2 - 1 = circ (K1 - 1) (K2 - 1) := partition_square K1 K2

end ThermoFused
