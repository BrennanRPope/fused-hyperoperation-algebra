/-
  The dial's potential: the exact facts behind the search
  ========================================================

  Companion: ../search.py (scores 71 guesses against measured numbers). Core Lean 4.37.0-rc1, no mathlib,
  no `sorry`. Imports Base only.

  What is proved here:
  1. The dial's two exact symmetries.  Flip: the dial −k is the dial k seen through a ↦ −a.  Rescaling:
     the dial k/s is the dial k seen through a ↦ s·a, so every non-zero dial is every other one.
  2. So a potential for the dial must be even (flip), and if it also respects rescaling it is a single
     power of k.  The double well (hump at 0, dips at ±k0) is even but cannot respect rescaling: some
     outside scale has to set k0.  That is the exact reason the search needed a coupling to something with
     units.
  3. The mirror well: k ✳ (−k) = −k², so (k ✳ (−k) + k0²)² = (k0² − k²)², the 1D sombrero, with its
     expansion at the dip (curvature 8k0²) and, in the usual normalisation, m² = 2λv² and
     hump height = m²v²/8.
  4. The rolling dial today: at the exponential potential's late-time fixed point, Ω_φ = 1 and
     w = λ²/3 − 1; for the conformal slopes λ² = 2/3 and 8/3 that is −7/9 and −1/9.
  5. Two guesses ruled out by exact arithmetic on the measured masses (units of 10 MeV): the
     Coleman–Weinberg Higgs has negative mass squared, and the Veltman mass is more than twice the
     measured Higgs mass.
  6. The conformal reading's held-out test: the factor that gives the slope √(2/3) also couples to
     matter with α0² = 1/6, which gives PPN γ = 5/7.
-/
import FusedAlgebra.Base

namespace DialPotential
open MultAdd
open Lean.Grind

set_option linter.unusedSectionVars false

/-! ## 1. The dial's exact symmetries -/

section Ring
variable {R : Type} [CommRing R]

/-- Flip: the dial −k is the dial k seen through a ↦ −a. -/
theorem dial_flip (k a b : R) : wAdd (-k) a b = -(wAdd k (-a) (-b)) := by
  unfold wAdd; grind

/-- Rescaling: if k'·s = k, the dial k' on (s·a, s·b) is s times the dial k on (a, b). -/
theorem dial_rescale (k k' s a b : R) (h : k' * s = k) :
    wAdd k' (s * a) (s * b) = s * wAdd k a b := by
  unfold wAdd; rw [← h]; grind

/-- At weight zero the dial is plain addition: the symmetric point. -/
theorem dial_zero (a b : R) : wAdd 0 a b = a + b := by
  unfold wAdd; grind

/-! ## 3. The mirror well -/

/-- A dial against its mirror: k ✳ (−k) = −k². -/
theorem multAdd_mirror (k : R) : multAdd k (-k) = -(k * k) := by
  unfold multAdd; grind

/-- The mirror well: the square of (k ✳ (−k) + k0²). -/
def well (k0 k : R) : R := (multAdd k (-k) + k0 * k0) * (multAdd k (-k) + k0 * k0)

theorem well_eq (k0 k : R) : well k0 k = (k0 * k0 - k * k) * (k0 * k0 - k * k) := by
  unfold well multAdd; grind

/-- The well is even: it respects the flip. -/
theorem well_even (k0 k : R) : well k0 (-k) = well k0 k := by
  unfold well multAdd; grind

/-- Zero at both dips. -/
theorem well_dips (k0 : R) : well k0 k0 = 0 ∧ well k0 (-k0) = 0 := by
  unfold well multAdd; constructor <;> grind

/-- The hump at the symmetric point k = 0 has height k0⁴. -/
theorem well_hump (k0 : R) : well k0 0 = k0 * k0 * (k0 * k0) := by
  unfold well multAdd; grind

/-- Expansion at the dip: no linear term, curvature 8k0² (twice the h² coefficient 4k0²). -/
theorem well_at_dip (k0 h : R) :
    well k0 (k0 + h) = 4 * k0 * k0 * (h * h) + 4 * k0 * (h * h * h) + h * h * h * h := by
  unfold well multAdd; grind

end Ring

/-! ## 2. Potentials that respect the dial's rescaling -/

/-- A potential that scales as V(s·k) = sⁿ V(k) is a single power: V(k) = V(1)·kⁿ. -/
theorem scale_covariant_is_power (V : Rat → Rat) (n : Nat)
    (h : ∀ s k : Rat, V (s * k) = s ^ n * V k) (k : Rat) : V k = V 1 * k ^ n := by
  have := h k 1
  rw [Rat.mul_one] at this
  rw [this, Rat.mul_comm]

theorem well_one_one : well (1 : Rat) 1 = 0 := by unfold well multAdd; grind
theorem well_one_two : well (1 : Rat) (2 * 1) = 9 := by unfold well multAdd; grind

/-- The double well cannot respect the rescaling for any power n: a scale from outside sets k0. -/
theorem well_breaks_scale (n : Nat) : ¬ (∀ s k : Rat, well 1 (s * k) = s ^ n * well 1 k) := by
  intro h
  have h2 := h 2 1
  rw [well_one_one, well_one_two, Rat.mul_zero] at h2
  grind

/-! ## 3b. The usual Higgs normalisation -/

/-- (λ/4)((v + x)² − v²)² = ½(2λv²)x² + λv x³ + (λ/4)x⁴: the mass squared is 2λv². -/
theorem quartic_mass (lam v x : Rat) :
    lam / 4 * (((v + x) * (v + x) - v * v) * ((v + x) * (v + x) - v * v))
      = 1 / 2 * (2 * lam * (v * v)) * (x * x) + lam * v * (x * x * x) + lam / 4 * (x * x * x * x) := by
  grind

/-- Hump height = m²v²/8 with m² = 2λv². -/
theorem hump_height (lam v : Rat) :
    lam / 4 * ((0 - v * v) * (0 - v * v)) = (2 * lam * (v * v)) * (v * v) / 8 := by
  grind

/-! ## 4. The rolling dial: the exponential potential's fixed point -/

/-- With x² = λ²/6 and y² = 1 − λ²/6 (L = λ²): Ω_φ = x² + y² = 1 and x² − y² = L/3 − 1, so w = L/3 − 1. -/
theorem exp_fixed_point (L : Rat) :
    L / 6 + (1 - L / 6) = 1 ∧ L / 6 - (1 - L / 6) = L / 3 - 1 := by
  constructor <;> grind

/-- Conformal slope n = 1 (λ² = 2/3): w = −7/9 at the fixed point. -/
theorem w_conformal_one : (2 / 3 : Rat) / 3 - 1 = -7 / 9 := by grind

/-- Conformal slope n = 2 (λ² = 8/3): w = −1/9 at the fixed point. -/
theorem w_conformal_two : (8 / 3 : Rat) / 3 - 1 = -1 / 9 := by grind

/-- The fixed point exists (y² ≥ 0) only for λ² ≤ 6. -/
theorem exp_fixed_point_exists (L : Rat) (h : L ≤ 6) : 0 ≤ 1 - L / 6 := by
  grind

/-! ## 5. Two guesses ruled out by the measured masses (units of 10 MeV) -/

/-- Coleman–Weinberg (Gildener–Weinberg) Higgs: 6mW⁴ + 3mZ⁴ − 12mt⁴ < 0, so m_h² < 0,
    with mW = 80.37, mZ = 91.19, mt = 172.57 GeV.  Holds with a wide margin (factor 22). -/
theorem cw_mass_negative : 6 * 8037 ^ 4 + 3 * 9119 ^ 4 < 12 * 17257 ^ 4 := by decide

/-- Even with the top 10 GeV lighter and W, Z 10 GeV heavier the sign stays negative. -/
theorem cw_mass_negative_robust : 6 * 9037 ^ 4 + 3 * 10119 ^ 4 < 12 * 16257 ^ 4 := by decide

/-- Veltman: 4mt² − 2mW² − mZ² > 4·mH² (mH = 125.20 GeV), so the Veltman Higgs is over twice as heavy. -/
theorem veltman_too_heavy : 4 * 12520 ^ 2 + 2 * 8037 ^ 2 + 9119 ^ 2 < 4 * 17257 ^ 2 := by decide

/-! ## 6. The conformal reading's held-out test -/

/-- A conformal dial couples to matter with α0² = 1/6 (the same factor that gives slope √(2/3));
    then PPN γ = 1 − 2α0²/(1 + α0²) = 5/7, against Cassini's 1 ± 2.3·10⁻⁵. -/
theorem conformal_gamma : (1 : Rat) - 2 * (1 / 6) / (1 + 1 / 6) = 5 / 7 := by grind

/-- Slope and coupling are tied: for V ∝ kⁿ the slope is n·√(2/3), so λ² = (2/3)n² = 4n²·α0². -/
theorem slope_coupling (n : Rat) : (2 / 3) * (n * n) = 4 * (n * n) * (1 / 6) := by grind

end DialPotential
