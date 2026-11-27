/-
  Reflection positivity for the weighted-circle chain (topics/why-i/reflection-positivity).
  Core Lean 4 only (no mathlib), no `sorry`. Library path: FusedAlgebra/Topics/ReflectionPositivity.lean.

  The exact facts behind the pictures:
  1. A chain of n equal slices has multiplier (1 + k·a)^n (wpow_lift, restated).
  2. Site reflection, two states, three times: the brute-force path sum <θF, F> equals a sum
     of squares when the slice is symmetric, so it is never negative.
  3. Link reflection, two states: <θF, F> = fᵀ T f is never negative when the 2×2 slice is
     positive semidefinite, and is negative for some F when a diagonal entry is negative.
  4. The explicit slice 1 + k·ε·E is nonnegative exactly on a window of dial values.
  5. The rewind slice (1 − k·ε·E)⁻¹ is positive for every step size ε exactly when k·E ≤ 0
     (the spectrum condition).
  6. Off the real axis the slice is not Hermitian, and the simplest pairing is not real;
     at a purely imaginary dial (Feynman's) it is real and strictly negative.
  7. The Wick turn: the real 2×2 form of z·h is linear in h, so an energy shift adds the same
     rotation for every h, and at z = −i the generator is antisymmetric and squares to −h².
  8. The positive-multiplier side of the wall −1/k is closed under ∘ₖ and its inverse, and the
     wall absorbs.

  First compiled in the shared library on Nov 27, as FusedAlgebra.Topics.ReflectionPositivity
  (namespace CircleOp.ReflectionPositivity, the library's convention; the topic copy's MultAdd,
  wAdd and wAdd_lift renamed to the library's CircleOp, circK and circK_lift; no proof changed).
-/
import FusedAlgebra.Base

namespace CircleOp
namespace ReflectionPositivity

open Lean.Grind

theorem sq_nn (x : Rat) : 0 ≤ x * x := by
  have h := OrderedRing.sq_nonneg (a := x)
  have e : x ^ 2 = x * x := by grind
  rw [e] at h; exact h

/-! ## 1. The chain's multiplier -/

/-- n equal slices a at dial k: 1 + k·(a ∘ₖ … ∘ₖ a) = (1 + k·a)^n. -/
theorem slices_lift (k a : Rat) (n : Nat) : 1 + k * wpow k a n = (1 + k * a) ^ n :=
  wpow_lift k a n

/-- Two different slices: the multiplier of the pair is the product. -/
theorem two_slices_lift (k a b : Rat) : 1 + k * circK k a b = (1 + k * a) * (1 + k * b) :=
  circK_lift k a b

/-! ## 2. Site reflection: a sum of squares -/

/-- Sum over the two states. -/
def sumB (f : Bool → Rat) : Rat := f false + f true

/-- <θF, F> for paths x0 → x1 → x2, mirror through x1, F looks at (x1, x2), free ends.
    θF looks at (x1, x0).  Weight of the path: T x1 x0 · T x2 x1. -/
def sitePair (T : Bool → Bool → Rat) (F : Bool → Bool → Rat) : Rat :=
  sumB fun x0 => sumB fun x1 => sumB fun x2 => F x1 x0 * F x1 x2 * (T x1 x0 * T x2 x1)

theorem sitePair_square (T F : Bool → Bool → Rat) (hT : ∀ a b, T a b = T b a) :
    sitePair T F = sumB fun x1 => (sumB fun x2 => F x1 x2 * T x2 x1) *
                                   (sumB fun x2 => F x1 x2 * T x2 x1) := by
  have h1 := hT false true
  simp only [sitePair, sumB]
  rw [h1]
  grind

/-- Site reflection positivity: a symmetric slice never gives a negative pairing. -/
theorem site_rp (T F : Bool → Bool → Rat) (hT : ∀ a b, T a b = T b a) : 0 ≤ sitePair T F := by
  rw [sitePair_square T F hT]
  simp only [sumB]
  have := sq_nn (F false false * T false false + F false true * T true false)
  have := sq_nn (F true false * T false true + F true true * T true true)
  grind

/-! ## 3. Link reflection: a positive slice -/

/-- <θF, F> across one link, two states: f₁²·t₁₁ + 2 f₁ f₂ t₁₂ + f₂²·t₂₂. -/
def linkPair (t11 t12 t22 f1 f2 : Rat) : Rat :=
  f1 * t11 * f1 + f1 * t12 * f2 + f2 * t12 * f1 + f2 * t22 * f2

theorem link_rp (t11 t12 t22 f1 f2 : Rat) (h11 : 0 ≤ t11) (h22 : 0 ≤ t22)
    (hdet : t12 * t12 ≤ t11 * t22) : 0 ≤ linkPair t11 t12 t22 f1 f2 := by
  unfold linkPair
  by_cases h0 : t11 = 0
  · subst h0
    have hs := sq_nn t12
    have ht : t12 * t12 = 0 := by grind
    have : t12 = 0 := by
      rcases Rat.mul_eq_zero.mp ht with h | h <;> exact h
    subst this
    have := Rat.mul_nonneg (Rat.mul_nonneg h22 (sq_nn f2)) (show (0:Rat) ≤ 1 by decide)
    grind
  · have hpos : 0 < t11 := by grind
    have key : t11 * (f1 * t11 * f1 + f1 * t12 * f2 + f2 * t12 * f1 + f2 * t22 * f2)
        = (t11 * f1 + t12 * f2) * (t11 * f1 + t12 * f2) + (t11 * t22 - t12 * t12) * (f2 * f2) := by
      grind
    have hA := sq_nn (t11 * f1 + t12 * f2)
    have hB := Rat.mul_nonneg (show 0 ≤ t11 * t22 - t12 * t12 by grind) (sq_nn f2)
    have hprod : 0 ≤ t11 * (f1 * t11 * f1 + f1 * t12 * f2 + f2 * t12 * f1 + f2 * t22 * f2) := by
      rw [key]; grind
    by_cases hneg : 0 ≤ f1 * t11 * f1 + f1 * t12 * f2 + f2 * t12 * f1 + f2 * t22 * f2
    · exact hneg
    · have hlt : f1 * t11 * f1 + f1 * t12 * f2 + f2 * t12 * f1 + f2 * t22 * f2 < 0 := by grind
      have := OrderedRing.mul_neg_of_pos_of_neg hpos hlt
      grind

/-- A negative diagonal entry breaks link positivity: F = indicator of that state. -/
theorem link_fails_negative_diag (t11 t12 t22 : Rat) (h : t11 < 0) :
    linkPair t11 t12 t22 1 0 < 0 := by
  unfold linkPair; grind

/-! ## 4. The explicit slice's window -/

/-- On an energy level E the explicit slice is 1 + k·ε·E; it is nonnegative exactly when
    k·ε·E ≥ −1 (the identity's side of the wall). -/
theorem slice_nonneg_iff (k e E : Rat) : 0 ≤ 1 + k * e * E ↔ -1 ≤ k * e * E := by
  constructor <;> intro h <;> grind

/-! ## 5. Rewind slices and the spectrum condition -/

/-- If k·E ≤ 0, the rewind slice's denominator 1 − k·ε·E is at least 1 for every ε ≥ 0. -/
theorem rewind_denominator_ge_one (k e E : Rat) (he : 0 ≤ e) (hkE : k * E ≤ 0) :
    1 ≤ 1 - k * e * E := by
  have := OrderedRing.mul_nonpos_of_nonneg_of_nonpos he hkE
  grind

/-- The rewind slice (1 − k·ε·E)⁻¹ is then positive for every step size. -/
theorem rewind_pos (k e E : Rat) (he : 0 ≤ e) (hkE : k * E ≤ 0) :
    0 < (1 - k * e * E)⁻¹ := by
  have h := rewind_denominator_ge_one k e E he hkE
  have hp : 0 < 1 - k * e * E := by grind
  exact Rat.inv_pos.mpr hp

/-- Conversely, if k·E > 0 the step ε = 2/(k·E) puts the denominator at −1:
    the rewind slice turns negative and link positivity fails. -/
theorem rewind_fails_if_kE_pos (k E : Rat) (hkE : 0 < k * E) :
    1 - k * (2 / (k * E)) * E = -1 := by
  have : k * E ≠ 0 := by grind
  grind

/-- So: positive rewind slices at every step size ⟺ k·E ≤ 0. -/
theorem rewind_all_steps_iff (k E : Rat) :
    (∀ e : Rat, 0 ≤ e → 0 < 1 - k * e * E) ↔ k * E ≤ 0 := by
  constructor
  · intro h
    by_cases hc : k * E ≤ 0
    · exact hc
    · have hpos : 0 < k * E := by grind
      have hinv : 0 < (k * E)⁻¹ := Rat.inv_pos.mpr hpos
      have hdiv : 2 / (k * E) = 2 * (k * E)⁻¹ := by grind
      have h2 := h (2 / (k * E)) (by
        rw [hdiv]; have := Rat.mul_pos (show (0:Rat) < 2 by decide) hinv; grind)
      rw [rewind_fails_if_kE_pos k E hpos] at h2
      grind
  · intro hkE e he
    have := rewind_denominator_ge_one k e E he hkE
    grind

/-! ## 6. Off the real axis -/

/-- Complex numbers as pairs (re, im). -/
def cmul (a b : Rat × Rat) : Rat × Rat := (a.1 * b.1 - a.2 * b.2, a.1 * b.2 + a.2 * b.1)

/-- The off-diagonal slice entry at dial k = kr + i·ki, with c = ε·h₁₂, is (kr·c, ki·c), and for
    a symmetric H the two off-diagonal entries are equal. Hermitian needs conj = itself:
    ki·c = −(ki·c). With c ≠ 0 that forces ki = 0. -/
theorem hermitian_iff_real_dial (ki c : Rat) (hc : c ≠ 0) : ki * c = -(ki * c) ↔ ki = 0 := by
  constructor
  · intro h
    have : ki * c = 0 := by grind
    rcases Rat.mul_eq_zero.mp this with h1 | h1
    · exact h1
    · exact absurd h1 hc
  · intro h; subst h; grind

/-- The simplest pairing (F = indicator of x1 = 0, x2 = 1, three times) is T₁₂·T₂₁ = (k·c)².
    Its imaginary part is 2·kr·ki·c². -/
theorem pairing_im (kr ki c : Rat) :
    (cmul (kr * c, ki * c) (kr * c, ki * c)).2 = 2 * kr * ki * (c * c) := by
  unfold cmul; grind

theorem pairing_not_real (kr ki c : Rat) (hr : kr ≠ 0) (hi : ki ≠ 0) (hc : c ≠ 0) :
    (cmul (kr * c, ki * c) (kr * c, ki * c)).2 ≠ 0 := by
  rw [pairing_im]
  intro h
  have h1 : kr * ki ≠ 0 := by
    intro h2; rcases Rat.mul_eq_zero.mp h2 with h3 | h3
    · exact hr h3
    · exact hi h3
  have h2 : c * c ≠ 0 := by
    intro h2; rcases Rat.mul_eq_zero.mp h2 with h3 | h3 <;> exact hc h3
  have h3 : kr * ki * (c * c) ≠ 0 := by
    intro h4; rcases Rat.mul_eq_zero.mp h4 with h5 | h5
    · exact h1 h5
    · exact h2 h5
  grind

/-- Feynman's dial (purely imaginary): the pairing is real and strictly negative. -/
theorem feynman_pairing_negative (ki c : Rat) (hi : ki ≠ 0) (hc : c ≠ 0) :
    (cmul (0, ki * c) (0, ki * c)).2 = 0 ∧ (cmul (0, ki * c) (0, ki * c)).1 < 0 := by
  have hk : ki * c ≠ 0 := by
    intro h; rcases Rat.mul_eq_zero.mp h with h1 | h1
    · exact hi h1
    · exact hc h1
  have hs := sq_nn (ki * c)
  have hne : ki * c * (ki * c) ≠ 0 := by
    intro h; rcases Rat.mul_eq_zero.mp h with h1 | h1 <;> exact hk h1
  unfold cmul
  constructor
  · grind
  · grind

/-! ## 7. The Wick turn: one rotation for every energy -/

/-- Real 2×2 form of the complex number z·h, entries (m11, m12, m21, m22). -/
def realForm (z : Rat × Rat) (h : Rat) : Rat × Rat × Rat × Rat :=
  (z.1 * h, -(z.2 * h), z.2 * h, z.1 * h)

def add4 (a b : Rat × Rat × Rat × Rat) : Rat × Rat × Rat × Rat :=
  (a.1 + b.1, a.2.1 + b.2.1, a.2.2.1 + b.2.2.1, a.2.2.2 + b.2.2.2)

def smul4 (c : Rat) (a : Rat × Rat × Rat × Rat) : Rat × Rat × Rat × Rat :=
  (c * a.1, c * a.2.1, c * a.2.2.1, c * a.2.2.2)

def mul4 (a b : Rat × Rat × Rat × Rat) : Rat × Rat × Rat × Rat :=
  (a.1 * b.1 + a.2.1 * b.2.2.1, a.1 * b.2.1 + a.2.1 * b.2.2.2,
   a.2.2.1 * b.1 + a.2.2.2 * b.2.2.1, a.2.2.1 * b.2.1 + a.2.2.2 * b.2.2.2)

/-- (A4) after the turn: shifting the energy by c adds c times the same rotation, whatever h is. -/
theorem wick_shift (z : Rat × Rat) (h c : Rat) :
    realForm z (h + c) = add4 (realForm z h) (smul4 c (realForm z 1)) := by
  simp only [realForm, add4, smul4, Prod.mk.injEq]
  refine ⟨?_, ?_, ?_, ?_⟩ <;> grind

/-- (A1) at z = −i: the generator is antisymmetric. -/
theorem wick_antisymm (h : Rat) :
    (realForm (0, -1) h).1 = 0 ∧ (realForm (0, -1) h).2.2.2 = 0 ∧
    (realForm (0, -1) h).2.1 = -(realForm (0, -1) h).2.2.1 := by
  unfold realForm; grind

/-- (A3) at z = −i: the generator squares to −h². -/
theorem wick_square (h : Rat) :
    mul4 (realForm (0, -1) h) (realForm (0, -1) h) = (-(h * h), 0, 0, -(h * h)) := by
  simp only [mul4, realForm, Prod.mk.injEq]
  refine ⟨?_, ?_, ?_, ?_⟩ <;> grind

/-- The rotation an energy shift adds is the same J = realForm (0, −1) 1 for every h. -/
theorem shift_rotation_shared (h₁ h₂ c : Rat) :
    add4 (realForm (0, -1) (h₁ + c)) (smul4 (-1) (realForm (0, -1) h₁)) =
    add4 (realForm (0, -1) (h₂ + c)) (smul4 (-1) (realForm (0, -1) h₂)) := by
  simp only [add4, smul4, realForm, Prod.mk.injEq]
  refine ⟨?_, ?_, ?_, ?_⟩ <;> grind

/-! ## 8. The identity's side of the wall -/

theorem positive_side_closed (k a b : Rat) (ha : 0 < 1 + k * a) (hb : 0 < 1 + k * b) :
    0 < 1 + k * circK k a b := by
  rw [circK_lift]; exact Rat.mul_pos ha hb

/-- The ∘ₖ inverse −a/(1 + k·a) stays on the positive side: its multiplier is 1/(1 + k·a). -/
theorem positive_side_inverse (k a : Rat) (ha : 0 < 1 + k * a) :
    1 + k * (-a / (1 + k * a)) = (1 + k * a)⁻¹ ∧ 0 < 1 + k * (-a / (1 + k * a)) := by
  have hne : 1 + k * a ≠ 0 := by grind
  have e : 1 + k * (-a / (1 + k * a)) = (1 + k * a)⁻¹ := by grind
  exact ⟨e, by rw [e]; exact Rat.inv_pos.mpr ha⟩

theorem wall_absorbs (k b : Rat) (hk : k ≠ 0) : circK k (-1 / k) b = -1 / k := by
  unfold circK; grind

end ReflectionPositivity
end CircleOp
