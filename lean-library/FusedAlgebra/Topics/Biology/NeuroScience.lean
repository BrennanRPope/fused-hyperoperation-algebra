/-
  Fused Hyperoperation Algebra: neuroscience (synaptic summation on passive dendrites).
  Core Lean 4.37 only (no mathlib), no `sorry`. Imports FusedAlgebra.Base only.
  Intended library path: FusedAlgebra/Topics/Biology/NeuroScience.lean

  1. One compartment: two conductance synapses with a shared reversal potential (units of
     the reversal, leak = 1) sum by the lens law (a + b − 2ab)/(1 − ab); in odds v/(1 − v)
     they simply add, and the reversal potential v = 1 absorbs.
  2. Any passive linear network: the soma response to two conductance synapses obeys the
     exact lens law V₁₂ = (a + b + k·ab)/(1 − q·ab) with
       k = −Z_ij(1/(E_j Z_si) + 1/(E_i Z_sj)),  q = Z_ij²/(E_i E_j Z_si Z_sj),
     proved here from the two-port equations in cleared form.
  3. Writing k = −(A + B), q = A·B, the law factorises: 1 − A·V₁₂ = (1 − Aa)(1 − Ab)/(1 − ABab),
     so T(v) = (1 − Av)/(1 − Bv) turns summation into multiplication. With q = 0 the lens law
     is the dial a + b + k·ab, the rule Hao et al. (2009) measured in CA1 pyramidal cells.
  4. Release from n sites, each with chance p, is wpow(−1).
  LIBRARY NOTE: first compiled in the shared library. 3 exact twins removed (star_def,
  release, weber); each is marked where it stood. No proof step or statement changed.
-/
import FusedAlgebra.Base

namespace NeuroScience
open CircleOp

set_option linter.unusedSectionVars false

open scoped CircleOp

section Ring
variable {R : Type} [Lean.Grind.CommRing R]

-- LIBRARY: star_def was an exact twin of MolecularBiology.star_def and is removed.

-- LIBRARY: release (Chance that at least one of n release sites fires: the "or" power) was an exact twin of Metabolism.limitation_pow and is removed.

-- LIBRARY: weber (Weber steps: n just-noticeable steps of relative size w multiply the stimulus by (1 + w)ⁿ) was an exact twin of MolecularBiology.mwc_R and is removed.

/-- In the soma-line case B = 1 the two constants sum to −1. -/
theorem soma_line (A : R) : -(A + 1) + A * 1 = -1 := by grind

/-- **The two-port theorem, cleared form.** Currents I₁ (synapse i alone), J₂ (synapse j alone)
    and I, J (both on) solve the network equations of a passive linear network with transfer
    impedances Z. Then the soma voltage V = Z_si I + Z_sj J obeys the lens law with the
    constants above, multiplied through by E_i E_j. -/
theorem twoport (Zsi Zsj Zii Zjj Zij gi gj Ei Ej I1 J2 I J : R)
    (h1 : I1 * (1 + gi * Zii) = gi * Ei)
    (h2 : J2 * (1 + gj * Zjj) = gj * Ej)
    (hI : I * (1 + gi * Zii) + gi * Zij * J = gi * Ei)
    (hJ : J * (1 + gj * Zjj) + gj * Zij * I = gj * Ej) :
    (Zsi * I + Zsj * J) * (Ei * Ej - Zij ^ 2 * I1 * J2)
      = Ei * Ej * (Zsi * I1 + Zsj * J2) - Zij * (Ei * Zsj + Ej * Zsi) * I1 * J2 := by
  grind

/-- Bookkeeping for the divided two-port theorem: clearing the lens numerator. -/
theorem clear_num (P Q x y EE r₁ r₂ : R) (t₁ : x * P * Q * EE = r₁) (t₂ : y * P * Q * EE = r₂) :
    (P + Q + -(x + y) * P * Q) * EE = EE * (P + Q) - (r₁ + r₂) := by
  rw [← t₁, ← t₂]; grind

/-- Bookkeeping for the divided two-port theorem: clearing the lens denominator. -/
theorem clear_den (P Q z EE r : R) (t : z * P * Q * EE = r) : (1 - z * P * Q) * EE = EE - r := by
  rw [← t]; grind

end Ring

section Field
variable {K : Type} [Lean.Grind.Field K]

/-- x / y · y = x for y ≠ 0. -/
theorem div_mul_cancel' (x y : K) (hy : y ≠ 0) : x / y * y = x := by
  rw [Lean.Grind.Field.div_eq_mul_inv, Lean.Grind.Semiring.mul_assoc,
    Lean.Grind.Field.inv_mul_cancel hy, Lean.Grind.Semiring.mul_one]

/-- The lens law: the dial a + b + k·ab, divided by 1 − q·ab. -/
def lens (k q a b : K) : K := (a + b + k * a * b) / (1 - q * a * b)

/-- With q = 0 the lens law is the dial: Hao et al.'s bilinear rule. -/
theorem lens_dial (k a b : K) : lens k 0 a b = circK k a b := by
  unfold lens circK; grind

/-- Steady-state voltage of one compartment (leak 1, reversal 1) with conductance g. -/
def vss (g : K) : K := g / (1 + g)

/-- **One compartment.** Two synapses sum by the lens law with k = −2, q = 1. -/
theorem one_compartment (g₁ g₂ : K) (h₁ : 1 + g₁ ≠ 0) (h₂ : 1 + g₂ ≠ 0)
    (h : 1 + (g₁ + g₂) ≠ 0) :
    vss (g₁ + g₂) = lens (-2) 1 (vss g₁) (vss g₂) := by
  unfold vss lens; grind

/-- In odds v/(1 − v) the conductance comes back: inputs add, dial 0. -/
theorem odds_vss (g : K) (h : 1 + g ≠ 0) : vss g / (1 - vss g) = g := by
  unfold vss; grind

/-- **Reversal absorbs.** A synapse holding the compartment at its reversal potential makes
    every other input irrelevant: v = 1 is the lens law's wall. -/
theorem reversal_absorbs (a : K) (h : 1 - a ≠ 0) : lens (-2) 1 a 1 = 1 := by
  unfold lens; grind

/-- **Factorisation.** With k = −(A + B), q = AB: 1 − A·V = (1 − Aa)(1 − Ab)/(1 − ABab). -/
theorem factor_A (A B a b : K) (h : 1 - A * B * a * b ≠ 0) :
    1 - A * lens (-(A + B)) (A * B) a b = (1 - A * a) * (1 - A * b) / (1 - A * B * a * b) := by
  unfold lens; grind

/-- The same factorisation with A and B exchanged. -/
theorem factor_B (A B a b : K) (h : 1 - A * B * a * b ≠ 0) :
    1 - B * lens (-(A + B)) (A * B) a b = (1 - B * a) * (1 - B * b) / (1 - A * B * a * b) := by
  unfold lens; grind

/-- **Summation becomes multiplication.** In T(v) = (1 − Av)/(1 − Bv) the lens law is a product. -/
theorem T_mul (A B a b : K) (h : 1 - A * B * a * b ≠ 0) (ha : 1 - B * a ≠ 0) (hb : 1 - B * b ≠ 0) :
    (1 - A * lens (-(A + B)) (A * B) a b) / (1 - B * lens (-(A + B)) (A * B) a b)
      = ((1 - A * a) / (1 - B * a)) * ((1 - A * b) / (1 - B * b)) := by
  rw [factor_A A B a b h, factor_B A B a b h]
  have i1 := Lean.Grind.Field.inv_mul_cancel ha
  have i2 := Lean.Grind.Field.inv_mul_cancel hb
  have hY : (1 - B * a) * (1 - B * b) ≠ 0 := by
    intro c
    have : (1 - B * a)⁻¹ * (1 - B * b)⁻¹ * ((1 - B * a) * (1 - B * b)) = 1 := by grind
    rw [c] at this; grind
  have iY := Lean.Grind.Field.inv_mul_cancel hY
  have iD := Lean.Grind.Field.inv_mul_cancel h
  have e : (1 - B * a) * (1 - B * b) / (1 - A * B * a * b) ≠ 0 := by
    intro c
    have := div_mul_cancel' ((1 - B * a) * (1 - B * b)) _ h
    rw [c] at this; grind
  have iE := Lean.Grind.Field.inv_mul_cancel e
  simp only [Lean.Grind.Field.div_eq_mul_inv] at *
  grind

/-- **Reversal at the soma absorbs, on any tree.** On the soma line (B = 1, reversal units)
    an input clamping the soma to its reversal fixes the response at 1. -/
theorem soma_reversal_absorbs (A a : K) (h : 1 - A * a ≠ 0) :
    lens (-(A + 1)) (A * 1) a 1 = 1 := by
  unfold lens; grind

/-- **The two-port theorem, divided form.** With nonzero reversals and transfer impedances,
    the soma voltage is exactly the lens law with the geometric constants k and q. -/
theorem twoport_lens (Zsi Zsj Zii Zjj Zij gi gj Ei Ej I1 J2 I J : K)
    (h1 : I1 * (1 + gi * Zii) = gi * Ei)
    (h2 : J2 * (1 + gj * Zjj) = gj * Ej)
    (hI : I * (1 + gi * Zii) + gi * Zij * J = gi * Ei)
    (hJ : J * (1 + gj * Zjj) + gj * Zij * I = gj * Ej)
    (hEi : Ei ≠ 0) (hEj : Ej ≠ 0) (hsi : Zsi ≠ 0) (hsj : Zsj ≠ 0)
    (hD : Ei * Ej - Zij ^ 2 * I1 * J2 ≠ 0) :
    Zsi * I + Zsj * J
      = lens (-(Zij / (Ej * Zsi) + Zij / (Ei * Zsj))) (Zij ^ 2 / (Ei * Ej * Zsi * Zsj))
          (Zsi * I1) (Zsj * J2) := by
  have c := twoport Zsi Zsj Zii Zjj Zij gi gj Ei Ej I1 J2 I J h1 h2 hI hJ
  have iEi := Lean.Grind.Field.inv_mul_cancel hEi
  have iEj := Lean.Grind.Field.inv_mul_cancel hEj
  have isi := Lean.Grind.Field.inv_mul_cancel hsi
  have isj := Lean.Grind.Field.inv_mul_cancel hsj
  -- the three pieces of the lens law, times Ei·Ej
  have t1 : Zij / (Ej * Zsi) * (Zsi * I1) * (Zsj * J2) * (Ei * Ej) = Zij * Ei * Zsj * I1 * J2 := by
    have : Zij / (Ej * Zsi) * (Zsi * I1) * (Zsj * J2) * (Ei * Ej)
        = Zij * Ei * Zsj * I1 * J2 * ((Ej⁻¹ * Ej) * (Zsi⁻¹ * Zsi)) := by
      rw [Lean.Grind.Field.div_eq_mul_inv, Lean.Grind.Field.inv_mul]; grind
    rw [this, iEj, isi]; grind
  have t2 : Zij / (Ei * Zsj) * (Zsi * I1) * (Zsj * J2) * (Ei * Ej) = Zij * Ej * Zsi * I1 * J2 := by
    have : Zij / (Ei * Zsj) * (Zsi * I1) * (Zsj * J2) * (Ei * Ej)
        = Zij * Ej * Zsi * I1 * J2 * ((Ei⁻¹ * Ei) * (Zsj⁻¹ * Zsj)) := by
      rw [Lean.Grind.Field.div_eq_mul_inv, Lean.Grind.Field.inv_mul]; grind
    rw [this, iEi, isj]; grind
  have t3 : Zij ^ 2 / (Ei * Ej * Zsi * Zsj) * (Zsi * I1) * (Zsj * J2) * (Ei * Ej) = Zij ^ 2 * I1 * J2 := by
    have : Zij ^ 2 / (Ei * Ej * Zsi * Zsj) * (Zsi * I1) * (Zsj * J2) * (Ei * Ej)
        = Zij ^ 2 * I1 * J2 * ((Ei⁻¹ * Ei) * (Ej⁻¹ * Ej) * (Zsi⁻¹ * Zsi) * (Zsj⁻¹ * Zsj)) := by
      rw [Lean.Grind.Field.div_eq_mul_inv, Lean.Grind.Field.inv_mul, Lean.Grind.Field.inv_mul,
        Lean.Grind.Field.inv_mul]; grind
    rw [this, iEi, iEj, isi, isj]; grind
  have eN : (Zsi * I1 + Zsj * J2 + -(Zij / (Ej * Zsi) + Zij / (Ei * Zsj)) * (Zsi * I1) * (Zsj * J2))
      * (Ei * Ej) = Ei * Ej * (Zsi * I1 + Zsj * J2) - Zij * (Ei * Zsj + Ej * Zsi) * I1 * J2 := by
    rw [clear_num _ _ _ _ _ _ _ t1 t2]; grind
  have eD : (1 - Zij ^ 2 / (Ei * Ej * Zsi * Zsj) * (Zsi * I1) * (Zsj * J2)) * (Ei * Ej)
      = Ei * Ej - Zij ^ 2 * I1 * J2 := by
    rw [clear_den _ _ _ _ _ t3]
  have hDn : 1 - Zij ^ 2 / (Ei * Ej * Zsi * Zsj) * (Zsi * I1) * (Zsj * J2) ≠ 0 := by
    intro z; rw [z] at eD; apply hD; grind
  have iDn := Lean.Grind.Field.inv_mul_cancel hDn
  unfold lens
  generalize hn : Zsi * I1 + Zsj * J2 + -(Zij / (Ej * Zsi) + Zij / (Ei * Zsj)) * (Zsi * I1) * (Zsj * J2) = N at *
  generalize hd : 1 - Zij ^ 2 / (Ei * Ej * Zsi * Zsj) * (Zsi * I1) * (Zsj * J2) = Dn at *
  rw [Lean.Grind.Field.div_eq_mul_inv]
  have : (Zsi * I + Zsj * J) * Dn * (Ei * Ej) = N * (Ei * Ej) := by grind
  have hEE : Ei * Ej ≠ 0 := by
    intro z
    have : Ei⁻¹ * Ej⁻¹ * (Ei * Ej) = 1 := by grind
    rw [z] at this; grind
  have iEE := Lean.Grind.Field.inv_mul_cancel hEE
  have h3 : (Zsi * I + Zsj * J) * Dn = N := by
    have e2 : (Zsi * I + Zsj * J) * Dn * ((Ei * Ej) * (Ei * Ej)⁻¹) = N * ((Ei * Ej) * (Ei * Ej)⁻¹) := by
      generalize (Ei * Ej)⁻¹ = w
      rw [← Lean.Grind.Semiring.mul_assoc, this, Lean.Grind.Semiring.mul_assoc]
    rw [Lean.Grind.Field.mul_inv_cancel hEE] at e2; grind
  grind

end Field

end NeuroScience
