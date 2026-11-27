/-
Fused Hyperoperation Algebra. Released under Apache 2.0.
-/
import FLTLeads.FreyMinimal
import Mathlib.Tactic.ComputeDegree
import Mathlib.Tactic.FieldSimp

/-!
# Split multiplicative reduction of the Frey model is our node law

Where the Frey model has multiplicative reduction and `a₄` reduces to `0` (that is, at the
primes dividing `AB`), the reduced curve is our nodal cubic `Y² + XY = X³ + a₂X²` with its node
at the origin. Its two tangent slopes are the roots of `T² + T - a₂`, and our port's
`node_chord_multAdd` shows that, once those slopes are in the residue field, the chord law of the
node is the circle operation `x ∘₁ y = x + y + xy` (called mult-add in the port's names): the
multiplicative group in shifted coordinates.

Mathlib's `HasSplitMultiplicativeReduction` asks that a certain quadratic built from the minimal
model split over the residue field. For the Frey model that quadratic is `c₄ · (T² + T - a₂)`,
so mathlib's notion and ours agree (`frey_hasSplit_iff`). Away from `2` this becomes the
classical test: split exactly when `B - A` is a square in the residue field
(`frey_hasSplit_iff_isSquare`); at a prime dividing `A` that is `B`, at one dividing `B` it is
`-A`.

This is known mathematics (Silverman VII.5, and the usual account of the Frey curve). The point
is that it now connects mathlib's reduction API to the circle-operation node law in Lean.
-/

open WeierstrassCurve IsDiscreteValuationRing IsDedekindDomain.HeightOneSpectrum Polynomial
open IsLocalRing (residue ResidueField)

namespace FreyMinimal

variable (R : Type*) [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
variable (K : Type*) [Field K] [Algebra R K] [IsFractionRing R K]

omit [IsDomain R] [IsDiscreteValuationRing R] in
/-- The integral model chosen by mathlib for `W⁄K` is `W` itself. -/
theorem integralModel_baseChange (W : WeierstrassCurve R) :
    integralModel R (W.baseChange K) = W := by
  have hinj := IsFractionRing.injective R K
  ext
  · exact hinj (by rw [integralModel_a₁_eq]; simp [baseChange])
  · exact hinj (by rw [integralModel_a₂_eq]; simp [baseChange])
  · exact hinj (by rw [integralModel_a₃_eq]; simp [baseChange])
  · exact hinj (by rw [integralModel_a₄_eq]; simp [baseChange])
  · exact hinj (by rw [integralModel_a₆_eq]; simp [baseChange])

/-- A monic quadratic `T² + T - a` splits exactly when it has a root. -/
theorem splits_quadratic_iff {k : Type*} [Field k] (a : k) :
    Splits (X ^ 2 + X - C a) ↔ ∃ t : k, t ^ 2 + t - a = 0 := by
  have hdeg : (X ^ 2 + X - C a : k[X]).natDegree = 2 := by compute_degree!
  constructor
  · intro h
    obtain ⟨t, ht⟩ := h.exists_eval_eq_zero (by
      rw [degree_eq_natDegree (by intro h0; simp [h0] at hdeg), hdeg]; decide)
    exact ⟨t, by simpa using ht⟩
  · rintro ⟨t, ht⟩
    exact Splits.of_natDegree_eq_two hdeg (by simpa using ht)

/-- Away from characteristic `2`, `T² + T - a` has a root exactly when `1 + 4a` is a square. -/
theorem exists_root_iff_isSquare {k : Type*} [Field k] (h2 : (2 : k) ≠ 0) (a : k) :
    (∃ t : k, t ^ 2 + t - a = 0) ↔ IsSquare (1 + 4 * a) := by
  constructor
  · rintro ⟨t, ht⟩
    exact ⟨2 * t + 1, by linear_combination -4 * ht⟩
  · rintro ⟨s, hs⟩
    refine ⟨(s - 1) / 2, ?_⟩
    field_simp
    linear_combination -hs

variable {R K}

/-- **Split multiplicative reduction of the Frey model is the split node.** Where the reduction
is multiplicative and `a₄` reduces to `0`, mathlib's split condition holds exactly when the
node's slope polynomial `T² + T - a₂` has a root in the residue field, which is when our port's
circle-operation chord law applies. -/
theorem frey_hasSplit_iff {a₂ a₄ : ℤ}
    (hm : HasMultiplicativeReduction R (((freyModel a₂ a₄).map (Int.castRingHom R)).baseChange K))
    (ha₄ : ((a₄ : ℤ) : R) ∈ IsLocalRing.maximalIdeal R) :
    HasSplitMultiplicativeReduction R (((freyModel a₂ a₄).map (Int.castRingHom R)).baseChange K)
      ↔ ∃ t : ResidueField R, t ^ 2 + t - (a₂ : ResidueField R) = 0 := by
  set W := (freyModel a₂ a₄).map (Int.castRingHom R) with hW
  have hWf : W = freyModel (a₂ : R) (a₄ : R) := by rw [hW, freyModel_map]; simp
  have hc : IsUnit W.c₄ := by
    have := hm.multiplicativeReduction
    rw [baseChange, map_c₄] at this
    exact (valuation_algebraMap_eq_one_iff R K _).1 this
  have hres0 : residue R (a₄ : R) = 0 := (IsLocalRing.residue_eq_zero_iff _).2 ha₄
  have hcres : residue R W.c₄ ≠ 0 := (hc.map (residue R)).ne_zero
  have key : (C W.c₄ * X ^ 2 + C (W.a₁ * W.c₄) * X -
      C (54 * W.b₆ - 3 * W.b₂ * W.b₄ + W.a₂ * W.c₄)).map (algebraMap R (ResidueField R)) =
      C (residue R W.c₄) * (X ^ 2 + X - C (a₂ : ResidueField R)) := by
    rw [IsLocalRing.ResidueField.algebraMap_eq]
    simp only [Polynomial.map_sub, Polynomial.map_add, Polynomial.map_mul, Polynomial.map_pow,
      Polynomial.map_C, Polynomial.map_X]
    rw [hWf] at hcres ⊢
    simp only [freyModel, b₂, b₄, b₆, map_mul, map_add, map_sub, map_ofNat, map_intCast,
      map_one, map_zero, map_pow] at hres0 ⊢
    have hpoly : ((a₄ : ℤ) : (ResidueField R)[X]) = 0 := by
      rw [← Polynomial.C_eq_intCast, hres0, C_0]
    rw [hpoly]
    ring
  have hsplit : Splits ((C W.c₄ * X ^ 2 + C (W.a₁ * W.c₄) * X -
      C (54 * W.b₆ - 3 * W.b₂ * W.b₄ + W.a₂ * W.c₄)).map (algebraMap R (ResidueField R))) ↔
      ∃ t : ResidueField R, t ^ 2 + t - (a₂ : ResidueField R) = 0 := by
    rw [key, splits_mul_iff_right (by simpa using hcres)
      (Splits.of_degree_le_zero degree_C_le), splits_quadratic_iff]
  constructor
  · intro h
    have := h.splitMultiplicativeReduction
    rw [integralModel_baseChange R K W] at this
    exact hsplit.1 this
  · intro h
    refine ⟨?_⟩
    show Splits _
    rw [integralModel_baseChange R K W]
    exact hsplit.2 h

/-- **Away from `2`, the Frey model splits exactly when `B - A` is a square mod the prime.** -/
theorem frey_hasSplit_iff_isSquare {A B a₂ a₄ : ℤ} (h2 : 4 * a₂ = B - A - 1)
    (h2R : (2 : ResidueField R) ≠ 0)
    (hm : HasMultiplicativeReduction R (((freyModel a₂ a₄).map (Int.castRingHom R)).baseChange K))
    (ha₄ : ((a₄ : ℤ) : R) ∈ IsLocalRing.maximalIdeal R) :
    HasSplitMultiplicativeReduction R (((freyModel a₂ a₄).map (Int.castRingHom R)).baseChange K)
      ↔ IsSquare ((B - A : ℤ) : ResidueField R) := by
  rw [frey_hasSplit_iff hm ha₄, exists_root_iff_isSquare h2R]
  have : (1 + 4 * a₂ : ℤ) = B - A := by linear_combination h2
  rw [← this]; push_cast; rfl

end FreyMinimal
