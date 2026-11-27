/-
The first part (`isMinimal_of_valuation_c₄_eq_one`) is copied from the Lean FLT project,
`FLT/Mathlib/AlgebraicGeometry/EllipticCurve/Reduction.lean`:
Copyright (c) 2026 Kevin Buzzard. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Kevin Buzzard, William Coram, Claude.
The rest is ours (Fused Hyperoperation Algebra), also offered under Apache 2.0.
-/
import Mathlib.AlgebraicGeometry.EllipticCurve.Reduction
import Mathlib.Tactic.LinearCombination
import Mathlib.RingTheory.Coprime.Lemmas
import Mathlib.RingTheory.Int.Basic

/-!
# The Frey model is minimal, with good or multiplicative reduction, over every DVR

Let `A + B = C` with `gcd A B = 1`, and let `a₂, a₄ ∈ ℤ` with `4a₂ = B - A - 1` and
`16a₄ = -AB` (the Frey normalisation `A ≡ 3 mod 4`, `16 ∣ B` provides them). The model
`Y² + XY = X³ + a₂X² + a₄X` has `c₄ = A² + AB + B²` and `256Δ = (ABC)²`.

For every discrete valuation ring `R` with fraction field `K`, the model over `K` is minimal and
has good or multiplicative reduction (`frey_hasGood_or_hasMultiplicative`, `frey_isMinimal`).
This covers every prime, `2` included, and the `q`-adic integers in particular.
When `2` is a unit of `R`, the reduction is multiplicative exactly when `ABC` lies in the
maximal ideal (`frey_hasMultiplicative_iff`).

None of this is new mathematics: it is Silverman VII.1.1 (a model with `v(c₄) = 0` is minimal)
and the standard semistability of the Frey curve. What is new is having it checked in Lean
against mathlib's reduction API. The minimality criterion is not in mathlib yet; the Lean FLT
project has it and we copy it here with credit.
-/

open WeierstrassCurve IsDiscreteValuationRing IsDedekindDomain.HeightOneSpectrum

namespace FreyMinimal

/-! ### From the Lean FLT project -/

section FLT

variable (R : Type*) [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
variable {K : Type*} [Field K] [Algebra R K] [IsFractionRing R K]

/-- A Weierstrass equation with integral coefficients and `v(c₄) = 0` is minimal.
Copied from the Lean FLT project (Buzzard, Coram, Claude). -/
theorem isMinimal_of_valuation_c₄_eq_one (W : WeierstrassCurve K) [IsIntegral R W]
    (hc₄ : valuation K (maximalIdeal R) W.c₄ = 1) : IsMinimal R W := by
  refine ⟨⟨by simpa using ‹IsIntegral R W›, ?_⟩⟩
  intro C hC _
  simp only [one_smul, ← Subtype.coe_le_coe, valuation_Δ_aux_eq_of_isIntegral R (C • W),
    valuation_Δ_aux_eq_of_isIntegral R W]
  have hint : valuation K (maximalIdeal R) (C • W).c₄ ≤ 1 := by
    simpa [← integralModel_c₄_eq R (C • W)] using valuation_le_one _ _
  rw [variableChange_c₄, map_mul, map_pow, hc₄, mul_one] at hint
  simpa [variableChange_Δ, map_mul, map_pow] using mul_le_of_le_one_left'
    (pow_le_one' ((pow_le_one_iff (by norm_num)).mp hint) 12)

end FLT

/-! ### Reduction of an integral model, read off from units -/

section Units

variable (R : Type*) [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
variable (K : Type*) [Field K] [Algebra R K] [IsFractionRing R K]

/-- The valuation of an element of `R` is `1` exactly when it is a unit. -/
theorem valuation_algebraMap_eq_one_iff (r : R) :
    valuation K (maximalIdeal R) (algebraMap R K r) = 1 ↔ IsUnit r := by
  have h1 := valuation_le_one (K := K) (maximalIdeal R) r
  have h2 := valuation_lt_one_iff_mem (K := K) (maximalIdeal R) r
  have h3 : r ∈ (maximalIdeal R).asIdeal ↔ ¬IsUnit r := by
    change r ∈ IsLocalRing.maximalIdeal R ↔ _
    rw [IsLocalRing.mem_maximalIdeal, mem_nonunits_iff]
  constructor
  · intro h
    by_contra hu
    exact (h2.2 (h3.2 hu)).ne h
  · intro hu
    exact le_antisymm h1 (not_lt.1 fun hlt ↦ h3.1 (h2.1 hlt) hu)

/-- The valuation of an element of `R` is below `1` exactly when it is not a unit. -/
theorem valuation_algebraMap_lt_one_iff (r : R) :
    valuation K (maximalIdeal R) (algebraMap R K r) < 1 ↔ ¬IsUnit r := by
  rw [← valuation_algebraMap_eq_one_iff R K,
    (valuation_le_one (K := K) (maximalIdeal R) r).lt_iff_ne]

instance isIntegral_baseChange (W : WeierstrassCurve R) : IsIntegral R (W.baseChange K) :=
  ⟨⟨W, rfl⟩⟩

/-- An integral model whose discriminant is a unit has good reduction (and is minimal). -/
theorem hasGoodReduction_of_isUnit_Δ (W : WeierstrassCurve R) (h : IsUnit W.Δ) :
    HasGoodReduction R (W.baseChange K) := by
  have hv : valuation K (maximalIdeal R) (W.baseChange K).Δ = 1 := by
    rw [baseChange, map_Δ]; exact (valuation_algebraMap_eq_one_iff R K W.Δ).2 h
  have hmin : IsMinimal R (W.baseChange K) := by
    refine ⟨⟨(by simpa using isIntegral_baseChange R K W), fun C _ _ ↦ ?_⟩⟩
    dsimp only
    rw [← Subtype.coe_le_coe, one_smul, valuation_Δ_aux_eq_of_isIntegral R (W.baseChange K), hv]
    exact (valuation_Δ_aux R (C • W.baseChange K)).2
  exact ⟨hv⟩

/-- An integral model whose `c₄` is a unit and whose discriminant is not has multiplicative
reduction (and is minimal). -/
theorem hasMultiplicativeReduction_of_isUnit_c₄ (W : WeierstrassCurve R) (hc : IsUnit W.c₄)
    (hΔ : ¬IsUnit W.Δ) : HasMultiplicativeReduction R (W.baseChange K) := by
  have hv : valuation K (maximalIdeal R) (W.baseChange K).c₄ = 1 := by
    rw [baseChange, map_c₄]; exact (valuation_algebraMap_eq_one_iff R K W.c₄).2 hc
  have := isMinimal_of_valuation_c₄_eq_one R _ hv
  refine ⟨?_, hv⟩
  rw [baseChange, map_Δ]; exact (valuation_algebraMap_lt_one_iff R K W.Δ).2 hΔ

/-- If `Δ` and `c₄` of an integral model are coprime, the model has good or multiplicative
reduction. -/
theorem hasGood_or_hasMultiplicative_of_isCoprime (W : WeierstrassCurve R)
    (h : IsCoprime W.Δ W.c₄) :
    HasGoodReduction R (W.baseChange K) ∨ HasMultiplicativeReduction R (W.baseChange K) := by
  by_cases hΔ : IsUnit W.Δ
  · exact .inl (hasGoodReduction_of_isUnit_Δ R K W hΔ)
  · obtain ⟨u, v, huv⟩ := h
    have hc : IsUnit W.c₄ := by
      rcases IsLocalRing.isUnit_or_isUnit_of_isUnit_add (huv ▸ isUnit_one) with hu | hu
      · exact absurd (isUnit_of_mul_isUnit_right hu) hΔ
      · exact isUnit_of_mul_isUnit_right hu
    exact .inr (hasMultiplicativeReduction_of_isUnit_c₄ R K W hc hΔ)

end Units

/-! ### The Frey model

`freyModel`, `frey_c4`, `frey_disc`, `frey_c4_not_dvd` and `frey_semistable` are the statements
of our port's `FreyPieces`, repeated here so that this folder builds on mathlib alone. -/

section Frey

/-- The model `Y² + XY = X³ + a₂X² + a₄X` (`a₁ = 1`, `a₃ = a₆ = 0`). -/
def freyModel {R : Type*} [CommRing R] (a₂ a₄ : R) : WeierstrassCurve R := ⟨1, a₂, 0, a₄, 0⟩

theorem freyModel_map {R S : Type*} [CommRing R] [CommRing S] (f : R →+* S) (a₂ a₄ : R) :
    (freyModel a₂ a₄).map f = freyModel (f a₂) (f a₄) := by
  ext <;> simp [freyModel]

theorem frey_c4 {R : Type*} [CommRing R] {A B a₂ a₄ : R} (h2 : 4 * a₂ = B - A - 1)
    (h4 : 16 * a₄ = -(A * B)) : (freyModel a₂ a₄).c₄ = A ^ 2 + A * B + B ^ 2 := by
  simp only [freyModel, WeierstrassCurve.c₄, WeierstrassCurve.b₂, WeierstrassCurve.b₄]
  linear_combination (4 * a₂ + B - A + 1) * h2 - 3 * h4

theorem frey_disc {R : Type*} [CommRing R] {A B C a₂ a₄ : R} (hC : A + B = C)
    (h2 : 4 * a₂ = B - A - 1) (h4 : 16 * a₄ = -(A * B)) :
    256 * (freyModel a₂ a₄).Δ = (A * B * C) ^ 2 := by
  subst hC
  have hb : 1 + 4 * a₂ = B - A := by linear_combination h2
  calc 256 * (freyModel a₂ a₄).Δ = (16 * a₄) ^ 2 * ((1 + 4 * a₂) ^ 2 - 4 * (16 * a₄)) := by
        simp only [freyModel, WeierstrassCurve.Δ, WeierstrassCurve.b₂, WeierstrassCurve.b₄,
          WeierstrassCurve.b₆, WeierstrassCurve.b₈]
        ring
    _ = (A * B * (A + B)) ^ 2 := by rw [hb, h4]; ring

theorem prime_dvd_gcd_false {q : ℕ} (hq : q.Prime) {A B : ℤ} (hcop : Int.gcd A B = 1)
    (hA : (q : ℤ) ∣ A) (hB : (q : ℤ) ∣ B) : False := by
  have h := Int.dvd_coe_gcd hA hB
  rw [hcop] at h
  exact hq.one_lt.ne' (Nat.dvd_one.1 (Int.natCast_dvd_natCast.1 h))

theorem frey_c4_not_dvd {A B C : ℤ} (hC : A + B = C) (hcop : Int.gcd A B = 1) {q : ℕ}
    (hq : q.Prime) (hd : (q : ℤ) ∣ A * B * C) : ¬(q : ℤ) ∣ A ^ 2 + A * B + B ^ 2 := by
  have hq' : Prime (q : ℤ) := Nat.prime_iff_prime_int.1 hq
  subst hC
  intro hc
  obtain (hA | hB) | hAB := (hq'.dvd_or_dvd hd).imp_left hq'.dvd_or_dvd
  · refine prime_dvd_gcd_false hq hcop hA (hq'.dvd_of_dvd_pow (n := 2) ?_)
    have e : B ^ 2 = A ^ 2 + A * B + B ^ 2 - A * (A + B) := by ring
    rw [e]; exact dvd_sub hc (hA.mul_right _)
  · refine prime_dvd_gcd_false hq hcop (hq'.dvd_of_dvd_pow (n := 2) ?_) hB
    have e : A ^ 2 = A ^ 2 + A * B + B ^ 2 - B * (A + B) := by ring
    rw [e]; exact dvd_sub hc (hB.mul_right _)
  · have hA : (q : ℤ) ∣ A := by
      refine hq'.dvd_of_dvd_pow (n := 2) ?_
      have e : A ^ 2 = A ^ 2 + A * B + B ^ 2 - B * (A + B) := by ring
      rw [e]; exact dvd_sub hc (hAB.mul_left _)
    exact prime_dvd_gcd_false hq hcop hA ((dvd_add_right hA).1 hAB)

theorem frey_semistable {A B C a₂ a₄ : ℤ} (hC : A + B = C) (hcop : Int.gcd A B = 1)
    (h2 : 4 * a₂ = B - A - 1) (h4 : 16 * a₄ = -(A * B)) {q : ℕ} (hq : q.Prime)
    (hd : (q : ℤ) ∣ (freyModel a₂ a₄).Δ) : ¬(q : ℤ) ∣ (freyModel a₂ a₄).c₄ := by
  rw [frey_c4 h2 h4]
  refine frey_c4_not_dvd hC hcop hq ((Nat.prime_iff_prime_int.1 hq).dvd_of_dvd_pow (n := 2) ?_)
  rw [← frey_disc hC h2 h4]
  exact hd.mul_left 256

/-- Over `ℤ`, `Δ` and `c₄` of the Frey model are coprime: semistability at every prime at once. -/
theorem frey_isCoprime {A B C a₂ a₄ : ℤ} (hC : A + B = C) (hcop : Int.gcd A B = 1)
    (h2 : 4 * a₂ = B - A - 1) (h4 : 16 * a₄ = -(A * B)) :
    IsCoprime (freyModel a₂ a₄).Δ (freyModel a₂ a₄).c₄ := by
  rw [Int.isCoprime_iff_gcd_eq_one]
  by_contra hg
  obtain ⟨q, hq, hqg⟩ := Nat.exists_prime_and_dvd hg
  have hqg' : (q : ℤ) ∣ (Int.gcd (freyModel a₂ a₄).Δ (freyModel a₂ a₄).c₄ : ℤ) :=
    Int.natCast_dvd_natCast.2 hqg
  exact frey_semistable hC hcop h2 h4 hq (hqg'.trans (Int.gcd_dvd_left _ _))
    (hqg'.trans (Int.gcd_dvd_right _ _))

variable {A B C a₂ a₄ : ℤ} (hC : A + B = C) (hcop : Int.gcd A B = 1)
  (h2 : 4 * a₂ = B - A - 1) (h4 : 16 * a₄ = -(A * B))
variable (R : Type*) [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
variable (K : Type*) [Field K] [Algebra R K] [IsFractionRing R K]

include hC hcop h2 h4

/-- **The Frey model has good or multiplicative reduction over every DVR.** -/
theorem frey_hasGood_or_hasMultiplicative :
    HasGoodReduction R (((freyModel a₂ a₄).map (Int.castRingHom R)).baseChange K) ∨
      HasMultiplicativeReduction R
        (((freyModel a₂ a₄).map (Int.castRingHom R)).baseChange K) := by
  refine hasGood_or_hasMultiplicative_of_isCoprime R K _ ?_
  rw [map_Δ, map_c₄]
  exact (frey_isCoprime hC hcop h2 h4).map _

/-- **The Frey model is minimal over every DVR**, the prime `2` included. -/
theorem frey_isMinimal :
    IsMinimal R (((freyModel a₂ a₄).map (Int.castRingHom R)).baseChange K) := by
  rcases frey_hasGood_or_hasMultiplicative hC hcop h2 h4 R K with h | h
  · exact h.toIsMinimal
  · exact h.toIsMinimal

omit hcop [IsDomain R] [IsDiscreteValuationRing R] in
/-- Away from `2`, `Δ` of the Frey model is a non-unit exactly when `ABC` is. -/
theorem frey_isUnit_Δ_iff (h2R : IsUnit (2 : R)) :
    IsUnit ((freyModel a₂ a₄).map (Int.castRingHom R)).Δ ↔ IsUnit ((A * B * C : ℤ) : R) := by
  have key := congrArg (Int.castRingHom R) (frey_disc hC h2 h4)
  simp only [map_mul, map_pow, map_ofNat, eq_intCast] at key
  have h256 : IsUnit (256 : R) := by
    have := h2R.pow 8; norm_num at this; exact this
  rw [map_Δ, eq_intCast, ← isUnit_pow_iff (n := 2) two_ne_zero (a := ((A * B * C : ℤ) : R)),
    Int.cast_mul, Int.cast_mul, ← key, IsUnit.mul_iff]
  simp [h256]

/-- Away from `2`, the reduction is multiplicative exactly when `ABC` lies in the maximal ideal,
and good otherwise. -/
theorem frey_hasMultiplicative_iff (h2R : IsUnit (2 : R)) :
    HasMultiplicativeReduction R (((freyModel a₂ a₄).map (Int.castRingHom R)).baseChange K) ↔
      ((A * B * C : ℤ) : R) ∈ IsLocalRing.maximalIdeal R := by
  rw [IsLocalRing.mem_maximalIdeal, mem_nonunits_iff,
    ← frey_isUnit_Δ_iff hC h2 h4 R h2R]
  constructor
  · intro hm hu
    have hv := hm.badReduction
    rw [baseChange, map_Δ] at hv
    exact (valuation_algebraMap_lt_one_iff R K _).1 hv hu
  · intro hu
    rcases frey_hasGood_or_hasMultiplicative hC hcop h2 h4 R K with hg | hm
    · have hv := hg.goodReduction
      rw [baseChange, map_Δ] at hv
      exact absurd ((valuation_algebraMap_eq_one_iff R K _).1 hv) hu
    · exact hm

end Frey

end FreyMinimal
