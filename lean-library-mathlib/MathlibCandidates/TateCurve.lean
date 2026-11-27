module

public import Mathlib.AlgebraicGeometry.EllipticCurve.Reduction
public import Mathlib.NumberTheory.ArithmeticFunction.Misc
public import Mathlib.RingTheory.LaurentSeries
public import Mathlib.RingTheory.PowerSeries.Inverse
public import Mathlib.Algebra.Polynomial.Splits

/-!
# The Tate curve

DRAFT, written as a candidate for Mathlib.

Over any commutative ring `R`, the **Tate curve** is the Weierstrass curve over `R⟦q⟧`
$$ E_q : Y^2 + XY = X^3 + a_4(q) X + a_6(q), $$
$$ a_4 = -5 s_3, \qquad a_6 = -\frac{5 s_3 + 7 s_5}{12}, \qquad
  s_k = \sum_{n \ge 1} \sigma_k(n) q^n. $$
The division by `12` is exact coefficient by coefficient, since `12 ∣ 5d³ + 7d⁵` for every `d`,
so `a₆` has coefficients in the image of `ℤ` and the curve is defined over every ring.

## Main definitions

* `TateCurve.s R k`: the power series `∑ σ_k(n) qⁿ`.
* `TateCurve.E₄ R`, `TateCurve.E₆ R`: the Eisenstein series `1 + 240 s₃` and `1 - 504 s₅`
  as power series in `q`.
* `TateCurve.curve R`: the Tate curve over `R⟦X⟧`, with `X` playing `q`.

## Main statements

* `TateCurve.curve_c₄`, `TateCurve.curve_c₆`: `c₄ = E₄` and `c₆ = -E₆`.
* `TateCurve.curve_c_relation`: `1728 Δ = E₄³ - E₆²`.
* `TateCurve.curve_Δ`: `Δ = q · u` with `u` a unit of `R⟦q⟧` (`Δ = q - 24q² + ⋯`).
* `TateCurve.curve_map_constantCoeff`: setting `q = 0` gives the nodal cubic `Y² + XY = X³`,
  whose tangents at the node, `Y = 0` and `Y = -X`, are defined over `R`.
* `TateCurve.hasSplitMultiplicativeReduction`: over a field `k`, the Tate curve over the
  Laurent series field `k⸨X⸩` has split multiplicative reduction, and `Δ` has valuation exactly
  that of `q` (`TateCurve.X_dvd_Δ`, `TateCurve.not_X_sq_dvd_Δ`).

## Not yet here

* The product formula `Δ = q ∏ (1 - qⁿ)²⁴` (equivalent to Jacobi's and Ramanujan's identity for
  the discriminant modular form).
* The `j`-invariant `j = 1/q + 744 + 196884 q + ⋯` as a Laurent series.
* Tate's uniformisation `K^× / q^ℤ ≅ E_q(K)` over a complete non-archimedean field.

## Implementation notes

The Lean FLT project (Buzzard et al.) has the same series over `ℤ⟦q⟧` (`TateCurve.a₄Formal`,
`TateCurve.a₆Formal`), the Tate curve over a non-archimedean field at a point `|q| < 1`, and the
Weierstrass identity behind the uniformisation (`TateCurve.weierstrass_equation`). This file is
the formal curve over `R⟦q⟧` itself, with its invariants and reduction type; the two should be
merged before either goes to Mathlib.

## References

* J. Silverman, *Advanced Topics in the Arithmetic of Elliptic Curves*, chapter V.
* J. Tate, *A review of non-archimedean elliptic functions* (1995).

## Tags

elliptic curve, Tate curve, split multiplicative reduction, Eisenstein series
-/

@[expose] public section

open PowerSeries
open scoped LaurentSeries

namespace TateCurve

/-! ### The integrality of `a₆` -/

theorem twelve_dvd (d : ℕ) : 12 ∣ 5 * d ^ 3 + 7 * d ^ 5 := by
  rw [Nat.dvd_iff_mod_eq_zero, Nat.add_mod, Nat.mul_mod, Nat.pow_mod, Nat.mul_mod 7,
    Nat.pow_mod d 5]
  have h := Nat.mod_lt d (show 12 > 0 by norm_num)
  generalize d % 12 = r at h ⊢
  interval_cases r <;> norm_num

theorem twelve_dvd_sigma (n : ℕ) :
    12 ∣ 5 * ArithmeticFunction.sigma 3 n + 7 * ArithmeticFunction.sigma 5 n := by
  simp only [ArithmeticFunction.sigma_apply, Finset.mul_sum, ← Finset.sum_add_distrib]
  exact Finset.dvd_sum fun d _ ↦ twelve_dvd d

/-- `(5σ₃(n) + 7σ₅(n)) / 12`, a natural number by `twelve_dvd_sigma`. -/
def t (n : ℕ) : ℕ := (5 * ArithmeticFunction.sigma 3 n + 7 * ArithmeticFunction.sigma 5 n) / 12

theorem twelve_mul_t (n : ℕ) :
    12 * t n = 5 * ArithmeticFunction.sigma 3 n + 7 * ArithmeticFunction.sigma 5 n :=
  Nat.mul_div_cancel' (twelve_dvd_sigma n)

@[simp] theorem t_zero : t 0 = 0 := by simp [t]

@[simp] theorem t_one : t 1 = 1 := by simp [t, ArithmeticFunction.sigma_apply]

/-! ### The coefficients -/

variable (R : Type*) [CommRing R]

/-- The power series `∑_{n ≥ 1} σ_k(n) qⁿ`. -/
noncomputable def s (k : ℕ) : R⟦X⟧ := mk fun n ↦ (ArithmeticFunction.sigma k n : R)

/-- The Eisenstein series `E₄ = 1 + 240 ∑ σ₃(n) qⁿ`. -/
noncomputable def E₄ : R⟦X⟧ := 1 + 240 * s R 3

/-- The Eisenstein series `E₆ = 1 - 504 ∑ σ₅(n) qⁿ`. -/
noncomputable def E₆ : R⟦X⟧ := 1 - 504 * s R 5

/-- The coefficient `a₄ = -5 s₃` of the Tate curve. -/
noncomputable def a₄ : R⟦X⟧ := -5 * s R 3

/-- The coefficient `a₆ = -(5 s₃ + 7 s₅) / 12` of the Tate curve, defined coefficientwise. -/
noncomputable def a₆ : R⟦X⟧ := -mk fun n ↦ (t n : R)

/-- **The Tate curve** `Y² + XY = X³ + a₄ X + a₆` over `R⟦q⟧`. -/
noncomputable def curve : WeierstrassCurve R⟦X⟧ := ⟨1, 0, 0, a₄ R, a₆ R⟩

theorem coeff_s (k n : ℕ) : coeff n (s R k) = (ArithmeticFunction.sigma k n : R) := coeff_mk _ _

@[simp] theorem constantCoeff_s (k : ℕ) : constantCoeff (s R k) = 0 := by
  simp [s, constantCoeff_mk]

theorem coeff_numeral_mul (m : ℕ) [m.AtLeastTwo] (n : ℕ) (f : R⟦X⟧) :
    coeff n ((OfNat.ofNat m : R⟦X⟧) * f) = (OfNat.ofNat m : R) * coeff n f := by
  rw [← map_ofNat (C (R := R)) m, coeff_C_mul]

/-- The defining relation of `a₆`: `12 a₆ = -(5 s₃ + 7 s₅)`. -/
theorem twelve_mul_a₆ : 12 * a₆ R = -(5 * s R 3 + 7 * s R 5) := by
  ext n
  rw [coeff_numeral_mul, map_neg, map_add, coeff_numeral_mul, coeff_numeral_mul, a₆, map_neg,
    coeff_mk, coeff_s, coeff_s]
  have h := congrArg (Nat.cast (R := R)) (twelve_mul_t n)
  push_cast at h
  linear_combination -h

@[simp] theorem constantCoeff_a₄ : constantCoeff (a₄ R) = 0 := by simp [a₄]

@[simp] theorem constantCoeff_a₆ : constantCoeff (a₆ R) = 0 := by simp [a₆, constantCoeff_mk]

theorem coeff_one_a₆ : coeff 1 (a₆ R) = -1 := by simp [a₆, coeff_mk]

/-! ### The invariants -/

theorem curve_b₂ : (curve R).b₂ = 1 := by simp [curve, WeierstrassCurve.b₂]

theorem curve_b₄ : (curve R).b₄ = 2 * a₄ R := by simp [curve, WeierstrassCurve.b₄]

theorem curve_b₆ : (curve R).b₆ = 4 * a₆ R := by simp [curve, WeierstrassCurve.b₆]

/-- `c₄` of the Tate curve is the Eisenstein series `E₄`. -/
theorem curve_c₄ : (curve R).c₄ = E₄ R := by
  simp only [WeierstrassCurve.c₄, curve_b₂, curve_b₄, E₄, a₄]
  ring

/-- `c₆` of the Tate curve is `-E₆`. -/
theorem curve_c₆ : (curve R).c₆ = -E₆ R := by
  simp only [WeierstrassCurve.c₆, curve_b₂, curve_b₄, curve_b₆, E₆, a₄]
  linear_combination (-72 : R⟦X⟧) * twelve_mul_a₆ R

/-- `1728 Δ = E₄³ - E₆²` for the Tate curve. -/
theorem curve_c_relation : 1728 * (curve R).Δ = E₄ R ^ 3 - E₆ R ^ 2 := by
  rw [(curve R).c_relation, curve_c₄, curve_c₆]
  ring

theorem isUnit_E₄ : IsUnit (E₄ R) := by
  rw [isUnit_iff_constantCoeff]
  simp [E₄]

/-- **The discriminant of the Tate curve is `q` times a unit**: `Δ = q - 24 q² + ⋯`. -/
theorem curve_Δ : ∃ u : R⟦X⟧, IsUnit u ∧ (curve R).Δ = X * u := by
  obtain ⟨A, hA⟩ := X_dvd_iff.2 (constantCoeff_a₄ R)
  obtain ⟨B, hB⟩ := X_dvd_iff.2 (constantCoeff_a₆ R)
  have hB1 : constantCoeff B = -1 := by
    have h := coeff_one_a₆ R
    rw [hB, show (1 : ℕ) = 0 + 1 from rfl, coeff_succ_X_mul,
      coeff_zero_eq_constantCoeff_apply] at h
    exact h
  refine ⟨-B + X * (A ^ 2 - 64 * X * A ^ 3 - 432 * B ^ 2 + 72 * A * B), ?_, ?_⟩
  · rw [isUnit_iff_constantCoeff]
    simp [hB1]
  · rw [WeierstrassCurve.Δ, curve_b₂, curve_b₄, curve_b₆]
    simp only [WeierstrassCurve.b₈, curve, hA, hB]
    ring

theorem X_dvd_Δ : X ∣ (curve R).Δ := by
  obtain ⟨u, -, h⟩ := curve_Δ R
  exact ⟨u, h⟩

/-- `q²` does not divide `Δ`: the discriminant has order exactly one in `q`. -/
theorem not_X_sq_dvd_Δ [Nontrivial R] : ¬X ^ 2 ∣ (curve R).Δ := by
  obtain ⟨u, hu, hΔ⟩ := curve_Δ R
  rintro ⟨g, hg⟩
  have h1 := congrArg (coeff (0 + 1)) hg
  rw [hΔ, coeff_succ_X_mul, pow_two, mul_assoc, coeff_succ_X_mul, coeff_zero_X_mul,
    coeff_zero_eq_constantCoeff_apply] at h1
  rw [isUnit_iff_constantCoeff, h1] at hu
  exact not_isUnit_zero hu

/-- Setting `q = 0` gives the nodal cubic `Y² + XY = X³`. -/
theorem curve_map_constantCoeff :
    (curve R).map (constantCoeff (R := R)) = ⟨1, 0, 0, 0, 0⟩ := by
  ext <;> simp [curve]

/-! ### Split multiplicative reduction over `k⸨X⸩` -/

section Reduction

open WeierstrassCurve IsDiscreteValuationRing IsDedekindDomain.HeightOneSpectrum

/-- A Weierstrass equation with integral coefficients and `v(c₄) = 0` is minimal.
This is `WeierstrassCurve.isMinimal_of_valuation_c₄_eq_one` from the Lean FLT project
(Kevin Buzzard, William Coram, Claude; Apache 2.0), repeated here until it reaches Mathlib. -/
theorem isMinimal_of_valuation_c₄_eq_one {A : Type*} [CommRing A] [IsDomain A]
    [IsDiscreteValuationRing A] {K : Type*} [Field K] [Algebra A K] [IsFractionRing A K]
    (W : WeierstrassCurve K) [IsIntegral A W]
    (hc₄ : valuation K (maximalIdeal A) W.c₄ = 1) : IsMinimal A W := by
  refine ⟨⟨by simpa using ‹IsIntegral A W›, ?_⟩⟩
  intro C hC _
  simp only [one_smul, ← Subtype.coe_le_coe, valuation_Δ_aux_eq_of_isIntegral A (C • W),
    valuation_Δ_aux_eq_of_isIntegral A W]
  have hint : valuation K (maximalIdeal A) (C • W).c₄ ≤ 1 := by
    simpa [← integralModel_c₄_eq A (C • W)] using valuation_le_one _ _
  rw [variableChange_c₄, map_mul, map_pow, hc₄, mul_one] at hint
  simpa [variableChange_Δ, map_mul, map_pow] using mul_le_of_le_one_left'
    (pow_le_one' ((pow_le_one_iff (by norm_num)).mp hint) 12)

variable (k : Type*) [Field k]

instance : IsIntegral k⟦X⟧ ((curve k).baseChange k⸨X⸩) := ⟨⟨curve k, rfl⟩⟩

theorem integralModel_curve : integralModel k⟦X⟧ ((curve k).baseChange k⸨X⸩) = curve k := by
  have hinj := IsFractionRing.injective k⟦X⟧ k⸨X⸩
  apply WeierstrassCurve.ext
  · exact hinj (by rw [integralModel_a₁_eq k⟦X⟧ ((curve k).baseChange k⸨X⸩)]; rfl)
  · exact hinj (by rw [integralModel_a₂_eq k⟦X⟧ ((curve k).baseChange k⸨X⸩)]; rfl)
  · exact hinj (by rw [integralModel_a₃_eq k⟦X⟧ ((curve k).baseChange k⸨X⸩)]; rfl)
  · exact hinj (by rw [integralModel_a₄_eq k⟦X⟧ ((curve k).baseChange k⸨X⸩)]; rfl)
  · exact hinj (by rw [integralModel_a₆_eq k⟦X⟧ ((curve k).baseChange k⸨X⸩)]; rfl)

variable {k} in
theorem valuation_algebraMap_eq_one_iff (r : k⟦X⟧) :
    valuation k⸨X⸩ (maximalIdeal k⟦X⟧) (algebraMap k⟦X⟧ k⸨X⸩ r) = 1 ↔ IsUnit r := by
  have h1 := valuation_le_one (K := k⸨X⸩) (maximalIdeal k⟦X⟧) r
  have h2 := valuation_lt_one_iff_mem (K := k⸨X⸩) (maximalIdeal k⟦X⟧) r
  have h3 : r ∈ (maximalIdeal k⟦X⟧).asIdeal ↔ ¬IsUnit r := by
    change r ∈ IsLocalRing.maximalIdeal k⟦X⟧ ↔ _
    rw [IsLocalRing.mem_maximalIdeal, mem_nonunits_iff]
  constructor
  · intro h
    by_contra hu
    exact (h2.2 (h3.2 hu)).ne h
  · intro hu
    exact le_antisymm h1 (not_lt.1 fun hlt ↦ h3.1 (h2.1 hlt) hu)

variable {k} in
theorem residue_eq_zero {r : k⟦X⟧} (h : constantCoeff r = 0) :
    IsLocalRing.residue k⟦X⟧ r = 0 := by
  rw [IsLocalRing.residue_eq_zero_iff, IsLocalRing.mem_maximalIdeal, mem_nonunits_iff,
    isUnit_iff_constantCoeff, h]
  exact not_isUnit_zero

/-- **Over a field `k`, the Tate curve over `k⸨q⸩` has split multiplicative reduction.** -/
theorem hasSplitMultiplicativeReduction :
    HasSplitMultiplicativeReduction k⟦X⟧ ((curve k).baseChange k⸨X⸩) := by
  have hc : valuation k⸨X⸩ (maximalIdeal k⟦X⟧) ((curve k).baseChange k⸨X⸩).c₄ = 1 := by
    rw [baseChange, map_c₄, valuation_algebraMap_eq_one_iff, curve_c₄]
    exact isUnit_E₄ k
  have hΔ : valuation k⸨X⸩ (maximalIdeal k⟦X⟧) ((curve k).baseChange k⸨X⸩).Δ < 1 := by
    rw [baseChange, map_Δ, (valuation_le_one (K := k⸨X⸩) _ _).lt_iff_ne,
      Ne, valuation_algebraMap_eq_one_iff, isUnit_iff_constantCoeff,
      ← coeff_zero_eq_constantCoeff_apply]
    obtain ⟨u, -, h⟩ := curve_Δ k
    rw [h, coeff_zero_X_mul]
    exact not_isUnit_zero
  have := isMinimal_of_valuation_c₄_eq_one (A := k⟦X⟧) _ hc
  have : HasMultiplicativeReduction k⟦X⟧ ((curve k).baseChange k⸨X⸩) := ⟨hΔ, hc⟩
  refine ⟨?_⟩
  show Polynomial.Splits _
  rw [integralModel_curve, IsLocalRing.ResidueField.algebraMap_eq]
  have ha₄ := residue_eq_zero (constantCoeff_a₄ k)
  have ha₆ := residue_eq_zero (constantCoeff_a₆ k)
  have key : (Polynomial.C (curve k).c₄ * Polynomial.X ^ 2 +
      Polynomial.C ((curve k).a₁ * (curve k).c₄) * Polynomial.X -
      Polynomial.C (54 * (curve k).b₆ - 3 * (curve k).b₂ * (curve k).b₄ +
        (curve k).a₂ * (curve k).c₄)).map (IsLocalRing.residue k⟦X⟧) =
      Polynomial.C (IsLocalRing.residue k⟦X⟧ (curve k).c₄) *
        (Polynomial.X * (Polynomial.X + Polynomial.C 1)) := by
    simp only [Polynomial.map_sub, Polynomial.map_add, Polynomial.map_mul, Polynomial.map_pow,
      Polynomial.map_C, Polynomial.map_X, curve_b₂, curve_b₄, curve_b₆, map_mul, map_add,
      map_sub, map_ofNat, map_one, ha₄, ha₆]
    simp only [curve, map_one, map_zero, one_mul, zero_mul, add_zero]
    ring
  rw [key]
  exact (Polynomial.Splits.X.mul (Polynomial.Splits.X_add_C 1)).C_mul _

end Reduction

end TateCurve
