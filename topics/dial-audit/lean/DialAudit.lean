/-
  FusedAlgebra.Topics.DialAudit: which dials are fundamental, which derived (DRAFT, Nov 26 2026)
  ==========================================================================================

  The exact identities behind topics/dial-audit/. Core Lean 4 only, no `sorry`.
  The unit change of the plane law itself is already proved as LandscapeMap.units_num and
  LandscapeMap.units_den (k ↦ λk, d ↦ λ²d); this file adds what the audit leans on.

  1. Unit-free numbers of the (k, d) plane: d/k² and Δ/k² (stated cross-multiplied).
  2. The odds absorb k: three odds fused at dial k are the k = 1 fusion of the odds k·t, over k.
  3. The order dial's counting shadow: the q-integer [n]_q is the n-fold ◯[q − 1] power of 1.
  4. A unit change moves the order dial's swap constant not at all (q is unit-free), while it
     moves k: the contrast the audit draws between q and k.
-/
import FusedAlgebra.Base

namespace CircleOp
namespace DialAudit

open Lean.Grind

section Ring
variable {R : Type} [CommRing R]

/-! ## 1. Unit-free numbers of the (k, d) plane -/

/-- Under the unit change k ↦ λk, d ↦ λ²d the ratio d/k² is unchanged (cross-multiplied form). -/
theorem ratio_unit_free (k d lam : R) : (lam * lam * d) * (k * k) = d * ((lam * k) * (lam * k)) := by
  grind

/-- The discriminant Δ = k² − 4d scales like unit², so Δ/k² is unit-free and Δ itself is not. -/
theorem discriminant_scales (k d lam : R) :
    (lam * k) * (lam * k) - 4 * (lam * lam * d) = lam * lam * (k * k - 4 * d) := by
  grind

/-- The curvature reading −Δ/4 = d − k²/4, written without division: 4·(d − k²/4) = 4d − k². -/
theorem curvature_reading (k d : R) : 4 * d - k * k = -(k * k - 4 * d) := by
  grind

/-! ## 2. The odds absorb the homothety dial -/

/-- Three odds fused at dial k, scaled by k, are the k = 1 fusion of the odds k·t. -/
theorem odds_absorb_k (k t₁ t₂ t₃ : R) :
    k * circK k (circK k t₁ t₂) t₃ = circ (circ (k * t₁) (k * t₂)) (k * t₃) := by
  unfold circK circ; grind

/-- The same for the lift: 1 + k·(fusion) is the product of the three factors 1 + k·t. -/
theorem odds_absorb_k_lift (k t₁ t₂ t₃ : R) :
    1 + k * circK k (circK k t₁ t₂) t₃ = (1 + k * t₁) * (1 + k * t₂) * (1 + k * t₃) := by
  unfold circK; grind

/-! ## 3. The order dial's counting shadow -/

/-- The q-integer [n]_q = 1 + q + … + q^(n−1). -/
def qint (q : R) : Nat → R
  | 0 => 0
  | n + 1 => qint q n + q ^ n

/-- The n-fold ◯[k] power of 1. -/
def powOne (k : R) : Nat → R
  | 0 => 0
  | n + 1 => circK k (powOne k n) 1

theorem qint_pow (q : R) (n : Nat) : 1 + (q - 1) * qint q n = q ^ n := by
  induction n with
  | zero => simp only [qint]; grind
  | succ n ih =>
    simp only [qint]
    have hp : q ^ (n + 1) = q ^ n * q := Lean.Grind.Semiring.pow_succ q n
    grind

/-- [n]_q is the n-fold ◯[q − 1] power of 1: the q of counting is the homothety dial k = q − 1. -/
theorem qint_is_circ_power (q : R) (n : Nat) : powOne (q - 1) n = qint q n := by
  induction n with
  | zero => rfl
  | succ n ih =>
    simp only [powOne, qint, ih]
    have h := qint_pow q n
    unfold circK
    grind

/-! ## 4. q is unit-free, k is not -/

/-- If two moves obey y·x = q·x·y, rescaling them (by scalars α, β) keeps the same q. -/
theorem swap_constant_unit_free (q x y α β : R) (h : y * x = q * (x * y)) :
    (β * y) * (α * x) = q * ((α * x) * (β * y)) := by
  grind

/-- k is not: in units λ times bigger the fused law's dial changes from k to k·λ
(the repo's units_change_weight, restated here next to the q statement it contrasts with). -/
theorem homothety_moves_with_units (k lam a b : R) :
    circK k (lam * a) (lam * b) = lam * circK (k * lam) a b := by
  unfold circK; grind

end Ring

end DialAudit
end CircleOp
