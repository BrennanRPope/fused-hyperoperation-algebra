/-
  Fused Hyperoperation Algebra: the Beraha cluster test.
  Core Lean 4.37 only (no mathlib), no `sorry`. Imports FusedAlgebra.Base.
  Intended library path: FusedAlgebra/Topics/LeeYangCluster.lean (after LeeYang).
  LIBRARY NOTE: copied unchanged from lee-yang/lean/LeeYangCluster.lean. No exact twins in the
  library. `circK_units` is the ring-general relative of Landscape.units_num (Int, stated on num).

  The colour-count reading of the dial (k = Q − 1) predicts special dial values
  k_n = 1 + 2·cos(2π/n): 1.618, 2, 2.247, 2.414, 2.532, 2.618, ..., 2.732 for n = 5..12.
  Two facts behind the cluster test on the landscape:

  1. Every predicted value except n = 6 (k = 2) is a root of an integer polynomial with no
     integer root, so a law at a whole-number dial value can only ever hit k = 2.
     (Each value is in fact irrational; the integer statement is what the test uses.)
  2. A change of units a ↦ λa turns weight λk into weight k. So on the q = 0 axis a positive
     dial value is not a property of a law: 2× leverage at k = 2 is plain compounding at
     k = 1 in levered units.
-/
import FusedAlgebra.Base

namespace LeeYangCluster
open CircleOp

set_option linter.unusedSectionVars false

/-- If a·b = c with c ≠ 0, then |a| ≤ |c|. -/
theorem small_of_mul (a b c : Int) (hc : c ≠ 0) (h : a * b = c) :
    -c.natAbs ≤ a ∧ a ≤ c.natAbs := by
  have hd : a.natAbs ∣ c.natAbs := Int.natAbs_dvd_natAbs.mpr ⟨b, h.symm⟩
  have hpos : 0 < c.natAbs := Int.natAbs_pos.mpr hc
  have := Nat.le_of_dvd hpos hd
  omega

/-- n = 5: k = golden ratio, k² = k + 1, has no integer solution. -/
theorem beraha5_not_int (k : Int) : k * k ≠ k + 1 := by
  intro h
  have hm : k * (k - 1) = 1 := by rw [Int.mul_sub, h]; omega
  have ⟨h1, h2⟩ := small_of_mul k (k - 1) 1 (by decide) hm
  have : k = -1 ∨ k = 0 ∨ k = 1 := by simp at h1 h2; omega
  rcases this with rfl | rfl | rfl <;> simp at h

/-- n = 7: k = 1 + 2cos(2π/7) is a root of k³ − 2k² − k + 1, which has no integer root. -/
theorem beraha7_not_int (k : Int) : k * k * k - 2 * (k * k) - k + 1 ≠ 0 := by
  intro h
  have hm : k * (2 * k + 1 - k * k) = 1 := by
    have : k * (2 * k + 1 - k * k) = 2 * (k * k) + k - k * k * k := by
      simp only [Int.mul_sub, Int.mul_add, Int.mul_one]; grind
    omega
  have ⟨h1, h2⟩ := small_of_mul k _ 1 (by decide) hm
  have : k = -1 ∨ k = 0 ∨ k = 1 := by simp at h1 h2; omega
  rcases this with rfl | rfl | rfl <;> simp at h

/-- n = 8: k = 1 + √2, k² = 2k + 1, has no integer solution. -/
theorem beraha8_not_int (k : Int) : k * k ≠ 2 * k + 1 := by
  intro h
  have hm : k * (k - 2) = 1 := by rw [Int.mul_sub, h]; omega
  have ⟨h1, h2⟩ := small_of_mul k (k - 2) 1 (by decide) hm
  have : k = -1 ∨ k = 0 ∨ k = 1 := by simp at h1 h2; omega
  rcases this with rfl | rfl | rfl <;> simp at h

/-- n = 9: k = 1 + 2cos(2π/9) is a root of k³ − 3k² + 3, which has no integer root. -/
theorem beraha9_not_int (k : Int) : k * k * k - 3 * (k * k) + 3 ≠ 0 := by
  intro h
  have hm : k * (k * (3 - k)) = 3 := by
    have : k * (k * (3 - k)) = 3 * (k * k) - k * k * k := by
      simp only [Int.mul_sub]; grind
    omega
  have ⟨h1, h2⟩ := small_of_mul k _ 3 (by decide) hm
  have : k = -3 ∨ k = -2 ∨ k = -1 ∨ k = 0 ∨ k = 1 ∨ k = 2 ∨ k = 3 := by simp at h1 h2; omega
  rcases this with rfl | rfl | rfl | rfl | rfl | rfl | rfl <;> simp at h

/-- n = 10: k = 1 + golden ratio, k² = 3k − 1, has no integer solution. -/
theorem beraha10_not_int (k : Int) : k * k ≠ 3 * k - 1 := by
  intro h
  have hm : k * (3 - k) = 1 := by rw [Int.mul_sub, h]; omega
  have ⟨h1, h2⟩ := small_of_mul k (3 - k) 1 (by decide) hm
  have : k = -1 ∨ k = 0 ∨ k = 1 := by simp at h1 h2; omega
  rcases this with rfl | rfl | rfl <;> simp at h

/-- n = 12: k = 1 + √3, k² = 2k + 2, has no integer solution. -/
theorem beraha12_not_int (k : Int) : k * k ≠ 2 * k + 2 := by
  intro h
  have hm : k * (k - 2) = 2 := by rw [Int.mul_sub, h]; omega
  have ⟨h1, h2⟩ := small_of_mul k (k - 2) 2 (by decide) hm
  have : k = -2 ∨ k = -1 ∨ k = 0 ∨ k = 1 ∨ k = 2 := by simp at h1 h2; omega
  rcases this with rfl | rfl | rfl | rfl | rfl <;> simp at h

/-- Units: rescaling the variable by λ turns weight λ·k into weight k. -/
theorem circK_units {R : Type} [Lean.Grind.CommRing R] (lam k a b : R) :
    lam * circK (lam * k) a b = circK k (lam * a) (lam * b) := by
  unfold circK; grind

/-- 2× leverage: weight 2 on unlevered returns is the plain circle operation (called mult-add before 1.0) on levered returns. -/
theorem leverage_is_compounding {R : Type} [Lean.Grind.CommRing R] (a b : R) :
    2 * circK 2 a b = circK 1 (2 * a) (2 * b) := by
  have h := circK_units (2 : R) 1 a b
  rw [show (2 : R) * 1 = 2 by grind] at h
  exact h

/-- The same move off the axis: the q-coordinate scales by λ², so k²/q is what stays. -/
theorem units_bend {R : Type} [Lean.Grind.CommRing R] (lam q a b : R) :
    1 - (lam * lam * q) * a * b = 1 - q * (lam * a) * (lam * b) := by grind

end LeeYangCluster
