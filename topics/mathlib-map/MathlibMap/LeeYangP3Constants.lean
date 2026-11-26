import Mathlib

/-!
# Lee–Yang at roundness 3: the certificate's constants (DRAFT)

The p = 3 certificate (`p3_certificate.py`) moves the shape parameter c of
G_c(s) = ∫₀¹ (1 - x³)^c cos(s x) dx across boxes of width 1/100. Its error terms need bounds on

* L2 = ∫₀¹ log²(1 - x³) dx   (second c-derivative), and
* L3 = ∫₀¹ |log(1 - x³)|³ dx (third c-derivative).

Earlier the certificate used the decimals 0.75 and 2.11 with a safety margin. It now uses the whole numbers
**L2 ≤ 2 = 2!** and **L3 ≤ 6 = 3!**, and still passes. The reason is one line: on [0, 1], 1 - x ≤ 1 - x³, so
|log(1 - x³)| ≤ |log(1 - x)|, and ∫₀¹ |log(1 - x)|^k dx = ∫₀¹ (-log u)^k du = k! (Euler's integral for Γ(k + 1)).

What Lean checks here:
* the pointwise comparison, for every x in [0, 1] and every power k;
* the integral comparison, given Euler's integral and integrability as hypotheses (classical; the change of
  variables to Mathlib's `Real.Gamma_eq_integral` is not done here);
* the exact fractions and whole numbers the certificate's error terms reduce to.

For reference (not used by the proof, not checked in Lean): from the positive series
-log(1 - x³) = Σ x^(3k)/k, rigorous partial sums plus tail bounds give
747/1000 < L2 < 748/1000 and 1051/500 < L3 < 263/125.
-/

namespace LeeYangP3Constants

/-- On [0, 1], 1 - x ≤ 1 - x³, because x - x³ = x(1 - x) + x²(1 - x). -/
theorem one_sub_le_one_sub_cube {x : ℝ} (h0 : 0 ≤ x) (h1 : x ≤ 1) : 1 - x ≤ 1 - x ^ 3 := by
  nlinarith [mul_nonneg h0 (sub_nonneg.mpr h1), mul_nonneg (mul_nonneg h0 h0) (sub_nonneg.mpr h1)]

/-- The log comparison behind both bounds. At x = 1 both sides are |log 0| = 0. -/
theorem abs_log_le {x : ℝ} (h0 : 0 ≤ x) (h1 : x ≤ 1) :
    |Real.log (1 - x ^ 3)| ≤ |Real.log (1 - x)| := by
  rcases eq_or_lt_of_le h1 with h | hlt
  · subst h; simp
  · have hpos : 0 < 1 - x := by linarith
    have hle := one_sub_le_one_sub_cube h0 h1
    have hcube : 0 ≤ x ^ 3 := pow_nonneg h0 3
    have hA : Real.log (1 - x ^ 3) ≤ 0 := Real.log_nonpos (by linarith) (by linarith)
    have hB : Real.log (1 - x) ≤ 0 := Real.log_nonpos hpos.le (by linarith)
    have hC : Real.log (1 - x) ≤ Real.log (1 - x ^ 3) := Real.log_le_log hpos hle
    rw [abs_of_nonpos hA, abs_of_nonpos hB]
    linarith

/-- The same comparison for every power k. -/
theorem pow_abs_log_le {x : ℝ} (h0 : 0 ≤ x) (h1 : x ≤ 1) (k : ℕ) :
    |Real.log (1 - x ^ 3)| ^ k ≤ |Real.log (1 - x)| ^ k := by
  have h := abs_log_le h0 h1
  gcongr

/-- ∫₀¹ |log(1 - x³)|^k ≤ k!, given Euler's integral ∫₀¹ |log(1 - x)|^k = k! and integrability. -/
theorem integral_le_factorial (k : ℕ)
    (hf : IntervalIntegrable (fun x : ℝ => |Real.log (1 - x ^ 3)| ^ k) MeasureTheory.volume 0 1)
    (hg : IntervalIntegrable (fun x : ℝ => |Real.log (1 - x)| ^ k) MeasureTheory.volume 0 1)
    (euler : ∫ x in (0 : ℝ)..1, |Real.log (1 - x)| ^ k = (k.factorial : ℝ)) :
    ∫ x in (0 : ℝ)..1, |Real.log (1 - x ^ 3)| ^ k ≤ (k.factorial : ℝ) := by
  rw [← euler]
  exact intervalIntegral.integral_mono_on zero_le_one hf hg (fun x hx => pow_abs_log_le hx.1 hx.2 k)

/-- The two constants the certificate uses are these factorials: L2 ≤ 2! = 2 and L3 ≤ 3! = 6. -/
theorem factorials_used : (Nat.factorial 2 : ℝ) = 2 ∧ (Nat.factorial 3 : ℝ) = 6 := by
  norm_num [Nat.factorial]

/-- Near-zone error terms in c, as exact fractions: box half-width δ = 1/200 gives δ² L2 / 2 = 1/40000, and the
central-difference step η = 1/10¹² gives η² L3 / 6 = 1/10²⁴. -/
theorem near_zone_error_terms :
    ((1 : ℚ) / 200) ^ 2 * 2 / 2 = 1 / 40000 ∧ ((1 : ℚ) / 10 ^ 12) ^ 2 * 6 / 6 = 1 / 10 ^ 24 := by
  norm_num

/-- The far zone's remainder constants are whole numbers: 6!, 7!, 9! and 10!. -/
theorem far_zone_whole_numbers :
    Nat.factorial 6 = 720 ∧ Nat.factorial 7 = 5040 ∧ Nat.factorial 9 = 362880 ∧
      Nat.factorial 10 = 3628800 := by
  norm_num [Nat.factorial]

end LeeYangP3Constants
