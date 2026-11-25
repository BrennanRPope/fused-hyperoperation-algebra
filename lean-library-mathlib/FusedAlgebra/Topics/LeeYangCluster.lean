module

public import FusedAlgebra.Basic

/-!
# The Beraha cluster test

The colour-count reading of the dial (`k = Q - 1`) predicts special dial values
`k_n = 1 + 2 * cos (2 * π / n)`: 1.618, 2, 2.247, 2.414, 2.532, 2.618, …, 2.732 for
`n = 5, …, 12`. Two facts are behind the cluster test on the landscape.

1. Every predicted value except `n = 6` (`k = 2`) is a root of an integer polynomial with no
   integer root (`beraha5_not_int`, …, `beraha12_not_int`), so a law at a whole-number dial value
   can only ever hit `k = 2`. Each value is in fact irrational (the rational root theorem; in
   mathlib the general tool is `IsIntegrallyClosed` of `ℤ`, e.g. `irrational_nrt_of_notint_nrt`
   for pure roots); the integer statement is what the test uses, and that is what is proved.
2. A change of units `a ↦ l * a` turns weight `l * k` into weight `k` (`wAdd_units`). So on the
   `q = 0` axis a positive dial value is not a property of a law: `2 ×` leverage at `k = 2` is
   plain compounding at `k = 1` in levered units (`leverage_is_compounding`).
-/

@[expose] public section

namespace LeeYangCluster

open MultAdd

/-! ### Beraha numbers are not integers -/

/-- `n = 5`: `k` is the golden ratio, `k ^ 2 = k + 1`, which has no integer solution. -/
theorem beraha5_not_int (k : ℤ) : k ^ 2 ≠ k + 1 := by
  intro h
  have hk := Int.natAbs_le_of_dvd_ne_zero (Dvd.intro (k - 1) (by linear_combination h))
    (one_ne_zero : (1 : ℤ) ≠ 0)
  have : -1 ≤ k ∧ k ≤ 1 := by omega
  obtain ⟨h1, h2⟩ := this
  interval_cases k <;> omega

/-- `n = 7`: `k = 1 + 2 * cos (2 * π / 7)` is a root of `k ^ 3 - 2 * k ^ 2 - k + 1`, which has no
integer root. -/
theorem beraha7_not_int (k : ℤ) : k ^ 3 - 2 * k ^ 2 - k + 1 ≠ 0 := by
  intro h
  have hk := Int.natAbs_le_of_dvd_ne_zero
    (Dvd.intro (2 * k + 1 - k ^ 2) (by linear_combination -h)) (one_ne_zero : (1 : ℤ) ≠ 0)
  have : -1 ≤ k ∧ k ≤ 1 := by omega
  obtain ⟨h1, h2⟩ := this
  interval_cases k <;> omega

/-- `n = 8`: `k = 1 + √2`, `k ^ 2 = 2 * k + 1`, which has no integer solution. -/
theorem beraha8_not_int (k : ℤ) : k ^ 2 ≠ 2 * k + 1 := by
  intro h
  have hk := Int.natAbs_le_of_dvd_ne_zero (Dvd.intro (k - 2) (by linear_combination h))
    (one_ne_zero : (1 : ℤ) ≠ 0)
  have : -1 ≤ k ∧ k ≤ 1 := by omega
  obtain ⟨h1, h2⟩ := this
  interval_cases k <;> omega

/-- `n = 9`: `k = 1 + 2 * cos (2 * π / 9)` is a root of `k ^ 3 - 3 * k ^ 2 + 3`, which has no
integer root. -/
theorem beraha9_not_int (k : ℤ) : k ^ 3 - 3 * k ^ 2 + 3 ≠ 0 := by
  intro h
  have hk := Int.natAbs_le_of_dvd_ne_zero
    (Dvd.intro (k * (3 - k)) (by linear_combination -h)) (by norm_num : (3 : ℤ) ≠ 0)
  have : -3 ≤ k ∧ k ≤ 3 := by omega
  obtain ⟨h1, h2⟩ := this
  interval_cases k <;> omega

/-- `n = 10`: `k` is `1` plus the golden ratio, `k ^ 2 = 3 * k - 1`, which has no integer
solution. -/
theorem beraha10_not_int (k : ℤ) : k ^ 2 ≠ 3 * k - 1 := by
  intro h
  have hk := Int.natAbs_le_of_dvd_ne_zero (Dvd.intro (3 - k) (by linear_combination -h))
    (one_ne_zero : (1 : ℤ) ≠ 0)
  have : -1 ≤ k ∧ k ≤ 1 := by omega
  obtain ⟨h1, h2⟩ := this
  interval_cases k <;> omega

/-- `n = 12`: `k = 1 + √3`, `k ^ 2 = 2 * k + 2`, which has no integer solution. -/
theorem beraha12_not_int (k : ℤ) : k ^ 2 ≠ 2 * k + 2 := by
  intro h
  have hk := Int.natAbs_le_of_dvd_ne_zero (Dvd.intro (k - 2) (by linear_combination h))
    (by norm_num : (2 : ℤ) ≠ 0)
  have : -2 ≤ k ∧ k ≤ 2 := by omega
  obtain ⟨h1, h2⟩ := this
  interval_cases k <;> omega

/-! ### Units -/

variable {R : Type*}

/-- Units: rescaling the variable by `l` turns weight `l * k` into weight `k`. -/
theorem wAdd_units [CommSemiring R] (l k a b : R) :
    l * wAdd (l * k) a b = wAdd k (l * a) (l * b) := by
  unfold wAdd; ring

/-- `2 ×` leverage: weight `2` on unlevered returns is plain mult-add on levered returns. -/
theorem leverage_is_compounding [CommSemiring R] (a b : R) :
    2 * wAdd 2 a b = (2 * a) ✳ (2 * b) :=
  mul_wAdd 2 a b

/-- The same move off the axis: the `q`-coordinate scales by `l ^ 2`, so `k ^ 2 / q` is what
stays. -/
theorem units_bend [CommRing R] (l q a b : R) :
    1 - (l ^ 2 * q) * a * b = 1 - q * (l * a) * (l * b) := by
  ring

end LeeYangCluster
