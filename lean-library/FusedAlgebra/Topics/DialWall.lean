/-
  Dial wall: exact facts behind the wall at Re s = 0 for Z_k(s) = ∏ (1 - k p^-s)^(-1/k).

  Z_k(s) = ∏_N ζ(Ns)^(M(k,N)/k), where M(k,N) counts aperiodic necklaces.
  This file proves, for N ≤ 6:
  * the necklace (word) identities at orders z^4, z^5, z^6;
  * k = +1 kills every exponent from N = 2 on, and k = -1 every one from N = 3 on
    (that is why only k = ±1 give finitely many ζ factors);
  * at prime length the exponent vanishes only at k ∈ {0, ±1};
  * the cancellation found numerically: at k = -1/2, M1 + M2 + M3 = 0 and 3 M1 + 4 M2 = 0,
    so rational k need more than the generic argument;
  * the copies points k = -1/m give N m^N M(-1/m, N) = (-1)^N mod m for N = 2, 3,
    and an integer of the form ±1 + m·A is never 0 when m ≥ 2 (rational root step);
  * on the critical line, ρ/N has real part 1/(2N), so N' ρ/N lies on it only for N' = N.
  Core Lean only; no sorry.
-/
import FusedAlgebra.Topics.DialSurvey

namespace CircleOp
namespace DialWall
open DialSurvey

def M4 (k : Rat) : Rat := (k * k * k * k - k * k) / 4
def M5 (k : Rat) : Rat := (k * k * k * k * k - k) / 5
def M6 (k : Rat) : Rat := (k * k * k * k * k * k - k * k * k - k * k + k) / 6

/-- Words of length 4 split by period: 4·M4 + 2·M2 + 1·M1 = k^4. -/
theorem words_4 (k : Rat) : 4 * M4 k + 2 * M2 k + M1 k = k * k * k * k := by
  unfold M4 M2 M1; grind

theorem words_5 (k : Rat) : 5 * M5 k + M1 k = k * k * k * k * k := by
  unfold M5 M1; grind

theorem words_6 (k : Rat) :
    6 * M6 k + 3 * M3 k + 2 * M2 k + M1 k = k * k * k * k * k * k := by
  unfold M6 M3 M2 M1; grind

/-- k = +1: every necklace count from length 2 on is 0 (only ζ(s) is left). -/
theorem bose_all_zero : M2 1 = 0 ∧ M3 1 = 0 ∧ M4 1 = 0 ∧ M5 1 = 0 ∧ M6 1 = 0 := by
  unfold M2 M3 M4 M5 M6; grind

/-- k = -1: from length 3 on every count is 0; length 2 gives 1 (the ζ(2s) of ζ(s)/ζ(2s)). -/
theorem fermi_from_three :
    M2 (-1) = 1 ∧ M3 (-1) = 0 ∧ M4 (-1) = 0 ∧ M5 (-1) = 0 ∧ M6 (-1) = 0 := by
  unfold M2 M3 M4 M5 M6; grind

theorem M2_zero_iff (k : Rat) : M2 k = 0 ↔ k = 0 ∨ k = 1 := by
  unfold M2
  constructor
  · intro h
    by_cases h0 : k = 0
    · exact Or.inl h0
    · right
      have : k * (k - 1) = 0 := by grind
      have : k - 1 = 0 := by
        rcases Rat.mul_eq_zero.mp this with h1 | h1
        · exact absurd h1 h0
        · exact h1
      grind
  · rintro (h | h) <;> subst h <;> grind

theorem M3_zero_iff (k : Rat) : M3 k = 0 ↔ k = 0 ∨ k = 1 ∨ k = -1 := by
  unfold M3
  constructor
  · intro h
    have e : k * ((k - 1) * (k + 1)) = 0 := by grind
    rcases Rat.mul_eq_zero.mp e with h1 | h1
    · exact Or.inl h1
    · rcases Rat.mul_eq_zero.mp h1 with h2 | h2
      · right; left; grind
      · right; right; grind
  · rintro (h | h | h) <;> subst h <;> grind

theorem M5_zero_iff (k : Rat) : M5 k = 0 ↔ k = 0 ∨ k = 1 ∨ k = -1 := by
  unfold M5
  constructor
  · intro h
    have e : k * ((k - 1) * ((k + 1) * (k * k + 1))) = 0 := by grind
    rcases Rat.mul_eq_zero.mp e with h1 | h1
    · exact Or.inl h1
    · rcases Rat.mul_eq_zero.mp h1 with h2 | h2
      · right; left; grind
      · rcases Rat.mul_eq_zero.mp h2 with h3 | h3
        · right; right; grind
        · have : 0 ≤ k * k := by
            by_cases hk : 0 ≤ k
            · exact Rat.mul_nonneg hk hk
            · have hk' : 0 ≤ -k := by grind
              have := Rat.mul_nonneg hk' hk'
              grind
          grind
  · rintro (h | h | h) <;> subst h <;> grind

/-- The cancellation: three contributors with one zero each cancel at k = -1/2. -/
theorem cancel_half : M1 (-1/2) + M2 (-1/2) + M3 (-1/2) = 0 := by
  unfold M1 M2 M3; grind

theorem cancel_half' : 3 * M1 (-1/2) + 4 * M2 (-1/2) = 0 := by
  unfold M1 M2; grind

/-- k = 1/2: the exponent at ρ/2 is M2/k = -1/4, not a whole number, so a branch point. -/
theorem half_exponent : M2 (1/2) / (1/2) = -1/4 := by
  unfold M2; grind

/-- Copies points, length 2: 2 m² M2(-1/m) = 1 + m. -/
theorem copies_len2 (m : Rat) (hm : m ≠ 0) : 2 * (m * m) * M2 (-1 / m) = 1 + m := by
  unfold M2; grind

/-- Copies points, length 3: 3 m³ M3(-1/m) = m² - 1. -/
theorem copies_len3 (m : Rat) (hm : m ≠ 0) : 3 * (m * m * m) * M3 (-1 / m) = m * m - 1 := by
  unfold M3; grind

/-- Rational root step: ±1 + m·A is never 0 for a whole m ≥ 2. -/
theorem unit_plus_multiple_ne_zero (m A : Int) (hm : 2 ≤ m) (s : Int) (hs : s = 1 ∨ s = -1) :
    s + m * A ≠ 0 := by
  intro h
  rcases Int.lt_trichotomy A 0 with hA | hA | hA
  · have : m * A ≤ m * (-1) := Int.mul_le_mul_of_nonneg_left (by omega) (by omega)
    omega
  · subst hA; omega
  · have : m * 1 ≤ m * A := Int.mul_le_mul_of_nonneg_left (by omega) (by omega)
    omega

/-- On the critical line ρ/N has real part 1/(2N); N'·ρ/N is back on the line only for N' = N. -/
theorem one_contributor (N N' : Nat) (hN : 0 < N) :
    ((N' : Rat) / (2 * N) = 1 / 2) → N' = N := by
  intro h
  have hN' : (N : Rat) ≠ 0 := by
    intro h0
    exact absurd (Rat.natCast_inj.mp (h0.trans (by rfl : (0:Rat) = ((0:Nat):Rat)))) (by omega)
  have : (N' : Rat) = N := by grind
  exact Rat.natCast_inj.mp this

end DialWall
end CircleOp
