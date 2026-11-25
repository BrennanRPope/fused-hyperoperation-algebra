/-
  Dial survey: the exact facts behind the quick pass on the five dials
  ===================================================================

  Part of Fused Hyperoperation Algebra. Companion page and checks:
  math-theory/dial-survey/ (index.html, checks.py). Core Lean 4.37 only
  (no mathlib), no `sorry`. Imports Base only.

  1. The k dial. In units λ times bigger, weight k reads as weight k·λ
     (`units_change_weight`), so on the bare operation only the sign of k
     survives. Once each element carries its own dial, the dial values fix
     the units and k's size matters (part 3).
  2. Per-element dials are odds. With p = t/(1 + t), "in A or in B" for
     independent random subsets has odds t ◯ s (`odds_union`). The uniform
     dial t = 2^a − 1 adds exponents: (2^a − 1) ◯ (2^b − 1) = 2^(a+b) − 1
     (`exponent_ruler`), so the per-element dial and the fractional exponent
     are one dial, read through the circle-operation ruler (called mult-add ruler before 1.0).
  3. k acting on per-element dials. One element with weight k contributes
     (1 − k·z)^(−1/k); its z^e coefficient is kRise k e / e!, where
     kRise k e = ∏_{j<e} (1 + j·k). Then:
       k = −1: no element twice (`kRise_fermi`): subsets, the power set.
       k =  0: kRise = 1, the coefficient is 1/e! (`kRise_poisson`).
       k = +1: kRise = e! (`kRise_bose`): multisets.
       k = −1/m: m^e · kRise = m(m−1)…(m−e+1) (`kRise_parafermi`), which
       is 0 for e > m (`fall_vanish`): at most m copies of each element.
  4. Necklaces. 1/(1 − kz) = ∏_N (1 − z^N)^(−M(k,N)), with M(k,N) the
     aperiodic necklaces on k colours, checked to z^3 for every k
     (`cyclotomic_z2`, `cyclotomic_z3`). Taking the 1/k-th power and putting
     z = p^(−s) on every prime gives Z_k(s) = ∏_N ζ(Ns)^(M(k,N)/k): k = −1 is
     ζ(s)/ζ(2s), k = 1 is ζ(s) (exponents checked in `fermi_exponents`).
-/
import FusedAlgebra.Base

namespace CircleOp
namespace DialSurvey

open Lean.Grind

/-! ## 1. Units change the weight -/

/-- In units λ times bigger, weight k becomes weight k·λ. -/
theorem units_change_weight {R : Type} [CommSemiring R] (k l a b : R) :
    circK k (l * a) (l * b) = l * circK (k * l) a b := by
  unfold circK; grind

/-- With units 1/k the weight becomes 1: every nonzero k is the circle operation in its own units. -/
theorem weight_to_one (k a b : Rat) (hk : k ≠ 0) :
    circK k (a / k) (b / k) = (1 / k) * circ a b := by
  unfold circK circ; grind

/-! ## 2. Per-element dials are odds -/

/-- Odds t ↦ probability t/(1+t); "in A or in B" (independent) has odds t ◯ s. -/
theorem odds_union (t s : Rat) (ht : 1 + t ≠ 0) (hs : 1 + s ≠ 0) :
    1 - circ t s / (1 + circ t s) = (1 - t / (1 + t)) * (1 - s / (1 + s)) := by
  have hl : 1 + circ t s = (1 + t) * (1 + s) := circ_lift t s
  have hm : circ t s = (1 + t) * (1 + s) - 1 := by rw [← hl]; grind
  rw [hm]; grind

/-- The complement of a subset with odds t has odds 1/t. -/
theorem odds_complement (t : Rat) (ht : t ≠ 0) (ht1 : 1 + t ≠ 0) :
    1 - t / (1 + t) = (1 / t) / (1 + 1 / t) := by
  grind

/-- The uniform dial 2^a − 1 is the exponent a read through the circle-operation ruler:
    dials combine by the circle operation exactly when exponents add. -/
theorem exponent_ruler (a b : Nat) :
    circ ((2 : Int) ^ a - 1) (2 ^ b - 1) = 2 ^ (a + b) - 1 := by
  unfold circ; rw [Int.pow_add]; grind

/-! ## 3. k acting on per-element dials: how many copies an element may have -/

/-- kRise k e = (1)(1 + k)(1 + 2k)…(1 + (e−1)k). -/
def kRise (k : Rat) : Nat → Rat
  | 0 => 1
  | e + 1 => kRise k e * (1 + (e : Rat) * k)

theorem kRise_succ (k : Rat) (e : Nat) : kRise k (e + 1) = kRise k e * (1 + (e : Rat) * k) := rfl

/-- k = −1: an element may appear at most once (subsets). -/
theorem kRise_fermi (e : Nat) : kRise (-1) (e + 2) = 0 := by
  induction e with
  | zero => simp only [kRise]; grind
  | succ n ih => rw [kRise_succ, ih]; grind

/-- k = 0: kRise is 1, so the coefficient is 1/e! (the exponential, Poisson). -/
theorem kRise_poisson (e : Nat) : kRise 0 e = 1 := by
  induction e with
  | zero => rfl
  | succ n ih => rw [kRise_succ, ih]; grind

/-- e!, defined here (core Lean has no factorial). -/
def fct : Nat → Nat
  | 0 => 1
  | n + 1 => (n + 1) * fct n

theorem fct_succ (n : Nat) : fct (n + 1) = (n + 1) * fct n := rfl

/-- k = +1: kRise is e!, so every multiplicity has weight 1 (multisets). -/
theorem kRise_bose (e : Nat) : kRise 1 e = (fct e : Rat) := by
  induction e with
  | zero => rfl
  | succ n ih =>
    rw [kRise_succ, ih]; simp only [fct]
    grind

/-- The falling product m(m−1)…(m−e+1) in Rat. -/
def fall (m : Rat) : Nat → Rat
  | 0 => 1
  | e + 1 => fall m e * (m - (e : Rat))

/-- k = −1/m: m^e · kRise = m(m−1)…(m−e+1). -/
theorem kRise_parafermi (m : Rat) (hm : m ≠ 0) (e : Nat) :
    m ^ e * kRise (-1 / m) e = fall m e := by
  induction e with
  | zero => simp only [kRise, fall, Rat.pow_zero]; grind
  | succ n ih =>
    rw [kRise_succ, Rat.pow_succ]
    simp only [fall]
    rw [← ih]
    have h : m * (1 + (n : Rat) * (-1 / m)) = m - n := by
      rw [Rat.div_def]; have := Rat.mul_inv_cancel _ hm; grind
    grind

/-- At m (a whole number) the falling product vanishes past e = m:
    with k = −1/m no element appears more than m times. -/
theorem fall_vanish (m e : Nat) : fall (m : Rat) (m + 1 + e) = 0 := by
  induction e with
  | zero =>
    show fall (m : Rat) m * ((m : Rat) - (m : Rat)) = 0
    grind
  | succ n ih =>
    show fall (m : Rat) (m + 1 + n) * ((m : Rat) - ((m + 1 + n : Nat) : Rat)) = 0
    rw [ih]; grind

/-- At e = m it does not vanish: m copies are allowed. -/
theorem fall_full (m : Nat) : fall (m : Rat) m = (fct m : Rat) := by
  -- fall m j · (m − j)! = m!, for every j ≤ m, written with m = r + j
  have key : ∀ j r : Nat, fall ((r + j : Nat) : Rat) j * (fct r : Rat) = (fct (r + j) : Rat) := by
    intro j
    induction j with
    | zero => intro r; simp only [fall, Nat.add_zero]; grind
    | succ n ih =>
      intro r
      have h := ih (r + 1)
      have e : r + 1 + n = r + (n + 1) := by omega
      rw [e, fct_succ r] at h
      simp only [fall]
      rw [Rat.natCast_mul] at h
      have c1 : ((r + 1 : Nat) : Rat) = (r : Rat) + 1 := by rw [Rat.natCast_add]; simp
      have c2 : ((r + (n + 1) : Nat) : Rat) = (r : Rat) + (n : Rat) + 1 := by
        rw [Rat.natCast_add, Rat.natCast_add]; simp; grind
      rw [c1] at h; rw [c2] at h ⊢
      generalize fall ((r : Rat) + (n : Rat) + 1) n = F at h ⊢
      generalize ((fct r : Nat) : Rat) = G at h ⊢
      rw [← h]; grind
  have := key m 0
  simp only [fct, Nat.zero_add] at this
  grind

/-! ## 4. Necklaces: the cyclotomic identity, low orders, for every k -/

/-- Aperiodic necklaces of length 1, 2, 3 on k colours. -/
def M1 (k : Rat) : Rat := k
def M2 (k : Rat) : Rat := (k * k - k) / 2
def M3 (k : Rat) : Rat := (k * k * k - k) / 3

/-- z² coefficient: (1 − z)^(−M1) gives M1(M1+1)/2, (1 − z²)^(−M2) gives M2;
    together k². Every two-letter word splits uniquely into Lyndon words. -/
theorem cyclotomic_z2 (k : Rat) : M1 k * (M1 k + 1) / 2 + M2 k = k * k := by
  unfold M1 M2; grind

/-- z³ coefficient: k(k+1)(k+2)/6 + k·M2 + M3 = k³. -/
theorem cyclotomic_z3 (k : Rat) :
    M1 k * (M1 k + 1) * (M1 k + 2) / 6 + M1 k * M2 k + M3 k = k * k * k := by
  unfold M1 M2 M3; grind

/-- k = −1: exponents M/k on ζ(Ns) are 1, −1, 0: Z = ζ(s)/ζ(2s).
    k = +1: exponents 1, 0, 0: Z = ζ(s). -/
theorem fermi_exponents :
    M1 (-1) / (-1) = 1 ∧ M2 (-1) / (-1) = -1 ∧ M3 (-1) / (-1) = 0 := by
  unfold M1 M2 M3; refine ⟨?_, ?_, ?_⟩ <;> grind

theorem bose_exponents : M1 1 = 1 ∧ M2 1 = 0 ∧ M3 1 = 0 := by
  unfold M1 M2 M3; refine ⟨?_, ?_, ?_⟩ <;> grind

/-- For whole k the necklace counts are whole: k² − k is even, k³ − k divisible by 3. -/
theorem necklace_whole (k : Int) : 2 ∣ k * k - k ∧ 3 ∣ k * k * k - k := by
  constructor
  · have hk : k = 2 * (k / 2) + k % 2 := by omega
    have hr : k % 2 = 0 ∨ k % 2 = 1 := by omega
    generalize k / 2 = q at hk
    rcases hr with h | h <;> rw [h] at hk <;> subst hk
    · exact ⟨2 * q * q - q, by grind⟩
    · exact ⟨2 * q * q + q, by grind⟩
  · have hk : k = 3 * (k / 3) + k % 3 := by omega
    have hr : k % 3 = 0 ∨ k % 3 = 1 ∨ k % 3 = 2 := by omega
    generalize k / 3 = q at hk
    rcases hr with h | h | h <;> rw [h] at hk <;> subst hk
    · exact ⟨9 * q * q * q - q, by grind⟩
    · exact ⟨9 * q * q * q + 9 * q * q + 2 * q, by grind⟩
    · exact ⟨9 * q * q * q + 18 * q * q + 11 * q + 2, by grind⟩

end DialSurvey
end CircleOp
