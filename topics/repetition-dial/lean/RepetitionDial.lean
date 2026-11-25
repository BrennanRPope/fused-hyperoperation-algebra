/-
  The repetition dial: how many times one element may be used
  ===========================================================

  Part of Fused Hyperoperation Algebra (1.0 checklist, item 1). Companion page
  and checks: math-theory/repetition-dial/ (README.md, checks.py, figures).
  Core Lean 4.37 only (no mathlib), no `sorry`. Imports Base only.

  An element ("atom") with twist k has the weight Z_k(x) = (1 + kx)^(1/k):
  k = 1 is "at most once" (1 + x, subsets), k = 0 is e^x (classical),
  k = −1 is "any number of times" (1/(1 − x), multisets).

  1. Atoms compose by the weighted mult-add: if 1 + ka = A and 1 + kb = B then
     1 + k(a ⊞ₖ b) = AB (`atoms_compose`). Read through the k-log (Z^k − 1)/k,
     every atom's value is its own weight x, for every k.
  2. The coefficient of xⁿ in Z_k is cellCoeff k n = ∏_{j<n}(1 − jk)/n!, a
     polynomial in k, so the dial is smooth straight through k = 0
     (`gbinom_cells`: it equals the binomial C(1/k, n)·kⁿ for k ≠ 0;
     `cellCoeff_zero_succ`: at k = 0 it is 1/n!).
  3. 1/k is the number of cells per element. k = 1/m caps every element at m
     copies (`cap_at_m`), and the count is exactly "choose n of m cells"
     (`cells_exclusive`, with `excl` the number of ways). k ≤ 0 never caps:
     every coefficient is positive (`unlimited_when_nonpos`), and bosonic cells
     always have a way to hold any number of copies (`bose_pos`). A positive k
     that is not 1/m gives a negative count (`two_thirds_negative`).
  4. Mean copies n = x/(1 + kx) and fluctuation n(1 − kn) = x/(1 + kx)²
     (`fluctuation`), so 1/n − 1/x = k (`inverse_gap`): k is read off any
     single measurement.
  5. Haldane's semions (g = 1/2: w(1 + w) = 1/x², n = 1/(w + 1/2)) have no
     single k: a constant 1/n − 1/x = k would hold at most at one density
     (`semion_one_density`, `semion_no_constant_k`).
  The capped count 1 + x + … + x^m = (1 − x^(m+1))/(1 − x) is
  SizeDialResidual.at_most_m; 1/(1 − x) as "at most once" at doubled sizes is
  SizeDialResidual.distinct_times_doubled.
-/
import FusedAlgebra.Base

namespace RepetitionDial

open MultAdd
open Lean.Grind

set_option linter.unusedSectionVars false

/-! ## 1. Atoms compose by the weighted mult-add -/

section Compose
variable {R : Type} [Lean.Grind.CommRing R]

/-- **Atoms compose by ⊞ₖ.** If the k-powers of two atoms are 1 + ka and
    1 + kb, their product is 1 + k(a ⊞ₖ b). -/
theorem atoms_compose (k a b A B : R) (hA : 1 + k * a = A) (hB : 1 + k * b = B) :
    1 + k * wAdd k a b = A * B := by
  rw [← hA, ← hB]; exact wAdd_lift k a b

/-- Three atoms: the order of composing does not matter. -/
theorem atoms_compose3 (k a b c : R) :
    1 + k * wAdd k (wAdd k a b) c = (1 + k * a) * (1 + k * b) * (1 + k * c) := by
  rw [wAdd_lift, wAdd_lift]

end Compose

/-! ## 2. The coefficients are polynomials in k -/

/-- Coefficient of xⁿ in (1 + kx)^(1/k): ∏_{j<n}(1 − jk)/n!. -/
def cellCoeff (k : Rat) : Nat → Rat
  | 0 => 1
  | n + 1 => cellCoeff k n * (1 - (n : Rat) * k) / ((n : Rat) + 1)

/-- Generalised binomial coefficient C(a, n). -/
def gbinom (a : Rat) : Nat → Rat
  | 0 => 1
  | n + 1 => gbinom a n * (a - (n : Rat)) / ((n : Rat) + 1)

theorem natCast_succ_pos (n : Nat) : (0 : Rat) < (n : Rat) + 1 := by
  have : (0 : Rat) ≤ (n : Rat) := by exact_mod_cast Nat.zero_le n
  grind

/-- **The binomial series is the cell polynomial.** For k ≠ 0,
    C(1/k, n)·kⁿ = ∏_{j<n}(1 − jk)/n!. -/
theorem gbinom_cells (k : Rat) (hk : k ≠ 0) (n : Nat) :
    gbinom (1 / k) n * k ^ n = cellCoeff k n := by
  induction n with
  | zero => simp [gbinom, cellCoeff]
  | succ n ih =>
    have hn : (n : Rat) + 1 ≠ 0 := by have := natCast_succ_pos n; grind
    simp only [gbinom, cellCoeff]
    rw [← ih, Rat.pow_succ]
    grind

/-- At k = 0 each step divides by n + 1: the coefficients are 1/n! (e^x). -/
theorem cellCoeff_zero_succ (n : Nat) :
    cellCoeff 0 (n + 1) = cellCoeff 0 n / ((n : Rat) + 1) := by
  simp only [cellCoeff]; grind

/-- k = 1: one exclusive cell, so no element twice (coefficient of x² is 0). -/
theorem cellCoeff_one_two : cellCoeff 1 2 = 0 := by decide +kernel

/-- k = −1: one bosonic cell, every number of copies counted once. -/
theorem cellCoeff_neg_one (n : Nat) : cellCoeff (-1) n = 1 := by
  induction n with
  | zero => rfl
  | succ n ih =>
    have hn : (n : Rat) + 1 ≠ 0 := by have := natCast_succ_pos n; grind
    simp only [cellCoeff, ih]
    grind

/-! ## 3. 1/k is the number of cells per element -/

/-- **k = 1/m caps each element at m copies.** -/
theorem cap_at_m (m : Nat) (hm : 0 < m) (n : Nat) (hn : m < n) :
    cellCoeff (1 / (m : Rat)) n = 0 := by
  induction n with
  | zero => omega
  | succ n ih =>
    simp only [cellCoeff]
    by_cases h : m = n
    · subst h
      have hm' : (m : Rat) ≠ 0 := by
        have : (0 : Rat) < (m : Rat) := by exact_mod_cast hm
        grind
      have : (1 : Rat) - (m : Rat) * (1 / (m : Rat)) = 0 := by grind
      rw [this]; grind
    · have : m < n := by omega
      rw [ih this]; grind

/-- Ways to put n copies in m exclusive cells (each cell at most one). -/
def excl : Nat → Nat → Nat
  | 0, 0 => 1
  | 0, _ + 1 => 0
  | _ + 1, 0 => 1
  | m + 1, n + 1 => excl m (n + 1) + excl m n

/-- Ways to put n copies in m bosonic cells (each cell any number). -/
def bose : Nat → Nat → Nat
  | 0, 0 => 1
  | 0, _ + 1 => 0
  | _ + 1, 0 => 1
  | m + 1, n + 1 => bose m (n + 1) + bose (m + 1) n

theorem excl_zero_right (m : Nat) : excl m 0 = 1 := by cases m <;> rfl

theorem excl_over (m n : Nat) (h : m < n) : excl m n = 0 := by
  induction m generalizing n with
  | zero => cases n with
    | zero => omega
    | succ n => rfl
  | succ m ih => cases n with
    | zero => omega
    | succ n => simp only [excl]; rw [ih (n + 1) (by omega), ih n (by omega)]

/-- **Bosonic cells never cap:** any number of copies has a way in. -/
theorem bose_pos (m n : Nat) : 0 < bose (m + 1) n := by
  induction n with
  | zero => cases m <;> simp [bose]
  | succ n ih => simp only [bose]; omega

/-- One bosonic cell: exactly one way for every n (that is 1/(1 − x)). -/
theorem bose_zero_succ (n : Nat) : bose 0 (n + 1) = 0 := by simp [bose]

theorem bose_one (n : Nat) : bose 1 n = 1 := by
  induction n with
  | zero => simp [bose]
  | succ n ih => simp [bose, ih]

/-- One exclusive cell: one way for 0 or 1 copy, none for 2 or more. -/
theorem excl_one_two (n : Nat) : excl 1 (n + 2) = 0 := excl_over 1 (n + 2) (by omega)

/-- Absorption: (n + 1)·C(m, n + 1) = (m − n)·C(m, n). -/
theorem excl_absorb (m n : Nat) :
    ((n : Int) + 1) * (excl m (n + 1) : Int) = ((m : Int) - n) * (excl m n : Int) := by
  induction m generalizing n with
  | zero =>
    cases n with
    | zero => simp [excl]
    | succ n => simp [excl]
  | succ m ih =>
    cases n with
    | zero =>
      have h1 := ih 0
      simp only [excl, excl_zero_right] at h1 ⊢
      push_cast at h1 ⊢
      grind
    | succ n =>
      have h1 := ih (n + 1)
      have h2 := ih n
      simp only [excl] at h1 h2 ⊢
      push_cast at h1 h2 ⊢
      grind

/-- **The k = 1/m coefficients are cell counts:** cellCoeff (1/m) n · mⁿ = C(m, n),
    the number of ways to fill n of m exclusive cells. -/
theorem cells_exclusive (m : Nat) (hm : 0 < m) (n : Nat) :
    cellCoeff (1 / (m : Rat)) n * (m : Rat) ^ n = (excl m n : Rat) := by
  have hm' : (m : Rat) ≠ 0 := by
    have : (0 : Rat) < (m : Rat) := by exact_mod_cast hm
    grind
  induction n with
  | zero => simp [cellCoeff, excl_zero_right]
  | succ n ih =>
    have hn : (n : Rat) + 1 ≠ 0 := by have := natCast_succ_pos n; grind
    have hab := excl_absorb m n
    have habQ : ((n : Rat) + 1) * (excl m (n + 1) : Rat) = ((m : Rat) - n) * (excl m n : Rat) := by
      have := congrArg (fun z : Int => (z : Rat)) hab
      push_cast at this; exact this
    simp only [cellCoeff]
    rw [Rat.pow_succ]
    have hx : (excl m (n + 1) : Rat) = ((m : Rat) - n) * (excl m n : Rat) / ((n : Rat) + 1) := by
      grind
    rw [hx, ← ih]
    grind

/-- **k ≤ 0 never caps:** every coefficient is positive. -/
theorem unlimited_when_nonpos (k : Rat) (hk : k ≤ 0) : ∀ n, 0 < cellCoeff k n := by
  intro n
  induction n with
  | zero => simp [cellCoeff]; decide
  | succ n ih =>
    simp only [cellCoeff]
    have h0 : (0 : Rat) ≤ (n : Rat) := by exact_mod_cast Nat.zero_le n
    have hnk : 0 ≤ (n : Rat) * (-k) := Rat.mul_nonneg h0 (by grind)
    have hf : 0 < 1 - (n : Rat) * k := by grind
    have hinv : 0 < ((n : Rat) + 1)⁻¹ := Rat.inv_pos.mpr (natCast_succ_pos n)
    have := Rat.mul_pos (Rat.mul_pos ih hf) hinv
    rw [Rat.div_def]; exact this

/-- **A positive k that is not 1/m gives a negative count:** at k = 2/3 the
    coefficient of x³ is −1/54. -/
theorem two_thirds_negative : cellCoeff (2 / 3) 3 = -1 / 54 := by decide +kernel

/-! ## 4. Mean copies and fluctuation -/

section Occupation
variable {K : Type} [Lean.Grind.Field K]

/-- **Fluctuation.** With n = x/(1 + kx): n(1 − kn) = x/(1 + kx)². Fermions
    (k = 1) quiet the count, bosons (k = −1) bunch it, k = 0 is Poisson. -/
theorem fluctuation (k x : K) (h : 1 + k * x ≠ 0) :
    (x / (1 + k * x)) * (1 - k * (x / (1 + k * x))) = x / ((1 + k * x) * (1 + k * x)) := by
  grind

/-- **k from one measurement:** 1/n − 1/x = k. -/
theorem inverse_gap (k x : K) (hx : x ≠ 0) (h : 1 + k * x ≠ 0) :
    (x / (1 + k * x))⁻¹ - x⁻¹ = k := by
  grind

end Occupation

/-! ## 5. Haldane's semions have no single k -/

/-- With w(1 + w) = 1/x² (u = 1/x) and n = 1/(w + 1/2), a constant gap
    1/n − 1/x = k means w + 1/2 − k = u, so −2kw + (1/2 − k)² = 0. -/
theorem semion_one_density (k w : Rat) (h : (w + 1/2 - k) ^ 2 = w * (1 + w)) :
    -2 * k * w + (1/2 - k) ^ 2 = 0 := by grind

/-- **No constant k for semions:** two different densities cannot share one gap k. -/
theorem semion_no_constant_k (k w₁ w₂ : Rat) (hne : w₁ ≠ w₂)
    (h₁ : (w₁ + 1/2 - k) ^ 2 = w₁ * (1 + w₁)) (h₂ : (w₂ + 1/2 - k) ^ 2 = w₂ * (1 + w₂)) :
    False := by
  have e₁ := semion_one_density k w₁ h₁
  have e₂ := semion_one_density k w₂ h₂
  have hk : k * (w₁ - w₂) = 0 := by grind
  have hd : w₁ - w₂ ≠ 0 := by grind
  have k0 : k = 0 := by
    rcases Rat.mul_eq_zero.mp hk with h | h
    · exact h
    · exact absurd h hd
  subst k0
  grind

end RepetitionDial
