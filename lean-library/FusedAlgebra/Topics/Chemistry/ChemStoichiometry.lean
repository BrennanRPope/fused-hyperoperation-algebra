/-
  Chemistry 5, stoichiometry, definite and multiple proportions
  =============================================================

  Companion page and checks: ../index.html, ../checks.py (14 checks).
  Core Lean 4.37.0-rc1, no mathlib, no `sorry`.

  Yields multiply, so losses combine by the weighted mult-add at k = −1
  (the library's `circK (-1)`, failure-union). n equal steps lose
  `wpow (-1) ℓ n = 1 − (1 − ℓ)^n`. A synthesis tree's yield along any path
  depends only on how many couplings that path passes through.

  Isotope patterns are generating functions: a molecule's pattern is the
  product of its atoms' polynomials, which is the mult-add lift.

  Proportions: combining two compounds adds their atom counts, so the
  combined formula is the mediant of the parents' ratios and lies between.

  LIBRARY NOTE: the local infix for circ was replaced by Base's scoped ◯ (`open scoped CircleOp`) in the 1.0 rename.
  Merged at Topics/Chemistry/ChemStoichiometry; it now uses Base's ◯. Two exact twins removed
  (see the LIBRARY comments).
-/
import FusedAlgebra.Base

namespace ChemStoichiometry
open CircleOp

open scoped CircleOp

/-- Two steps: overall loss is the k = −1 mult-add of the step losses. -/
theorem loss_two (l1 l2 : Rat) : 1 - (1 - l1) * (1 - l2) = circK (-1) l1 l2 := by
  rw [circK_either]

/-- n equal steps: overall loss 1 − (1 − ℓ)^n is the k = −1 power. -/
theorem loss_n (l : Rat) (n : Nat) : 1 - (1 - l) ^ n = wpow (-1) l n := by
  rw [wpow_neg_one]

-- LIBRARY: yield_lift (The lift: yields are 1 − loss and multiply) was an exact twin of StatFused.de_morgan (also InfoFused.survival) and is removed.

/-- Losses compose in any grouping (associativity of the fold). -/
theorem loss_assoc (a b c : Rat) : circK (-1) (circK (-1) a b) c = circK (-1) a (circK (-1) b c) :=
  circK_assoc (-1) a b c

/-- Ten steps at 90%: (9/10)^10 overall, so 1 − (9/10)^10 is lost. -/
theorem ten_steps : wpow (-1) (1/10 : Rat) 10 = 1 - (9/10) ^ 10 := by
  rw [← loss_n]; grind

/-- Depth of a balanced coupling tree on 2^d fragments is d; a linear chain on n fragments has depth n − 1. -/
def treeDepth (d : Nat) : Nat := d
def chainDepth (n : Nat) : Nat := n - 1
/-- Eight fragments: three levels convergent, seven linear. -/
theorem eight_fragments : treeDepth 3 = 3 ∧ chainDepth 8 = 7 := ⟨rfl, rfl⟩

/-- A two-isotope atom's polynomial p + q x; two atoms multiply to p² + 2pq x + q² x². -/
theorem two_atoms (p q x : Rat) : (p + q * x) * (p + q * x) = p * p + 2 * p * q * x + q * q * x * x := by
  grind

/-- Chlorine at 3 : 1: Cl₂ is 9 : 6 : 1. -/
theorem chlorine_961 (x : Rat) : (3 + x) * (3 + x) = 9 + 6 * x + x * x := by grind

/-- Bromine at 1 : 1: Br₂ is 1 : 2 : 1. -/
theorem bromine_121 (x : Rat) : (1 + x) * (1 + x) = 1 + 2 * x + x * x := by grind

-- LIBRARY: any_heavy (At least one heavy atom among n: 1 − (1 − p)^n, the k = −1 power of p) was an exact twin of loss_n, earlier in this file and is removed.

/-- Mediant: combining a:b with c:d gives (a+c):(b+d), which lies between them. -/
theorem mediant_between (a b c d : Nat) (h : a * d < c * b) :
    a * (b + d) < (a + c) * b ∧ (a + c) * d < c * (b + d) := by
  constructor
  · rw [Nat.mul_add, Nat.add_mul]; omega
  · rw [Nat.add_mul, Nat.mul_add]; omega

/-- Fe₃O₄ = FeO + Fe₂O₃: Fe:O 3:4 from 1:1 and 2:3. -/
theorem magnetite : (1 + 2, 1 + 3) = (3, 4) := rfl
/-- 2/3 < 3/4 < 1/1 for Fe:O, by cross-multiplication. -/
theorem magnetite_between : 2 * 4 < 3 * 3 ∧ 3 * 1 < 1 * 4 := by decide

/-- Dalton: oxygen atoms per two nitrogens in N₂O, NO, N₂O₃, NO₂, N₂O₅ are 1, 2, 3, 4, 5. -/
def oxygenPerTwoN : List (Nat × Nat) := [(2, 1), (1, 1), (2, 3), (1, 2), (2, 5)]
theorem multiple_proportions :
    oxygenPerTwoN.map (fun p => p.2 * 2 / p.1) = [1, 2, 3, 4, 5] := by decide

/-- Haber balances: N 2 = 2·1, H 3·2 = 2·3. -/
theorem haber_balanced : 1 * 2 = 2 * 1 ∧ 3 * 2 = 2 * 3 := by decide

/-- Respiration balances C, H, O: C₆H₁₂O₆ + 6 O₂ → 6 CO₂ + 6 H₂O. -/
theorem respiration_balanced : 6 = 6 * 1 ∧ 12 = 6 * 2 ∧ 6 + 6 * 2 = 6 * 2 + 6 * 1 := by decide

/-- Nitrogenase's N₂ + 8 H → 2 NH₃ + H₂ balances, and so does N₂ + 6 H → 2 NH₃:
    stoichiometry alone leaves the H₂ count free. -/
theorem nitrogenase_family : (2 = 2 * 1 ∧ 8 = 2 * 3 + 2) ∧ (2 = 2 * 1 ∧ 6 = 2 * 3 + 0) := by decide

end ChemStoichiometry
