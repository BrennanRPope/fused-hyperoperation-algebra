/-
  Chemistry 4, bonding and valence
  ================================

  Companion page and checks: ../index.html, ../checks.py (19 checks).
  Core Lean 4.37.0-rc1, no mathlib, no `sorry`.

  Overlap: |a φa + b φb|² = a² + b² + 2ab·S, the same law of cosines as the
  invariant mass (chemistry 2), with the overlap S as the cosine.

  Hückel: the chain recurrence P(n+2) = k·P(n+1) − q·P(n) is the landscape
  point (k, q). Its Cassini invariant P(n+1)² − P(n+2)P(n) is q^(n+1); a
  plain carbon chain has q = 1 and k = x (energy in units of β). For x² < 4
  the characteristic polynomial t² − x t + 1 has no real root (the turning
  country: Bloch waves, the band); x² > 4 is the gap (two walls).

  LIBRARY NOTE: `local infixl " ✳ "` for multAdd; delete it if Base declares ✳.
-/
import FusedAlgebra.Base

namespace ChemBonding
open MultAdd

local infixl:65 " ✳ " => multAdd

/-- Squared length of a φa + b φb when the overlap is S. -/
def norm2 (a b S : Rat) : Rat := a * a + b * b + 2 * a * b * S

/-- S = 1 is the border: identical orbitals, coefficients add. -/
theorem overlap_border (a b : Rat) : norm2 a b 1 = (a + b) * (a + b) := by
  unfold norm2; grind

/-- S = −1: opposite phases, coefficients subtract. -/
theorem overlap_opposite (a b : Rat) : norm2 a b (-1) = (a - b) * (a - b) := by
  unfold norm2; grind

/-- S = 0: orthogonal orbitals, Pythagoras. -/
theorem overlap_orthogonal (a b : Rat) : norm2 a b 0 = a * a + b * b := by
  unfold norm2; grind

/-- The equal-weight combination φa + φb has length² 2 + 2S = (1+S)·2, the lift of S at weight 1. -/
theorem bonding_norm (S : Rat) : norm2 1 1 S = 2 * (1 + S) := by
  unfold norm2; grind

/-- Asymmetric splitting: with bonding (α+β)/(1+S) and antibonding (α−β)/(1−S),
    rise·(1−S) and fall·(1+S) are both αS − β, so rise/fall = (1+S)/(1−S).
    Stated cleared of denominators: the two numerators agree. -/
theorem split_asymmetry (α β S : Rat) :
    (α - β) - α * (1 - S) = α * S - β ∧ α * (1 + S) - (α + β) = α * S - β := by
  constructor <;> grind

/-- 1ⁿ = 1 over Rat (core has no lemma by this name). -/
theorem one_pow' (n : Nat) : (1 : Rat) ^ n = 1 := by
  induction n with
  | zero => rfl
  | succ m ih => rw [Rat.pow_succ, ih]; grind

/-- Squares are nonnegative over Rat. -/
theorem sq_nonneg' (y : Rat) : 0 ≤ y * y := by
  rcases (Rat.le_total : (0 : Rat) ≤ y ∨ y ≤ 0) with h | h
  · exact Rat.mul_nonneg h h
  · have h2 : 0 ≤ -y := by grind
    have := Rat.mul_nonneg h2 h2
    grind

/-- The chain recurrence at landscape point (k, q): P 0 = 1, P 1 = k. -/
def chainP (k q : Rat) : Nat → Rat
  | 0 => 1
  | 1 => k
  | n + 2 => k * chainP k q (n + 1) - q * chainP k q n

/-- Cassini step: the invariant is multiplied by q at each site. -/
theorem cassini_step (k q : Rat) (n : Nat) :
    chainP k q (n + 2) * chainP k q (n + 2) - chainP k q (n + 3) * chainP k q (n + 1)
      = q * (chainP k q (n + 1) * chainP k q (n + 1) - chainP k q (n + 2) * chainP k q n) := by
  simp only [chainP]; grind

/-- Cassini: P(n+1)² − P(n+2)·P(n) = q^(n+1). -/
theorem cassini (k q : Rat) (n : Nat) :
    chainP k q (n + 1) * chainP k q (n + 1) - chainP k q (n + 2) * chainP k q n = q ^ (n + 1) := by
  induction n with
  | zero => simp only [chainP]; grind
  | succ m ih => rw [cassini_step, ih, Rat.pow_succ]; grind

/-- A plain carbon chain (q = 1) has invariant 1 at every length. -/
theorem huckel_unimodular (x : Rat) (n : Nat) :
    chainP x 1 (n + 1) * chainP x 1 (n + 1) - chainP x 1 (n + 2) * chainP x 1 n = 1 := by
  rw [cassini]; exact one_pow' (n + 1)

/-- Butadiene (4 carbons): P 4 = x⁴ − 3x² + 1, whose roots are ±φ, ±1/φ (golden ratio). -/
theorem butadiene (x : Rat) : chainP x 1 4 = x ^ 4 - 3 * x ^ 2 + 1 := by
  simp only [chainP]; grind

/-- Allyl (3 carbons): P 3 = x (x² − 2): levels 0, ±√2. -/
theorem allyl (x : Rat) : chainP x 1 3 = x * (x * x - 2) := by
  simp only [chainP]; grind

/-- Completing the square: 4(t² − x t + 1) = (2t − x)² + (4 − x²). -/
theorem band_square (t x : Rat) : 4 * (t * t - x * t + 1) = (2 * t - x) * (2 * t - x) + (4 - x * x) := by
  grind

/-- Turning country: inside the band (x² < 4) the transfer polynomial has no rational root. -/
theorem no_root_in_band (t x : Rat) (h : x * x < 4) : t * t - x * t + 1 ≠ 0 := by
  intro h0
  have hs := band_square t x
  have hn : 0 ≤ (2 * t - x) * (2 * t - x) := sq_nonneg' _
  have : 0 < 4 - x * x := by grind
  grind

/-- Border: at the band edge x = 2 the root t = 1 is double. -/
theorem band_edge : (1 : Rat) * 1 - 2 * 1 + 1 = 0 := by grind

/-- Gap (two walls): x = 5/2 has the reciprocal real pair t = 2, 1/2. -/
theorem gap_pair : (2 : Rat) * 2 - 5/2 * 2 + 1 = 0 ∧ (1/2 : Rat) * (1/2) - 5/2 * (1/2) + 1 = 0 := by
  constructor <;> grind

/-- Hückel's rule: closing a ring's lowest level plus m degenerate pairs takes 2 + 4m electrons. -/
def huckelClosed (m : Nat) : Nat := 2 + 4 * m
/-- Möbius rings have only pairs: 4m electrons. -/
def mobiusClosed (m : Nat) : Nat := 4 * m

theorem huckel_rule (m : Nat) : huckelClosed m % 4 = 2 := by
  unfold huckelClosed; omega
theorem mobius_rule (m : Nat) : mobiusClosed m % 4 = 0 := by
  unfold mobiusClosed; omega

/-- Benzene: 6 = 4·1 + 2 is a Hückel closed shell. -/
theorem benzene_aromatic : huckelClosed 1 = 6 := rfl

/-- MO bond order (bonding − antibonding)/2 for the second-row diatomics, in electrons. -/
def bondOrder2 (bonding anti : Nat) : Nat := (bonding - anti) / 2
/-- N₂: 8 bonding, 2 antibonding valence electrons → triple bond. -/
theorem n2_triple : bondOrder2 8 2 = 3 := rfl
/-- O₂: 8 bonding, 4 antibonding → double. -/
theorem o2_double : bondOrder2 8 4 = 2 := rfl
/-- F₂: 8 bonding, 6 antibonding → single. -/
theorem f2_single : bondOrder2 8 6 = 1 := rfl

/-- The dial applied three times, (1 + k e)³ = 1 + k·(triple): the associative ladder the
    checks test against bond energies (it fits C–C to 3% and misses N–N by 13%). -/
theorem dial_ladder (k e : Rat) :
    (1 + k * e) * (1 + k * e) * (1 + k * e)
      = 1 + k * (3 * e + 3 * k * e * e + k * k * e * e * e) := by
  grind

end ChemBonding
