import FusedAlgebra.Base

/-!
# Riemann, round 9: the Hilbert–Pólya string probe (DRAFT)

Small exact facts behind the Gram/Rosser analysis and the string's settings.
At a Gram point g_m, θ(g_m) = mπ, so N(g_m) = m + 1 + S(g_m) with S an integer there.
-/

namespace RiemannRound9

/-- Zeros in a Gram interval: 1 plus the jump of S. -/
theorem gram_interval_count (m S1 S2 : Int) :
    (m + 1 + 1 + S2) - (m + 1 + S1) = 1 + (S2 - S1) := by omega

/-- Gram's law holds in an interval exactly when S does not jump across it. -/
theorem gram_law_iff (S1 S2 : Int) : 1 + (S2 - S1) = 1 ↔ S2 = S1 := by omega

/-- An empty Gram interval means S falls by one across it. -/
theorem empty_interval (S1 S2 : Int) (h : 1 + (S2 - S1) = 0) : S2 = S1 - 1 := by omega

/-- Zeros in a Gram block of length l: l plus the jump of S. -/
theorem rosser_block (m l S1 S2 : Int) :
    (m + l + 1 + S2) - (m + 1 + S1) = l + (S2 - S1) := by omega

/-- Rosser's rule holds in a block exactly when S takes the same value at its two ends. -/
theorem rosser_iff (l S1 S2 : Int) : l + (S2 - S1) = l ↔ S2 = S1 := by omega

/-- Two intervals, each with at least one zero, and two zeros in all: one each. -/
theorem two_exact (a b : Nat) (ha : 1 ≤ a) (hb : 1 ≤ b) (h : a + b = 2) :
    a = 1 ∧ b = 1 := by omega

/-- The same for three intervals: the counting step behind a "perfect rule". -/
theorem three_exact (a b c : Nat) (ha : 1 ≤ a) (hb : 1 ≤ b) (hc : 1 ≤ c)
    (h : a + b + c = 3) : a = 1 ∧ b = 1 ∧ c = 1 := by omega

/-- Turing's squeeze: zeros found on the line ≤ all zeros, and all ≤ found, so equal. -/
theorem turing_squeeze (on total : Nat) (h1 : on ≤ total) (h2 : total ≤ on) :
    on = total := by omega

/-- A good Gram point (S even) stays good when S moves by two. -/
theorem good_parity (S : Int) (h : S % 2 = 0) : (S + 2) % 2 = 0 := by omega

/-- The Fermi primon gas (k = −1): (1 + u)(1 − u) = 1 − u², i.e. ζ(s)/ζ(2s). -/
theorem fermi_factor (u : Int) : (1 + u) * (1 - u) = 1 - u * u := by grind

/-- Curvature dial: roots a, b of 1 − x + d x² (a + b = 1, ab = d) give (1 − a)(1 − b) = d. -/
theorem curvature_roots (a b d : Int) (h1 : a + b = 1) (h2 : a * b = d) :
    (1 - a) * (1 - b) = d := by grind

/-- Time reversal pairs each orbit with its reverse: the weight doubles (A = 2). -/
theorem time_reversal_double (w : Int) : w + w = 2 * w := by omega

end RiemannRound9
