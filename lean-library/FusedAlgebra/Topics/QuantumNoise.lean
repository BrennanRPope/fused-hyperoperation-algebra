/-
  Quantum noise through the mult-add cube
  =======================================

  Part of Fused Hyperoperation Algebra. Companion page: index.html.
  Numbers: checks.py. Pictures: figures.py.

  STATUS: core Lean 4 only (no mathlib), no `sorry`. Uses core `Rat` and the
  `grind` tactic for ring arithmetic and `decide` for exact rational checks.
  Run it with `lean QuantumNoise.lean`. Checked with Lean 4.37.0-rc1.

  Mult-add with a weight k:  a ⊕ₖ b = a + b + k·a·b.
  k = 1 is mult-add itself. The noise channels of a qubit use
    k = −1    amplitude damping (T1 decay), "any fault counts"
    k = −2    bit-flip and phase-flip, "only an odd number of faults counts"
    k = −4/3  depolarizing noise on one qubit

  What Lean does NOT check here: that these formulas describe the quantum
  channels (that is the density-matrix simulation in checks.py), anything
  with real numbers, limits or continuity, and the threshold as a number
  (Lean only checks that the threshold lies between two given rationals).

  LIBRARY NOTE: weighted mult-add, its power and its laws come from FusedAlgebra.Base; this module is now wrapped in namespace QuantumNoise.
-/
import FusedAlgebra.Base

namespace QuantumNoise

open CircleOp

/-! ## 1. Every weight is the same mult-add, scaled
  (commutative, associative, identity 0, scaling, the lift and the wall are in Base:
  circK_comm, circK_assoc, circK_zero, circK_scale, circK_lift, circK_wall) -/

/-! ## 2. The three channel rules -/

theorem amplitude_damping (g h : Rat) : 1 - circK (-1) g h = (1 - g) * (1 - h) := by unfold circK; grind
theorem bit_flip (p q : Rat) : 1 - 2 * circK (-2) p q = (1 - 2 * p) * (1 - 2 * q) := by unfold circK; grind
theorem depolarizing (p q : Rat) : 1 - 4/3 * circK (-4/3) p q = (1 - 4/3 * p) * (1 - 4/3 * q) := by unfold circK; grind
/-- The rule from the call, p + q − pq, is the k = −1 rule. It is right for amplitude
    damping. For bit-flips it is off by exactly p·q. -/
theorem call_rule_gap (p q : Rat) : circK (-1) p q - circK (-2) p q = p * q := by unfold circK; grind

/-! ## 3. Stacking n channels climbs the cube -/

/-- n stacked bit-flips: the chance of an odd number of flips P satisfies 1 − 2P = (1 − 2p)ⁿ. -/
theorem stacked_bit_flips (p : Rat) (n : Nat) : 1 - 2 * wpow (-2) p n = (1 - 2 * p) ^ n := by
  have := wpow_lift (-2) p n; grind

/-! ## 4. Pauli channels: the tetrahedron in the cube -/

/-- A Pauli channel: chances of doing nothing, X, Y or Z. -/
structure Pauli where
  i : Rat
  x : Rat
  y : Rat
  z : Rat

/-- Its point in the cube [−1,1]³ (how much of each Bloch direction survives). -/
def Pauli.lx (p : Pauli) : Rat := p.i + p.x - p.y - p.z
def Pauli.ly (p : Pauli) : Rat := p.i - p.x + p.y - p.z
def Pauli.lz (p : Pauli) : Rat := p.i - p.x - p.y + p.z

/-- Doing one channel then another: X·Y = Z (up to phase) and so on. -/
def Pauli.comp (p q : Pauli) : Pauli :=
  { i := p.i*q.i + p.x*q.x + p.y*q.y + p.z*q.z
    x := p.i*q.x + p.x*q.i + p.y*q.z + p.z*q.y
    y := p.i*q.y + p.y*q.i + p.x*q.z + p.z*q.x
    z := p.i*q.z + p.z*q.i + p.x*q.y + p.y*q.x }

theorem comp_lx (p q : Pauli) : (p.comp q).lx = p.lx * q.lx := by
  simp only [Pauli.comp, Pauli.lx]; grind
theorem comp_ly (p q : Pauli) : (p.comp q).ly = p.ly * q.ly := by
  simp only [Pauli.comp, Pauli.ly]; grind
theorem comp_lz (p q : Pauli) : (p.comp q).lz = p.lz * q.lz := by
  simp only [Pauli.comp, Pauli.lz]; grind
theorem comp_total (p q : Pauli) :
    (p.comp q).i + (p.comp q).x + (p.comp q).y + (p.comp q).z = (p.i+p.x+p.y+p.z) * (q.i+q.x+q.y+q.z) := by
  simp only [Pauli.comp]; grind

/-- The four faces of the tetrahedron are the four chances (times 4).
    So "inside the tetrahedron" is exactly "all chances are nonnegative". -/
theorem faces (p : Pauli) (h : p.i + p.x + p.y + p.z = 1) :
    1 + p.lx + p.ly + p.lz = 4 * p.i ∧ 1 + p.lx - p.ly - p.lz = 4 * p.x ∧
    1 - p.lx + p.ly - p.lz = 4 * p.y ∧ 1 - p.lx - p.ly + p.lz = 4 * p.z := by
  simp only [Pauli.lx, Pauli.ly, Pauli.lz]; refine ⟨?_, ?_, ?_, ?_⟩ <;> grind

/-! ## 5. Codes: the 3-qubit repetition code -/

/-- Majority vote of 3 fails when 2 or 3 copies flip (the top two layers of the 3-cube). -/
def maj3 (p : Rat) : Rat := 3 * p^2 * (1 - p) + p^3
/-- The vote in the polarization coordinate λ = 1 − 2p. -/
def f (l : Rat) : Rat := (3 * l - l^3) / 2

theorem maj3_polar (p : Rat) : 1 - 2 * maj3 p = f (1 - 2 * p) := by unfold maj3 f; grind
/-- Phase-flips go straight through the bit-flip code: an odd number of them is a
    logical phase-flip, so the polarization is cubed. -/
theorem odd3_polar (q : Rat) : 1 - 2 * (3 * q * (1 - q)^2 + q^3) = (1 - 2 * q)^3 := by grind

/-- The vote only has fixed points at −1, 0 and 1 (the clean corner, the wall, the mirror). -/
theorem f_fixed (l : Rat) : f l - l = l * (1 - l) * (1 + l) / 2 := by unfold f; grind

/-- Between the wall and the clean corner the vote always helps. -/
theorem f_helps (l : Rat) (h0 : 0 < l) (h1 : l < 1) : l < f l := by
  have a : 0 < 1 - l := by grind
  have b : 0 < 1 + l := by grind
  have c : 0 < l * (1 - l) := Rat.mul_pos h0 a
  have d : 0 < l * (1 - l) * (1 + l) := Rat.mul_pos c b
  have e := f_fixed l
  grind

/-! ## 6. The Shor code's thresholds, bracketed exactly
  Two levels of the concatenated Shor code send the bit-flip polarization x to
  f(f(x)⁹) and the phase-flip polarization z to f(f(z³))³. A point that moves up
  is inside the safe box; one that moves down is outside. -/

def gx (x : Rat) : Rat := f (f x ^ 9)
def gz (z : Rat) : Rat := f (f (z ^ 3)) ^ 3

theorem gx_below : gx (806/1000) < 806/1000 := by decide +kernel
theorem gx_above : gx (807/1000) > 807/1000 := by decide +kernel
theorem gz_below : gz (850/1000) < 850/1000 := by decide +kernel
theorem gz_above : gz (851/1000) > 851/1000 := by decide +kernel
-- In flip chances: bit-flips cross between 9.65% and 9.70%, phase-flips between 7.45% and 7.50%.

/-! ## 7. MacWilliams: switching to the dual code is flipping the cube's signs -/

theorem macwilliams_rep3 (x y : Rat) : ((x + y)^3 + (x - y)^3) / 2 = x^3 + 3 * x * y^2 := by grind

/-- The chance a 3-repetition code misses an error (all three flip), computed two ways:
    on the code's own cube, and from the dual code's layers in the polarization 1 − 2p. -/
theorem undetected_two_ways (p : Rat) :
    ((1 - p)^3 + p^3) - (1 - p)^3 = (2/8) * (1 + 3 * (1 - 2*p)^2) - (1 - p)^3 := by grind

/-! ## 8. The echo: composing with the mult-add inverse -/

theorem circ_inverse (a : Rat) (h : 1 + a ≠ 0) : circ a (-a / (1 + a)) = 0 :=
  circ_circInv a h

end QuantumNoise

