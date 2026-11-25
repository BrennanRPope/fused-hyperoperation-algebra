/-
  Rule breakers: which four-wave weight arrays no line of points produces
  ========================================================================
  Companion checks: ../waves/ (breakers.py, breakers2.py, breakers_graphs.json).
  Core Lean 4 only (no mathlib), no `sorry`. Imports FusedAlgebra.Topics.SeveralWaves
  (for `rule4` and `crossR`). Intended path: FusedAlgebra/Topics/RuleBreakers.lean.
  Not compiled in the session that wrote it; every statement was checked with sympy.

  `rule4` (SeveralWaves) is the one polynomial rule that the six pair weights of four
  point-waves obey. Its zero set is exactly the closure of the arrays that come from points,
  so an array that obeys it is a point-wave array or a limit of them, and one that breaks it
  is not even a limit.

  What the file proves:
    A. Graph-shaped arrays (weight on the pairs of a graph, 1 on the others). Trees, the
       triangle with a lone or pendant wave, and every smaller graph obey the rule for all
       weights. The 4-cycle gives rule4 = ∏ (w − 1)²: it breaks the rule unless one of its
       weights is 1. The 4-cycle with a diagonal breaks it unless k = 0.
    B. One weight k on every pair (a single mult-add weight shared by all waves):
       rule4 k k k k k k = −k³ (k − 1)⁸ (k + 2), so four waves carry it only at k = −2
       (besides the trivial 0 and 1). Three waves carry any k: the third wave sits at a root
       of t² − (k + 1)t + 1 (the identity below). Four waves force the two roots to be the
       cube roots of unity, which is k = −2; a fifth wave has no root left.
    C. The k = −2 configuration exists: waves (0, ∞), (1, −1/2), (ω, −ω/2), (ω², −ω²/2)
       with ω a cube root of unity, checked in ℚ(ω). On the sphere these are the four axes of
       a regular tetrahedron, each vertex paired with its antipode.
-/
import FusedAlgebra.Topics.SeveralWaves

namespace RuleBreakers
open SeveralWaves

set_option maxRecDepth 100000

/-! ## A. Graph-shaped arrays (argument order: 01 02 03 12 13 23) -/

/-- a path of three pairs obeys the rule for all weights -/
theorem path_obeys (x y z : Rat) : rule4 x 1 1 y 1 z = 0 := by
  unfold rule4; grind

/-- a star obeys the rule for all weights -/
theorem star_obeys (x y z : Rat) : rule4 x y z 1 1 1 = 0 := by
  unfold rule4; grind

/-- a triangle with a lone fourth wave obeys the rule -/
theorem triangle_obeys (x y z : Rat) : rule4 x y 1 z 1 1 = 0 := by
  unfold rule4; grind

/-- a triangle with a pendant wave obeys the rule for all weights -/
theorem paw_obeys (x y z w : Rat) : rule4 x y 1 z 1 w = 0 := by
  unfold rule4; grind

/-- the 4-cycle 0–1–2–3–0 with diagonals 1: the rule is the product of (w − 1)² over its pairs -/
theorem cycle_rule (x y z w : Rat) :
    rule4 x 1 w y 1 z = (x - 1) * (x - 1) * ((y - 1) * (y - 1)) * ((z - 1) * (z - 1))
      * ((w - 1) * (w - 1)) := by
  unfold rule4; grind

/-- the 4-cycle with one diagonal, one weight k -/
theorem diamond_rule (k : Rat) :
    rule4 k k k k k 1 = k * ((k - 1) * (k - 1) * (k - 1) * (k - 1)) * ((k - 1) * (k - 1) * (k - 1) * (k - 1)) := by
  unfold rule4; grind

/-- the cycle breaks the rule when its weights avoid 1: an example with weight 2 -/
theorem cycle_two_breaks : rule4 2 1 2 2 1 2 = 1 := by
  unfold rule4; decide +kernel

/-! ## B. One weight on every pair -/

/-- every pair of four waves at weight k -/
theorem uniform_rule (k : Rat) :
    rule4 k k k k k k = -(k * k * k) * ((k - 1) * (k - 1) * (k - 1) * (k - 1))
      * ((k - 1) * (k - 1) * (k - 1) * (k - 1)) * (k + 2) := by
  unfold rule4; grind

/-- so k = 3 (or any k other than 0, 1, −2) is not carried by four point-waves -/
theorem uniform_three_breaks : rule4 3 3 3 3 3 3 = -34560 := by
  unfold rule4; decide +kernel

/-- … and k = −2 is -/
theorem uniform_minus_two_obeys : rule4 (-2) (-2) (-2) (-2) (-2) (-2) = 0 := by
  unfold rule4; decide +kernel

/-- with waves (0, ∞) and (1, 1/k), a third wave (t, t/k) has weight k to both exactly when
    t is a root of t² − (k + 1)t + 1 (numerator of the cross ratio minus k) -/
theorem third_wave (k t : Rat) :
    (1 - t) * (k - k * t) - k * ((k - t) * (1 - k * t)) = -(k * (k - 1)) * (t * t - (k + 1) * t + 1) := by
  grind

/-- four waves need the pair of roots to be the cube roots of unity (t² + t + 1 = 0);
    against the quadratic that leaves (k + 2)·t = 0, so k = −2 -/
theorem fourth_wave_forces (k t : Rat) :
    (t * t - (k + 1) * t + 1) - (t * t + t + 1) = -((k + 2) * t) := by
  grind

/-! ## C. The k = −2 configuration in ℚ(ω), ω² = −1 − ω -/

/-- a + b·ω with ω a primitive cube root of unity -/
structure Eis where
  re : Rat
  om : Rat
  deriving DecidableEq

namespace Eis
def mul (x y : Eis) : Eis := ⟨x.re * y.re - x.om * y.om, x.re * y.om + x.om * y.re - x.om * y.om⟩
def sub (x y : Eis) : Eis := ⟨x.re - y.re, x.om - y.om⟩
def add (x y : Eis) : Eis := ⟨x.re + y.re, x.om + y.om⟩
def two : Eis := ⟨2, 0⟩
def zero : Eis := ⟨0, 0⟩
/-- numerator of the cross ratio plus 2 × denominator: zero means weight −2 -/
def minusTwo (ai bi aj bj : Eis) : Eis :=
  add (mul (sub ai aj) (sub bi bj)) (mul two (mul (sub ai bj) (sub bi aj)))
def den (ai bi aj bj : Eis) : Eis := mul (sub ai bj) (sub bi aj)
end Eis

open Eis in
/-- ω² really is −1 − ω in this arithmetic: ω·ω + ω + 1 = 0 -/
theorem omega_cube_root : add (add (mul ⟨0, 1⟩ ⟨0, 1⟩) ⟨0, 1⟩) ⟨1, 0⟩ = zero := by
  decide +kernel

open Eis in
/-- the waves (1, −1/2), (ω, −ω/2), (ω², −ω²/2): every pair has weight −2 (denominators
    nonzero), and each has weight a/b = −2 to the wave (0, ∞) -/
theorem tetrahedral_waves :
    minusTwo ⟨1, 0⟩ ⟨-1/2, 0⟩ ⟨0, 1⟩ ⟨0, -1/2⟩ = zero ∧
    minusTwo ⟨1, 0⟩ ⟨-1/2, 0⟩ ⟨-1, -1⟩ ⟨1/2, 1/2⟩ = zero ∧
    minusTwo ⟨0, 1⟩ ⟨0, -1/2⟩ ⟨-1, -1⟩ ⟨1/2, 1/2⟩ = zero ∧
    den ⟨1, 0⟩ ⟨-1/2, 0⟩ ⟨0, 1⟩ ⟨0, -1/2⟩ ≠ zero ∧
    den ⟨1, 0⟩ ⟨-1/2, 0⟩ ⟨-1, -1⟩ ⟨1/2, 1/2⟩ ≠ zero ∧
    den ⟨0, 1⟩ ⟨0, -1/2⟩ ⟨-1, -1⟩ ⟨1/2, 1/2⟩ ≠ zero ∧
    add ⟨1, 0⟩ (mul two ⟨-1/2, 0⟩) = zero ∧
    add ⟨0, 1⟩ (mul two ⟨0, -1/2⟩) = zero ∧
    add ⟨-1, -1⟩ (mul two ⟨1/2, 1/2⟩) = zero := by
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩ <;> decide +kernel

end RuleBreakers
