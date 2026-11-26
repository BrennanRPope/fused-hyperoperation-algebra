/-
Riemann round 8 (DRAFT): small exact identities behind the round's pictures.
Core Lean + grind, in the style of earlier rounds.

First compiled in the shared library on Nov 26, as FusedAlgebra.Topics.RiemannRound8
(namespace CircleOp.RiemannRound8, the library's convention; no proof changed).
-/
import FusedAlgebra.Base

namespace CircleOp.RiemannRound8

/-- One rung of the twisted ladder: a glide's even and odd ladders multiply to the cover's rung. -/
theorem rung (x : Int) : (1 - x) * (1 + x) = 1 - x * x := by grind

/-- The orientation double cover of k crosscaps is orientable of genus k − 1 (χ doubles). -/
theorem cover_chi (k : Int) : 2 * (2 - k) = 2 - 2 * (k - 1) := by grind

/-- Gauss–Bonnet at curvature −1: area/π = −2χ = 4g − 4 for genus g. -/
theorem area_genus (g : Int) : -2 * (2 - 2 * g) = 4 * g - 4 := by grind

/-- Klein bottle: even sector (2a − b) plus odd sector (2a + b) gives the torus 4a. -/
theorem klein_sectors (a b : Int) : (2 * a - b) + (2 * a + b) = 4 * a := by grind

/-- Two-particle heat flow: velocities ±2/d give relative speed 4/d, and −2·4·d = −8·d
    (the cleared form of z̈ = −8/d³). -/
theorem cm_two (d : Int) : -2 * (2 - (-2)) * d = -8 * d := by grind

/-- The 2×2 Calogero–Moser Lax matrix (x, x; −x, −x) squares to zero: nilpotent. -/
theorem lax_two_nilpotent (x : Int) :
    x * x + x * (-x) = 0 ∧ x * x + x * (-x) = 0 ∧
    (-x) * x + (-x) * (-x) = 0 ∧ (-x) * x + (-x) * (-x) = 0 := by grind

/-- Rewinding z² − a: the zeros meet when a + 2t = 0; with half-gap² a = d²/4 that is t = −d²/8. -/
theorem pair_collision (a d t : Int) (h1 : 4 * a = d * d) (h2 : 8 * t = -(d * d)) :
    a + 2 * t = 0 := by grind

/-- Fusion is multiplication after the shift a ↦ 1 + k a. -/
theorem fused_shift (k a b : Int) :
    (1 + k * a) * (1 + k * b) = 1 + k * (a + b + k * a * b) := by grind

/-- The fused centre absorbs: if k c = −1 then a ∘ₖ c = c. -/
theorem fused_absorb (k a c : Int) (h : k * c = -1) : a + c + k * a * c = c := by
  have : k * a * c = -a := by
    calc k * a * c = a * (k * c) := by grind
      _ = -a := by rw [h]; grind
  grind

/-- Fermi end of the fused residue: (1 + 1/p)(1 − 1/p) = 1 − 1/p², cleared. -/
theorem fermi_residue (p : Int) : (p + 1) * (p - 1) = p * p - 1 := by grind

/-- Möbius re-centring w = 1/(s − 1): |w + 1| = 1 ⇔ |s| = |s − 1| ⇔ 2 Re s = 1. -/
theorem mobius_circle (x y : Int) :
    x * x + y * y = (x - 1) * (x - 1) + y * y ↔ 2 * x = 1 := by
  constructor <;> intro h <;> grind

/-- Flat glide weight: the transverse map of a glide is −1, and |1 − (−1)| = 2. -/
theorem glide_det : (1 : Int) - (-1) = 2 := by decide

end CircleOp.RiemannRound8
