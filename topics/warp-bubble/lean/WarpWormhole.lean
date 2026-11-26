/-
  Warp bubbles, wormholes and artificial gravity: the elementary algebra (DRAFT)
  ================================================================================

  Companion checks: ../checks.py (26 of 26 pass). Core Lean 4.37.0-rc1, no mathlib, no `sorry`.
  Imports Base only. NOT YET COMPILED in this session (no Lean toolchain here); sent to the
  shared library for its first compile, same style as GRBridges (Rat + grind).

  1. Zero expansion forces a non-positive cross term (Natário's zero-expansion drive).
  2. A one-direction river with sideways shear has cross term −(s₁² + s₂²): negative energy.
  3. Light in a river moves at v + 1 and v − 1 (plain sum across the river, dial k = 0).
  4. The two-sided cross-ratio of the thermofield double: numerator and denominator in closed form,
     so its odds are (sinh / cosh)².
  5. Gluing regions next to each other: the cross-ratio cocycle, so link probabilities combine by
     the union rule ∘₋₁.
-/
import FusedAlgebra.Base

namespace WarpWormhole
open CircleOp

/-! ## 1. Zero expansion -/

/-- If the three strain rates sum to zero (no volume change), twice their cross term is minus the
sum of their squares, so the Eulerian energy density e2/8π can never be positive. -/
theorem zero_expansion_cross_term (a b c : Rat) (h : a + b + c = 0) :
    2 * (a * b + b * c + c * a) = -(a * a + b * b + c * c) := by
  have hc : c = -a - b := by grind
  subst hc; grind

/-! ## 2. Sideways shear of a one-direction river -/

/-- Strain of the river v = (v_x(x,y,z), 0, 0): K = [[p, s₁, s₂], [s₁, 0, 0], [s₂, 0, 0]].
Its cross term (sum of the three principal 2×2 minors) is −(s₁² + s₂²). -/
theorem shear_cross_term (p s1 s2 : Rat) :
    (p * 0 - s1 * s1) + (p * 0 - s2 * s2) + (0 * 0 - 0 * 0) = -(s1 * s1 + s2 * s2) := by
  grind

/-! ## 3. Light in the river -/

/-- Null directions of −dt² + (dx − v dt)²: (w − v)² − 1 factors as (w − (v+1))(w − (v−1)). -/
theorem river_light (v w : Rat) : (w - v) * (w - v) - 1 = (w - (v + 1)) * (w - (v - 1)) := by
  grind

/-! ## 4. Two-sided cross-ratio in the thermofield double -/

/-- With p = e^{k x₁}, q = e^{k x₂}, u = e^{k t} (k = 2π/β), the left points are p/u, q/u and the
right points −q u, −p u. Cleared of the common factor u², the cross-ratio's denominator
(p + q u²)(q + p u²) splits as (p − q)² u² + p q (1 + u²)²: the first part is the numerator
z₁₂ z₃₄ (times u²), the second is 1 − η's numerator. So the odds η/(1−η) = (p−q)² u² / (p q (1+u²)²)
= sinh²(πℓ/β) / cosh²(2πt/β). -/
theorem tfd_cross_ratio_split (p q u : Rat) :
    (p + q * u * u) * (q + p * u * u) = (p - q) * (p - q) * (u * u) + p * q * ((1 + u * u) * (1 + u * u)) := by
  grind

/-! ## 5. Gluing next-door regions: the cross-ratio cocycle -/

/-- One minus the cross-ratio of A = [a₁, a₂] and B = [b₁, b₂] is N/D with
N = (a₁ − b₂)(a₂ − b₁), D = (a₁ − b₁)(a₂ − b₂). For next-door B₁ = [b₁, b₂], B₂ = [b₂, b₃]:
N₁ N₂ D₁₃ = N₁₃ D₁ D₂, i.e. (1 − x₁)(1 − x₂) = 1 − x₁₃. -/
theorem cross_ratio_cocycle (a1 a2 b1 b2 b3 : Rat) :
    ((a1 - b2) * (a2 - b1)) * ((a1 - b3) * (a2 - b2)) * ((a1 - b1) * (a2 - b3)) =
    ((a1 - b3) * (a2 - b1)) * ((a1 - b1) * (a2 - b2)) * ((a1 - b2) * (a2 - b3)) := by
  grind

/-- And (1 − x)(1 − y) = 1 − z says z = x ∘₋₁ y: the union rule. -/
theorem glue_is_union (x y z : Rat) (h : (1 - x) * (1 - y) = 1 - z) : z = circK (-1) x y := by
  rw [circK_def]; grind

end WarpWormhole
