/-
  Affine Composition Theory: the algebra behind the Yang–Mills gap DRAFT.
  Core Lean 4.37 only (no mathlib), no `sorry`. Imports FusedAlgebra.Base.
  Intended library path: FusedAlgebra/Topics/YangMillsGap.lean

  STATUS: first compiled in the shared library on Nov 26 2026, as FusedAlgebra.Topics.YangMillsGap.
  Written in the style of the compiled LeeYang.lean (CommRing + grind, Rat.mul_nonneg).

  What is stated (the analysis parts: circumference 2π/√(−Δ), Laplacian spectrum, zeros,
  transfer-matrix gaps, are numerical, in checks.py, and are not stated here):
  1. Completing the square for the generator v(x) = 1 + kx + dx²:
       4d·v(x) = (2dx + k)² − Δ,  Δ = k² − 4d.
     So for d > 0 and Δ < 0 the generator never vanishes: no wall, the circle class.
  2. (The gap reading −Δ = 4 · (−Δ/4) is bookkeeping, not stated.)
  3. The parity dial k = −2 on {0, 1} is XOR: the Z₂ gauge group as a dial law.
  4. Abelian plaquette: with a shared dial the four links compose in any order.
  5. Twist ladder: levels (n + t)²; at half a rung (t = 1/2) rungs n and −1−n tie,
     and for 0 ≤ t ≤ 1/2 the gap (1 − t)² − t² = 1 − 2t.
  6. Rung costs: Casimir j(j+1) ratios, C₂(1)/C₂(1/2) = 8/3; U(1) rung n costs n².
  7. Lift form of the plaquette: (1 + kA)(1 + kB) − (1 + kB)(1 + kA) = k²(AB − BA)
     in any ring: the order gap is k² times the commutator (non-abelian links).
-/
import FusedAlgebra.Base

namespace YangMillsGap
open CircleOp

set_option linter.unusedSectionVars false

section CommRingPart
variable {R : Type} [Lean.Grind.CommRing R]

/-- 1. Completing the square for the dial generator. -/
theorem generator_square (k d x : R) :
    4 * d * (1 + k * x + d * x * x) = (2 * d * x + k) * (2 * d * x + k) - (k * k - 4 * d) := by
  grind

/-- 5a. At half a rung the rungs n and −1−n tie (twist t = 1/2, written with 2t = 1). -/
theorem half_rung_tie (n t : R) (ht : 2 * t = 1) :
    (n + t) * (n + t) = ((-1 - n) + t) * ((-1 - n) + t) := by
  grind

/-- 5b. The twisted gap between rung −1 and rung 0: (1 − t)² − t² = 1 − 2t. -/
theorem twisted_gap (t : R) : (1 - t) * (1 - t) - t * t = 1 - 2 * t := by
  grind

/-- 4. Abelian plaquette: shared dial k, any order of the four links gives the same holonomy. -/
theorem abelian_plaquette_order_free (k a b c e : R) :
    circK k (circK k (circK k a b) c) e = circK k (circK k (circK k e c) b) a := by
  unfold circK; grind

end CommRingPart

section NonCommutative
variable {R : Type} [Lean.Grind.Ring R]

/-- 7. Non-abelian links: the order gap of the lifts is k² times the commutator.
    Stated for an arbitrary (possibly non-commutative) ring via `Lean.Grind.Ring`.
    (Compiled Nov 26: the products are expanded by hand before `grind`, which
    does not normalize non-commutative products.) -/
theorem lift_order_gap (k A B : R) (hk : ∀ X : R, k * X = X * k) :
    (1 + k * A) * (1 + k * B) - (1 + k * B) * (1 + k * A) = k * k * (A * B - B * A) := by
  have hA : k * A * (k * B) = k * k * (A * B) := by
    rw [Lean.Grind.Semiring.mul_assoc, ← Lean.Grind.Semiring.mul_assoc A, ← hk A,
      Lean.Grind.Semiring.mul_assoc, Lean.Grind.Semiring.mul_assoc]
  have hB : k * B * (k * A) = k * k * (B * A) := by
    rw [Lean.Grind.Semiring.mul_assoc, ← Lean.Grind.Semiring.mul_assoc B, ← hk B,
      Lean.Grind.Semiring.mul_assoc, Lean.Grind.Semiring.mul_assoc]
  simp only [Lean.Grind.Semiring.left_distrib, Lean.Grind.Semiring.right_distrib,
    Lean.Grind.Semiring.mul_one, Lean.Grind.Semiring.one_mul, hA, hB]
  grind
end NonCommutative

/-- Squares are nonnegative over Rat (same proof as DialPointProcesses.sq_nonneg'). -/
theorem rat_sq_nonneg (s : Rat) : 0 ≤ s * s := by
  rcases (@Rat.le_total 0 s) with h | h
  · exact Rat.mul_nonneg h h
  · have h1 : 0 ≤ -s := by grind
    have := Rat.mul_nonneg h1 h1
    grind

/-- 1'. Over the rationals: d > 0 and Δ < 0 force v(x) > 0 for every x (no wall, circle class). -/
theorem circle_class_no_wall (k d x : Rat) (hd : 0 < d) (hD : k * k - 4 * d < 0) :
    0 < 1 + k * x + d * x * x := by
  have sq := rat_sq_nonneg (2 * d * x + k)
  have key : 4 * d * (1 + k * x + d * x * x) = (2 * d * x + k) * (2 * d * x + k) - (k * k - 4 * d) := by
    grind
  have pos : 0 < 4 * d * (1 + k * x + d * x * x) := by grind
  rcases (@Rat.le_total (1 + k * x + d * x * x) 0) with hle | hge
  · -- v ≤ 0 would make 4d·v ≤ 0
    have h4 : 0 ≤ 4 * d := by grind
    have hn : 0 ≤ -(1 + k * x + d * x * x) := by grind
    have := Rat.mul_nonneg h4 hn
    grind
  · rcases (Rat.le_iff_lt_or_eq.mp hge) with hlt | heq
    · exact hlt
    · rw [← heq] at pos; grind

/-- 3. The parity dial k = −2 on {0, 1} is XOR (the Z₂ gauge group). -/
theorem parity_dial_xor :
    circK (-2 : Int) 0 0 = 0 ∧ circK (-2 : Int) 0 1 = 1 ∧
    circK (-2 : Int) 1 0 = 1 ∧ circK (-2 : Int) 1 1 = 0 := by
  unfold circK; decide

/-- 6. Rung costs. SU(2): C₂(j) = j(j+1); C₂(1)/C₂(1/2) = 8/3. U(1): rung n costs n². -/
theorem casimir_ratio : ((1 : Rat) * (1 + 1)) = (8 / 3) * ((1 / 2) * (1 / 2 + 1)) := by
  grind

theorem u1_rung_ratio : ((2 : Int) * 2) = 4 * (1 * 1) := by decide

end YangMillsGap
