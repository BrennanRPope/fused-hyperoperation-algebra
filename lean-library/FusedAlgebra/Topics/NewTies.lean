/-
  New ties, round 2: the union fractal. Core Lean 4.37 only (no mathlib), no `sorry`.
  Self-contained (no imports), so it can sit in the shared library as FusedAlgebra/Topics/NewTies.lean.

  1. A union in a box of size M is a stretch about the point M (the box full).
  2. Stretches about one centre commute; about two centres they differ by (c₁ − c₂)(1 − s₁)(1 − s₂).
  3. So two stretches commute exactly when they share a centre (slopes ≠ 1).
  4. The key identity of the lower bound: (1 − P)(1 − Q) = PQ ↔ P + Q = 1.
  5. With one map per box the dimension equation is Moran's.
  6. A small exact count: three integer stretches, two sharing a centre, give 3, 8, 21 distinct composites
     (the shared-centre count), not 3, 9, 27.
  7. Plane laws: a law's maps fix w exactly when d·w² + k·w + 1 = 0, and that line is tangent to the border.

  First compiled in the shared library on Nov 27, as FusedAlgebra.Topics.NewTies
  (namespace CircleOp.NewTies, the library's convention; no proof changed).
-/
namespace CircleOp.NewTies

/-- A stretch about the centre c by the factor s: x ↦ c + s(x − c). -/
def hom (c s x : Int) : Int := c + s * (x - c)

/-- Every stretch fixes its centre. -/
theorem hom_centre (c s : Int) : hom c s c = c := by simp only [hom]; grind

/-- 1. A union in a box of size M, scaled by M: M·(a + x − ax/M) = M·M + (M − a)(x − M),
    i.e. the stretch about M (the full box) by the factor 1 − a/M. -/
theorem union_is_stretch (M a x : Int) : M * a + M * x - a * x = M * M + (M - a) * (x - M) := by grind

/-- 2a. Stretches about one centre commute. -/
theorem same_centre_commute (c s t x : Int) : hom c s (hom c t x) = hom c t (hom c s x) := by
  simp only [hom]; grind

/-- 2b. The order gap between stretches about two centres. -/
theorem order_gap (c₁ s₁ c₂ s₂ x : Int) :
    hom c₁ s₁ (hom c₂ s₂ x) - hom c₂ s₂ (hom c₁ s₁ x) = (c₁ - c₂) * ((1 - s₁) * (1 - s₂)) := by
  simp only [hom]; grind

/-- 3. Two stretches with slopes ≠ 1 commute exactly when they share a centre. -/
theorem commute_iff (c₁ s₁ c₂ s₂ : Int) (h₁ : s₁ ≠ 1) (h₂ : s₂ ≠ 1) :
    (∀ x, hom c₁ s₁ (hom c₂ s₂ x) = hom c₂ s₂ (hom c₁ s₁ x)) ↔ c₁ = c₂ := by
  constructor
  · intro h
    have g := order_gap c₁ s₁ c₂ s₂ 0
    rw [h 0, Int.sub_self] at g
    have z : (c₁ - c₂) * ((1 - s₁) * (1 - s₂)) = 0 := g.symm
    rcases Int.mul_eq_zero.mp z with e | e
    · grind
    · rcases Int.mul_eq_zero.mp e with e' | e'
      · exact absurd (by grind) h₁
      · exact absurd (by grind) h₂
  · intro e x; subst e; exact same_centre_commute c₁ s₁ s₂ x

/-- 4. The lower bound's key step: the first-return generators have total weight 1 exactly when
    the shared-centre dimension equation holds, (1/P − 1)(1/Q − 1) = 1 ↔ P + Q = 1, cleared of fractions. -/
theorem lower_bound_key (P Q : Int) : (1 - P) * (1 - Q) = P * Q ↔ P + Q = 1 := by
  constructor <;> intro h <;> grind

/-- 5. One map per box: 1 + Σ[(1 − wᵢ) − 1] = 1 − Σ wᵢ, Moran's equation. -/
theorem moran_special (ws : List Int) :
    1 + (ws.map (fun w => (1 - w) - 1)).sum = 1 - ws.sum := by
  induction ws with
  | nil => simp
  | cons w ws ih => simp only [List.map_cons, List.sum_cons]; grind

/-! ## 6. An exact count -/

/-- A map x ↦ s·x + b as a pair (s, b); composing p after q. -/
def after (p q : Int × Int) : Int × Int := (p.1 * q.1, p.1 * q.2 + p.2)

/-- Three integer stretches: about 0 by 2, about 0 by 3, about 1 by 5 (x ↦ 5x − 4). -/
def gens : List (Int × Int) := [(2, 0), (3, 0), (5, -4)]

/-- All composites of length n (with repeats). -/
def words : Nat → List (Int × Int)
  | 0 => [(1, 0)]
  | n + 1 => (words n).flatMap (fun w => gens.map (fun g => after g w))

/-- Distinct composites of length 1, 2, 3: 3, 8, 21, the shared-centre count (not 3, 9, 27). -/
theorem count_one : (words 1).eraseDups.length = 3 := by decide
theorem count_two : (words 2).eraseDups.length = 8 := by decide
theorem count_three : (words 3).eraseDups.length = 21 := by decide

/-! ## 7. Plane laws: sharing a wall is a tangent line -/

/-- The plane-law map x ↦ (x + a + kax)/(1 − dax) fixes w exactly when d·w² + k·w + 1 = 0 (cleared of the fraction). -/
theorem fixes_wall (k d a w : Int) :
    w * (1 - d * a * w) - (w + a + k * a * w) = -(a * (d * w * w + k * w + 1)) := by grind

/-- On the border 4d = k², the wall line d·w² + k·w + 1 = 0 becomes (k·w + 2)² = 0 (times 4): a double root,
    so the laws sharing the wall w form the tangent line to the border at k = −2/w. -/
theorem wall_line_tangent (k w : Int) : k * k * w * w + 4 * (k * w) + 4 = (k * w + 2) * (k * w + 2) := by grind

end CircleOp.NewTies
