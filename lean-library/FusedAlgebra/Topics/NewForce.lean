/-
DRAFT (Nov 26, 2026). Does the theory predict an unknown force? The algebra behind
topics/force-reach/new-force/README.md, compiled in the shared library on Nov 26 as FusedAlgebra.Topics.NewForce.

The session that wrote this stated it against Mathlib over ℝ (topics/force-reach/new-force/lean/NewForce.lean,
fourteen statements; twenty from Nov 27; thirty-five with parts J and K). Here all but one are proved in core Lean over any field of characteristic 0 or, for the two
flavour bounds and γ from the scalar strength, any ordered field. Statements and hypotheses are as written, except that parts J and K state π as any c (see that section).
-/
import FusedAlgebra.Base

open Lean.Grind

namespace CircleOp.NewForce

section Field
variable {K : Type} [Field K] [IsCharP K 0]

theorem tower_sum (q : K) (hq : q ≠ 1) :
    1 + 2 * q / (1 - q) = (1 + q) / (1 - q) := by
  have h : (1 - q) ≠ 0 := by intro h; apply hq; grind
  grind

theorem tower_leading_yukawa (q : K) (hq : q ≠ 1) :
    (1 + q) / (1 - q) - 1 - 2 * q = 2 * q ^ 2 / (1 - q) := by
  have h : (1 - q) ≠ 0 := by intro h; apply hq; grind
  grind

theorem dual_exponent_mond : (3 / 2 : K) / (3 / 2 - 1) = 3 := by grind
theorem holder_pair : (1 : K) / (3 / 2) + 1 / 3 = 1 := by grind
omit [IsCharP K 0] in
theorem vector_over_tensor : (2 : K) / (1 / 2) = 4 := by grind
theorem scalar_over_tensor : (3 / 2 : K) / (1 / 2) = 3 := by grind
theorem like_charges_cancel : (1 / 2 : K) + 3 / 2 - 2 = 0 := by grind
theorem opposite_charges_add : (1 / 2 : K) + 3 / 2 + 2 = 4 := by grind
theorem klein_strength (α : K) (hα : α ≠ 0) : 4 / (4 / α) = α := by grind
theorem one_circle_third : (1 : K) / (1 + 2) = 1 / 3 := by grind
theorem one_circle_gamma : (1 : K) / (1 + 1) - 1 = -1 / 2 := by grind
end Field

section Ordered
variable {K : Type} [LE K] [LT K] [Field K] [Std.IsLinearOrder K] [Std.LawfulOrderLT K] [OrderedRing K]

theorem su3_flavour_bound :
    (0 : K) < 11 * 3 / 3 - 2 * 16 / 3 ∧ 11 * 3 / 3 - 2 * 17 / 3 < (0 : K) := by
  constructor <;> grind

theorem su2_flavour_bound :
    (0 : K) < 11 * 2 / 3 - 2 * 10 / 3 ∧ 11 * 2 / 3 - 2 * 11 / 3 = (0 : K) := by
  constructor <;> grind

theorem gamma_from_strength (n : K) (hn : 0 ≤ n) :
    (1 - n / (n + 2)) / (1 + n / (n + 2)) = 1 / (n + 1) := by
  have h2 : n + 2 ≠ 0 := by grind
  have h1 : n + 1 ≠ 0 := by grind
  have h3 : 2 * n + 2 ≠ 0 := by grind
  grind
end Ordered

section Hold
variable {K : Type} [Field K] [IsCharP K 0]

/-! ### Part H: holding the circle (vacuum energy plus a 5D vacuum energy); added Nov 27 with sections 7 and 8 -/

/-- With V = a/R + C R*²/R⁶ and a = −6C/R*³, R* is a turning point: V′(R*) = 0. -/
theorem hold_turning_point (C R : K) (hR : R ≠ 0) :
    -(-6 * C / R ^ 3) / R ^ 2 - 6 * C * R ^ 2 / R ^ 7 = 0 := by
  grind

/-- Its curvature there: V″(R*) = 30C/R*⁶ (a minimum when C > 0). -/
theorem hold_curvature (C R : K) (hR : R ≠ 0) :
    2 * (-6 * C / R ^ 3) / R ^ 3 + 42 * C * R ^ 2 / R ^ 8 = 30 * C / R ^ 6 := by
  grind

/-- The vacuum energy left over: V(R*) = −5C/R*⁴. -/
theorem hold_leftover (C R : K) (hR : R ≠ 0) :
    (-6 * C / R ^ 3) / R + C * R ^ 2 / R ^ 6 = -5 * C / R ^ 4 := by
  grind

/-- Radion mass: (2/3)·30 = 20, so m² = 20C/(M̄²R⁴). -/
theorem hold_mass_factor : (2 / 3 : K) * 30 = 20 := by grind

/-! ### Part S: strength from shape and material -/

/-- Sphere: ⟨|K|²⟩ = 2R²/3, so the shape factor 4R²/⟨|K|²⟩ is 6 (a circle gives 4). -/
theorem sphere_shape_factor : (4 : K) / (2 / 3) = 6 := by grind

omit [IsCharP K 0] in
/-- Strong force with 5 flavours: b₀ = 11·3/3 − 2·5/3 = 23/3. -/
theorem strong_b0_five : (11 * 3 / 3 : K) - 2 * 5 / 3 = 23 / 3 := by grind
end Hold

/-! ### Part I: the new force up close (added Nov 27 from the topic file's Part I; same
statements over any field of characteristic 0, proved by grind in place of norm_num) -/
section PartI
variable {K : Type} [Field K] [IsCharP K 0]

omit [IsCharP K 0] in
/-- Each rung between static masses: 5D weight 2/3 over 4D weight 1/2 is 4/3. -/
theorem every_rung_four_thirds : (2 / 3 : K) / (1 / 2) = 4 / 3 := by grind

/-- The zero rung's 4/3 is gravity (1) plus the new force (1/3). -/
theorem zero_rung_split : (1 : K) + 1 / 3 = 4 / 3 := by grind

omit [IsCharP K 0] in
/-- Copies sit on the extreme-black-hole line: charge²/(G mass²) = 1 + a² = 1 + 3 = 4. -/
theorem copies_extreme : (1 : K) + 3 = 4 := by grind

/-- Near a static mass the circle swells: h₅₅/h₀₀ = (1/3)/(2/3) = 1/2. -/
theorem circle_swell_ratio : (1 / 3 : K) / (2 / 3) = 1 / 2 := by grind

omit [IsCharP K 0] in
/-- Electric strength dips twice as much: −2 · (2/3) = −4/3. -/
theorem alpha_dip : (-2 : K) * (2 / 3) = -4 / 3 := by grind

/-- The circle's potential tail e^{−√(2/3) σ} has slope parameter (1/2)(2/3) = 1/3. -/
theorem inflaton_tail_eps : (1 / 2 : K) * (2 / 3) = 1 / 3 := by grind
end PartI

/-! ### Parts J and K: the inverse, passing for gravity, and the other unknown forces (added Nov 27
from the topic file; over any field of characteristic 0, proved by grind in place of norm_num and
field_simp; where the topic file has Real.pi, the statement here holds for any nonzero c, π included.
`square_torus_no_mixing` (cos(π/2) = 0) needs the real cosine, so it stays in the topic file only,
built by CI's report-only Mathlib step) -/
section PartJK
variable {K : Type} [Field K] [IsCharP K 0]

omit [IsCharP K 0] in
/-- Beyond the size-carrier's range two like copies feel gravity (+1) and the vector (−4): a push of 3. -/
theorem like_copies_far_repel : (1 : K) + 0 - 4 = -3 := by grind

/-- Matter (size charge² 1/3) and a copy (size charge² 3): the cross weight is (1/3)·3 = 1. -/
theorem matter_copy_size : (1 / 3 : K) * 3 = 1 := by grind

/-- At the cap R₂ = R₁³/(4c²) the first circle equals the 6D fundamental length: 4c²R₁R₂ = R₁⁴
    (c = π in the write-up). -/
theorem second_circle_cap (c R₁ : K) (hc : c ≠ 0) :
    4 * c ^ 2 * R₁ * (R₁ ^ 3 / (4 * c ^ 2)) = R₁ ^ 4 := by grind

/-- The cap as a ratio of the two circles: (R₁³/(4c²))/R₁ = R₁²/(4c²) (c = π in the write-up). -/
theorem cap_ratio (c R₁ : K) (hR : R₁ ≠ 0) :
    (R₁ ^ 3 / (4 * c ^ 2)) / R₁ = R₁ ^ 2 / (4 * c ^ 2) := by grind

/-- A 1/3 force on matter but not light: light bends by (1 + γ)/2 = 3/4 of the expected amount. -/
theorem light_bends_three_quarters :
    (1 / 2 : K) * (1 + (1 - 1 / 3) / (1 + 1 / 3)) = 3 / 4 := by grind

omit [IsCharP K 0] in
/-- SU(3) with 16 flavours: b₀ = 11 − 32/3 = 1/3, still asymptotically free. -/
theorem bz_b0_sixteen : (11 * 3 / 3 : K) - 2 * 16 / 3 = 1 / 3 := by grind

/-- SU(3) with 16 flavours: b₁ = 102 − 38·16/3 = −302/3 < 0, so a fixed point exists. -/
theorem bz_b1_sixteen : (102 : K) - 38 * 16 / 3 = -302 / 3 := by grind

/-- The frozen strength: α* = −4c b₀/b₁ = 4c/302 at 16 flavours (c = π in the write-up). -/
theorem bz_alpha_sixteen (c : K) : -4 * c * (1 / 3) / (-302 / 3 : K) = 4 * c / 302 := by grind
end PartJK

end CircleOp.NewForce
