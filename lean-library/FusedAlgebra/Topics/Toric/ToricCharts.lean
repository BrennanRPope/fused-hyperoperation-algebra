/-
  Gluing tesseract charts into P¹ × P¹
  ====================================

  Part of Fused Mult-Add Algebra. Companion to TesseractCoords.lean (same folder). Python check: gluing_check.py.

  STATUS: CHECKED. This file uses core Lean 4 only (no mathlib) and was
  compiled with Lean 4.37.0-rc1 with no errors and no `sorry`. Run it with
  `lean ToricCharts.lean`. It works over any field K: anything with
  Lean's built-in field structure (`Lean.Grind.Field`), which includes the
  rationals, and every mathlib field.

  What it proves, in pictures:

  * A chart is ordinary n-space. Its tesseract corner (which switches are
    on) is the set of nonzero coordinates. Scaling by the torus never
    changes the corner, so corners label the torus orbits.

  * The projective line P¹ is two charts glued by x ↦ 1/x. In mult-add
    coordinates (x = 1 + u) that gluing map is exactly the mult-add
    inverse u ↦ -u/(1+u). The one point the inverse can't reach, u = -1
    (x = 0), is the point the other chart sees as infinity. The 1/x wall
    is where you change charts.

  * P¹ × P¹ is four charts, one per quadrant of its fan (rays +e₁, -e₁,
    +e₂, -e₂). Each chart is injective, together they cover every point,
    the transition maps invert the coordinates that switch, and the
    transitions agree on triple overlaps (the cocycle condition).

  * Each chart's origin is a fixed point that lies in that chart only;
    points of the torus lie in all four charts. Over a field with q
    elements this is the count (q-1)² + 4(q-1) + 4 = (q+1)²: one torus,
    four circles, four fixed points.

  This is gluing at the level of points (charts as injective maps that
  cover, with transition maps). It is not a construction of P¹ × P¹ as a
  scheme or variety, which would need mathlib.

  LIBRARY NOTE: mult-add, its inverse and the wall now come from FusedAlgebra.Base.
-/
import FusedAlgebra.Base

namespace MultAdd.Toric

open Lean.Grind

set_option linter.unusedSectionVars false

variable {K : Type} [Field K] [DecidableEq K]

/-! ## 1. Mult-add and its inverse: the gluing map
  (multAdd, multAddInv, multAdd_multAddInv, inv_one_add and no_multAddInv_at_wall
  are proved once in FusedAlgebra.Base) -/


/-! ## 2. A chart is n-space; its tesseract corner is the set of nonzero coordinates -/

/-- The tesseract corner of a point: switch i is on when coordinate i is nonzero. -/
def corner {n : Nat} (v : Fin n → K) : Fin n → Bool := fun i => decide (v i ≠ 0)

/-- The top corner (all switches on) is the torus. -/
theorem corner_all_on {n : Nat} (v : Fin n → K) :
    (∀ i, corner v i = true) ↔ ∀ i, v i ≠ 0 := by
  simp [corner]

/-- Scaling by the torus (nonzero factors) never moves a point to another
corner, so the corners label the torus orbits. -/
theorem corner_scale {n : Nat} (t v : Fin n → K) (ht : ∀ i, t i ≠ 0) :
    corner (fun i => t i * v i) = corner v := by
  funext i
  have hti := ht i
  have key : t i * v i = 0 ↔ v i = 0 := by
    constructor
    · intro h; grind
    · intro h; rw [h]; grind
  simp [corner, key]

/-! ## 3. The projective line: two charts glued by x ↦ 1/x -/

/-- A point of P¹: `some x` is the number x, `none` is the point at infinity. -/
abbrev P1 (K : Type) := Option K

/-- Chart `false` sees x. Chart `true` sees 1/x, so its origin is infinity.
These are the two cones (rays +e and -e) of the fan of P¹. -/
def chart1 : Bool → K → P1 K
  | false, x => some x
  | true, x => if x = 0 then none else some x⁻¹

theorem chart1_injective (b : Bool) (x y : K) (h : chart1 b x = chart1 b y) : x = y := by
  cases b with
  | false => simpa [chart1] using h
  | true =>
    simp only [chart1] at h
    by_cases hx : x = 0 <;> by_cases hy : y = 0 <;> simp [hx, hy] at h
    · rw [hx, hy]
    · have : x = (x⁻¹)⁻¹ := by grind
      rw [this, h]; grind

/-- The two charts cover P¹. -/
theorem chart1_cover (p : P1 K) : ∃ b x, chart1 b x = p := by
  cases p with
  | none => exact ⟨true, 0, by simp [chart1]⟩
  | some x => exact ⟨false, x, rfl⟩

/-- A point lies in chart b. -/
def InChart1 (b : Bool) (p : P1 K) : Prop := ∃ x, chart1 b x = p

/-- Points of the circle (nonzero, finite) lie in both charts. -/
theorem torus_in_both (x : K) (hx : x ≠ 0) (b : Bool) : InChart1 b (some x) := by
  cases b with
  | false => exact ⟨x, rfl⟩
  | true =>
    refine ⟨x⁻¹, ?_⟩
    have h1 : x⁻¹ ≠ 0 := by grind
    have h2 : (x⁻¹)⁻¹ = x := by grind
    simp [chart1, h1, h2]

/-- The origin of chart c lies in chart c and in no other chart. -/
theorem origin_in_one_chart (c d : Bool) : InChart1 d (chart1 c (0 : K)) ↔ d = c := by
  cases c <;> cases d <;> simp only [InChart1, chart1] <;> simp
  intro x hx h
  have : x * x⁻¹ = 1 := Field.mul_inv_cancel hx
  rw [h] at this; grind

/-- The transition from chart c to chart d on P¹: invert when the chart changes. -/
def trans1 (c d : Bool) (x : K) : K := if c = d then x else x⁻¹

theorem chart1_trans (c d : Bool) (x : K) (h : c ≠ d → x ≠ 0) :
    chart1 d (trans1 c d x) = chart1 c x := by
  cases c <;> cases d <;> simp [trans1, chart1]
  · have hx := h (by decide)
    have h1 : x⁻¹ ≠ 0 := by grind
    have h2 : (x⁻¹)⁻¹ = x := by grind
    simp [h1, h2]
  · have hx := h (by decide)
    simp [hx]

/-! ## 4. P¹ × P¹ from four charts, one per quadrant of the fan -/

/-- The four charts of P¹ × P¹. The chart (c₁, c₂) is the quadrant of the
fan spanned by the rays ±e₁ and ±e₂ that c₁ and c₂ pick. -/
def chart2 (c : Bool × Bool) (v : K × K) : P1 K × P1 K := (chart1 c.1 v.1, chart1 c.2 v.2)

theorem chart2_injective (c : Bool × Bool) (v w : K × K) (h : chart2 c v = chart2 c w) :
    v = w := by
  simp only [chart2, Prod.mk.injEq] at h
  have h1 := chart1_injective c.1 v.1 w.1 h.1
  have h2 := chart1_injective c.2 v.2 w.2 h.2
  cases v; cases w; simp_all

/-- The four charts cover P¹ × P¹. -/
theorem chart2_cover (p : P1 K × P1 K) : ∃ c v, chart2 c v = p := by
  obtain ⟨b1, x1, h1⟩ := chart1_cover p.1
  obtain ⟨b2, x2, h2⟩ := chart1_cover p.2
  exact ⟨(b1, b2), (x1, x2), by simp [chart2, h1, h2]⟩

/-- A point lies in chart c exactly when each factor lies in its chart. -/
def InChart2 (c : Bool × Bool) (p : P1 K × P1 K) : Prop := ∃ v, chart2 c v = p

theorem inChart2_iff (c : Bool × Bool) (p : P1 K × P1 K) :
    InChart2 c p ↔ InChart1 c.1 p.1 ∧ InChart1 c.2 p.2 := by
  constructor
  · rintro ⟨v, hv⟩
    rw [← hv]; exact ⟨⟨v.1, rfl⟩, ⟨v.2, rfl⟩⟩
  · rintro ⟨⟨x, hx⟩, ⟨y, hy⟩⟩
    exact ⟨(x, y), by simp [chart2, hx, hy]⟩

/-- Points of the torus lie in all four charts. -/
theorem torus_in_all (x y : K) (hx : x ≠ 0) (hy : y ≠ 0) (c : Bool × Bool) :
    InChart2 c (some x, some y) :=
  (inChart2_iff c _).2 ⟨torus_in_both x hx c.1, torus_in_both y hy c.2⟩

/-- Each chart's origin is one of the four fixed points, and it lies in
that chart only. -/
theorem fixed_point_in_one_chart (c d : Bool × Bool) :
    InChart2 d (chart2 c ((0 : K), (0 : K))) ↔ d = c := by
  rw [inChart2_iff]
  simp only [chart2, origin_in_one_chart]
  cases c; cases d; simp

/-- Transition between charts: invert each coordinate whose factor switches. -/
def trans2 (c d : Bool × Bool) (v : K × K) : K × K :=
  (trans1 c.1 d.1 v.1, trans1 c.2 d.2 v.2)

/-- Where two charts overlap, the transition carries one chart's
coordinates to the other's. -/
theorem chart2_trans (c d : Bool × Bool) (v : K × K)
    (h1 : c.1 ≠ d.1 → v.1 ≠ 0) (h2 : c.2 ≠ d.2 → v.2 ≠ 0) :
    chart2 d (trans2 c d v) = chart2 c v := by
  simp [chart2, trans2, chart1_trans c.1 d.1 v.1 h1, chart1_trans c.2 d.2 v.2 h2]

/-- The cocycle condition: c → d → e equals c → e, so the gluing is
consistent on every triple overlap. -/
theorem trans1_cocycle (c d e : Bool) (x : K) : trans1 d e (trans1 c d x) = trans1 c e x := by
  cases c <;> cases d <;> cases e <;> simp [trans1] <;> grind

theorem trans2_cocycle (c d e : Bool × Bool) (v : K × K) :
    trans2 d e (trans2 c d v) = trans2 c e v := by
  simp [trans2, trans1_cocycle]

/-! ## 5. The count -/

/-- With a = q - 1 points on the circle: one torus (a²), four circles
(a each), four fixed points. In ordinary coordinates it is (q + 1)². -/
theorem fan_count (a : Nat) : a ^ 2 + 4 * a + 4 = (a + 2) ^ 2 := by grind

end MultAdd.Toric
