module

public import FusedAlgebra.Basic

/-!
# Gluing tesseract charts into `P¹ × P¹`

Gluing at the level of points, over a field `K`: charts are injective maps that cover, with
transition maps that satisfy the cocycle condition.

* A chart is ordinary `n`-space. Its tesseract `corner` (which switches are on) is the set of
  nonzero coordinates, that is `Function.support` (`corner_eq_true_iff`). Scaling by the torus
  never changes the corner (`corner_scale`), so corners label the torus orbits.
* The projective line `P1 K` is two charts glued by `x ↦ x⁻¹` (`chart1`, `trans1`). In mult-add
  coordinates `x = 1 + u` the gluing map is the mult-add inverse (`trans1_one_add`), and the one
  point it cannot reach, the wall `u = -1`, is the point the other chart sees as infinity.
* `P¹ × P¹` is four charts, one per quadrant of its fan (`chart2`, `trans2`): each chart is
  injective, together they cover, the transitions invert the coordinates that switch and agree
  on triple overlaps (`trans2_cocycle`).
* Each chart's origin is a fixed point that lies in that chart only; the torus lies in all four.
  Over a field with `q` elements this is the count `(q - 1) ^ 2 + 4 * (q - 1) + 4 = (q + 1) ^ 2`
  (`fan_count`, `card_P1_prod_P1`).

## Relation to mathlib

The projective line here is `P1 K = OnePoint K` (the core-Lean original used `Option K`, which
is the same type). Mathlib identifies it with the projectivization of `K²`
(`OnePoint.equivProjectivization`), restated here as `P1.equivProjectivization`; under it chart
`false` is `x ↦ [x : 1]` and chart `true` is `x ↦ [1 : x]` (`equivProjectivization_chart1_false`,
`equivProjectivization_chart1_true`). So the two charts are the standard affine cover of
`Projectivization K (Fin 2 → K)`.

This file does not build `P¹ × P¹` as a scheme. Mathlib's scheme-level projective space is
`AlgebraicGeometry.ProjectiveSpace`; the `K`-points of `ℙ(1; Spec K)` are the points of `P1 K`,
but no statement here goes through schemes.
-/

@[expose] public section

namespace MultAdd.Toric

open OnePoint

/-! ### A chart is `n`-space; its tesseract corner is the set of nonzero coordinates -/

section Corner

variable {M : Type*} [MulZeroClass M] [DecidableEq M] {n : ℕ}

/-- The tesseract corner of a point: switch `i` is on when coordinate `i` is nonzero. -/
def corner (v : Fin n → M) : Fin n → Bool := fun i ↦ decide (v i ≠ 0)

/-- The corner is the support of the point. -/
theorem corner_eq_true_iff (v : Fin n → M) (i : Fin n) :
    corner v i = true ↔ i ∈ Function.support v := by
  simp [corner]

/-- The top corner (all switches on) is the torus. -/
theorem corner_all_on (v : Fin n → M) : (∀ i, corner v i = true) ↔ ∀ i, v i ≠ 0 := by
  simp [corner]

/-- Scaling by the torus (nonzero factors) never moves a point to another corner, so the corners
label the torus orbits. -/
theorem corner_scale [NoZeroDivisors M] (t v : Fin n → M) (ht : ∀ i, t i ≠ 0) :
    corner (fun i ↦ t i * v i) = corner v := by
  funext i
  simp [corner, ht i]

end Corner

variable {K : Type*} [Field K] [DecidableEq K]

/-! ### The projective line: two charts glued by `x ↦ x⁻¹` -/

/-- A point of `P¹`: `(x : P1 K)` is the number `x`, and `∞` is the point at infinity. This is
mathlib's one-point compactification `OnePoint K`. -/
abbrev P1 (K : Type*) := OnePoint K

/-- The projective line is the projectivization of `K²`: mathlib's
`OnePoint.equivProjectivization`. -/
def P1.equivProjectivization : P1 K ≃ Projectivization K (Fin 2 → K) :=
  OnePoint.equivProjectivization K

/-- Chart `false` sees `x`. Chart `true` sees `x⁻¹`, so its origin is infinity. These are the two
cones (rays `+e` and `-e`) of the fan of `P¹`. -/
def chart1 : Bool → K → P1 K
  | false, x => (x : OnePoint K)
  | true, x => if x = 0 then ∞ else ((x⁻¹ : K) : OnePoint K)

theorem chart1_false (x : K) : chart1 false x = (x : OnePoint K) := rfl

theorem chart1_true_zero : chart1 true (0 : K) = ∞ := ite_eq_left rfl

theorem chart1_true_of_ne_zero {x : K} (hx : x ≠ 0) :
    chart1 true x = ((x⁻¹ : K) : OnePoint K) := ite_eq_right hx

/-- In homogeneous coordinates chart `false` is `x ↦ [x : 1]`. -/
theorem equivProjectivization_chart1_false (x : K) :
    P1.equivProjectivization (chart1 false x) = Projectivization.mk K ![x, 1] (by simp) :=
  rfl

/-- In homogeneous coordinates chart `true` is `x ↦ [1 : x]`. -/
theorem equivProjectivization_chart1_true (x : K) :
    P1.equivProjectivization (chart1 true x) = Projectivization.mk K ![1, x] (by simp) := by
  by_cases hx : x = 0
  · subst hx; rw [chart1_true_zero]; rfl
  · rw [chart1_true_of_ne_zero hx, P1.equivProjectivization,
      OnePoint.equivProjectivization_apply_coe, Projectivization.mk_eq_mk_iff']
    exact ⟨x⁻¹, by ext i; fin_cases i <;> simp [hx]⟩

theorem chart1_injective (b : Bool) : Function.Injective (chart1 (K := K) b) := by
  intro x y h
  cases b with
  | false => exact OnePoint.coe_injective h
  | true =>
    by_cases hx : x = 0 <;> by_cases hy : y = 0
    · rw [hx, hy]
    · rw [hx, chart1_true_zero, chart1_true_of_ne_zero hy] at h
      exact absurd h.symm (OnePoint.coe_ne_infty _)
    · rw [hy, chart1_true_zero, chart1_true_of_ne_zero hx] at h
      exact absurd h (OnePoint.coe_ne_infty _)
    · rw [chart1_true_of_ne_zero hx, chart1_true_of_ne_zero hy] at h
      exact inv_injective (OnePoint.coe_injective h)

/-- The two charts cover `P¹`. -/
theorem chart1_cover (p : P1 K) : ∃ b x, chart1 b x = p := by
  induction p using OnePoint.rec with
  | infty => exact ⟨true, 0, chart1_true_zero⟩
  | coe x => exact ⟨false, x, rfl⟩

/-- A point lies in chart `b`. -/
def InChart1 (b : Bool) (p : P1 K) : Prop := ∃ x, chart1 b x = p

/-- Points of the circle (nonzero, finite) lie in both charts. -/
theorem torus_in_both (x : K) (hx : x ≠ 0) (b : Bool) : InChart1 b (x : OnePoint K) := by
  cases b with
  | false => exact ⟨x, rfl⟩
  | true => exact ⟨x⁻¹, by rw [chart1_true_of_ne_zero (inv_ne_zero hx), inv_inv]⟩

/-- The origin of chart `c` lies in chart `c` and in no other chart. -/
theorem origin_in_one_chart (c d : Bool) : InChart1 d (chart1 c (0 : K)) ↔ d = c := by
  constructor
  · rintro ⟨x, hx⟩
    cases c <;> cases d <;> try rfl
    · rw [chart1_false] at hx
      by_cases h : x = 0
      · rw [h, chart1_true_zero] at hx
        exact absurd hx.symm (OnePoint.coe_ne_infty _)
      · rw [chart1_true_of_ne_zero h] at hx
        exact absurd (OnePoint.coe_injective hx) (inv_ne_zero h)
    · rw [chart1_true_zero, chart1_false] at hx
      exact absurd hx (OnePoint.coe_ne_infty _)
  · rintro rfl
    exact ⟨0, rfl⟩

/-- The transition from chart `c` to chart `d` on `P¹`: invert when the chart changes. -/
def trans1 (c d : Bool) (x : K) : K := if c = d then x else x⁻¹

theorem chart1_trans (c d : Bool) (x : K) (h : c ≠ d → x ≠ 0) :
    chart1 d (trans1 c d x) = chart1 c x := by
  cases c <;> cases d
  · simp [trans1]
  · have hx : x ≠ 0 := h (by decide)
    simp [trans1, chart1_true_of_ne_zero (inv_ne_zero hx), chart1_false]
  · have hx : x ≠ 0 := h (by decide)
    simp [trans1, chart1_true_of_ne_zero hx, chart1_false]
  · simp [trans1]

omit [DecidableEq K] in
/-- In mult-add coordinates `x = 1 + u` the gluing map of `P¹` is the mult-add inverse. The wall
`u = -1` is the point the other chart sees as infinity. -/
theorem trans1_one_add (u : K) (h : 1 + u ≠ 0) :
    trans1 false true (1 + u) = 1 + multAddInv u := by
  simp [trans1, inv_one_add u h]

/-! ### `P¹ × P¹` from four charts, one per quadrant of the fan -/

/-- The four charts of `P¹ × P¹`. The chart `(c₁, c₂)` is the quadrant of the fan spanned by the
rays `±e₁` and `±e₂` that `c₁` and `c₂` pick. -/
def chart2 (c : Bool × Bool) (v : K × K) : P1 K × P1 K := (chart1 c.1 v.1, chart1 c.2 v.2)

theorem chart2_injective (c : Bool × Bool) : Function.Injective (chart2 (K := K) c) :=
  fun _ _ h ↦ Prod.ext (chart1_injective c.1 (congrArg Prod.fst h))
    (chart1_injective c.2 (congrArg Prod.snd h))

/-- The four charts cover `P¹ × P¹`. -/
theorem chart2_cover (p : P1 K × P1 K) : ∃ c v, chart2 c v = p := by
  obtain ⟨b₁, x₁, h₁⟩ := chart1_cover p.1
  obtain ⟨b₂, x₂, h₂⟩ := chart1_cover p.2
  exact ⟨(b₁, b₂), (x₁, x₂), Prod.ext h₁ h₂⟩

/-- A point lies in chart `c`. -/
def InChart2 (c : Bool × Bool) (p : P1 K × P1 K) : Prop := ∃ v, chart2 c v = p

/-- A point lies in chart `c` exactly when each factor lies in its chart. -/
theorem inChart2_iff (c : Bool × Bool) (p : P1 K × P1 K) :
    InChart2 c p ↔ InChart1 c.1 p.1 ∧ InChart1 c.2 p.2 := by
  constructor
  · rintro ⟨v, rfl⟩
    exact ⟨⟨v.1, rfl⟩, ⟨v.2, rfl⟩⟩
  · rintro ⟨⟨x, hx⟩, ⟨y, hy⟩⟩
    exact ⟨(x, y), Prod.ext hx hy⟩

/-- Points of the torus lie in all four charts. -/
theorem torus_in_all (x y : K) (hx : x ≠ 0) (hy : y ≠ 0) (c : Bool × Bool) :
    InChart2 c ((x : OnePoint K), (y : OnePoint K)) :=
  (inChart2_iff c _).2 ⟨torus_in_both x hx c.1, torus_in_both y hy c.2⟩

/-- Each chart's origin is one of the four fixed points, and it lies in that chart only. -/
theorem fixed_point_in_one_chart (c d : Bool × Bool) :
    InChart2 d (chart2 c ((0 : K), (0 : K))) ↔ d = c := by
  rw [inChart2_iff, chart2, origin_in_one_chart, origin_in_one_chart, Prod.ext_iff]

/-- Transition between charts: invert each coordinate whose factor switches. -/
def trans2 (c d : Bool × Bool) (v : K × K) : K × K :=
  (trans1 c.1 d.1 v.1, trans1 c.2 d.2 v.2)

/-- Where two charts overlap, the transition carries one chart's coordinates to the other's. -/
theorem chart2_trans (c d : Bool × Bool) (v : K × K)
    (h₁ : c.1 ≠ d.1 → v.1 ≠ 0) (h₂ : c.2 ≠ d.2 → v.2 ≠ 0) :
    chart2 d (trans2 c d v) = chart2 c v := by
  simp [chart2, trans2, chart1_trans c.1 d.1 v.1 h₁, chart1_trans c.2 d.2 v.2 h₂]

omit [DecidableEq K] in
/-- The cocycle condition on `P¹`: `c → d → e` equals `c → e`. -/
theorem trans1_cocycle (c d e : Bool) (x : K) : trans1 d e (trans1 c d x) = trans1 c e x := by
  cases c <;> cases d <;> cases e <;> simp [trans1]

omit [DecidableEq K] in
/-- The cocycle condition: `c → d → e` equals `c → e`, so the gluing is consistent on every
triple overlap. -/
theorem trans2_cocycle (c d e : Bool × Bool) (v : K × K) :
    trans2 d e (trans2 c d v) = trans2 c e v := by
  simp [trans2, trans1_cocycle]

/-! ### The count -/

/-- With `a = q - 1` points on the circle: one torus (`a ^ 2`), four circles (`a` each), four
fixed points. In ordinary coordinates it is `(q + 1) ^ 2`. -/
theorem fan_count {R : Type*} [CommSemiring R] (a : R) : a ^ 2 + 4 * a + 4 = (a + 2) ^ 2 := by
  ring

omit [Field K] [DecidableEq K] in
/-- Over a finite field with `q` elements, `P¹ × P¹` has `(q + 1) ^ 2` points. -/
theorem card_P1_prod_P1 [Fintype K] :
    Fintype.card (P1 K × P1 K) = (Fintype.card K + 1) ^ 2 := by
  rw [Fintype.card_prod, sq]
  congr 1 <;> exact Fintype.card_option

end MultAdd.Toric
