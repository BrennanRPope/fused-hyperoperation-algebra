/-
  The projective plane P² from three toric charts
  ================================================

  Part of Fused Mult-Add Algebra. Companion to ToricCharts.lean (P¹ × P¹,
  same folder). Python check: p2_check.py.

  STATUS: CHECKED. Core Lean 4 only (no mathlib), compiled with Lean
  4.37.0-rc1 with no errors and no `sorry`; `#print axioms` shows only the
  standard three (propext, Classical.choice, Quot.sound). Run it with
  `lean ProjectivePlane.lean`. It works over any field K with Lean's
  built-in field structure (`Lean.Grind.Field`).

  What is new compared with P¹ × P¹:

  * The fan of P² has three rays, e₁, e₂ and -e₁-e₂, and three cones. The
    cones are not all quadrants, so P² is not a product and its charts do
    not just invert one coordinate at a time.

  * Every chart change is a monomial map: a 2×2 integer matrix acting on
    exponents. c0 → c1 is (a, b) ↦ (1/a, b/a), which mixes the two
    coordinates. In mult-add coordinates (a = 1 + u, b = 1 + w) that is
    (u, w) ↦ (u⁻, w ◯ u⁻): take the mult-add inverse of u, then mult-add
    it onto w. So the gluing of P² is mult-add division.

  * Windings x ↦ xⁿ commute with every chart change, because monomial maps
    commute with powers. So the winding [x:y:z] ↦ [xⁿ:yⁿ:zⁿ] is one
    well-defined map of the whole plane, read the same way in every chart.
    For n = q on a field with q elements it is the Frobenius, and its fixed
    points are the q² + q + 1 points of P².

  * Points: the torus lies in all three charts, each circle (a coordinate
    line minus its two fixed points) in two, each fixed point in one. With
    a = q - 1 that is a² + 3a + 3 = q² + q + 1.

  As with ToricCharts.lean this is gluing at the level of points (charts as
  injective maps that cover, transition maps, the cocycle condition), not a
  construction of P² as a scheme, which would need mathlib.

  LIBRARY NOTE: mult-add and its inverse now come from FusedAlgebra.Base.
-/
import FusedAlgebra.Base

namespace CircleOp.Toric.P2

open Lean.Grind

set_option linter.unusedSectionVars false

variable {K : Type} [Field K] [DecidableEq K]

/-- A point [x:y:z] of P². `inl (a, b)` is [a:b:1]; `inr (some u)` is [1:u:0];
`inr none` is [0:1:0]. -/
abbrev P2 (K : Type) := Sum (K × K) (Option K)

/-- The three cones of the fan: c0 = ⟨e₁, e₂⟩, c1 = ⟨e₂, -e₁-e₂⟩, c2 = ⟨-e₁-e₂, e₁⟩. -/
inductive Cone | c0 | c1 | c2
  deriving DecidableEq, Repr

open Cone

/-- The three charts, one per cone of the fan of P² (rays e₁, e₂, -e₁-e₂).
c0 sees (x/z, y/z), c1 sees (z/x, y/x), c2 sees (z/y, x/y). -/
def chart : Cone → K × K → P2 K
  | c0, (a, b) => .inl (a, b)
  | c1, (v, u) => if v = 0 then .inr (some u) else .inl (v⁻¹, u * v⁻¹)
  | c2, (t, s) => if t = 0 then (if s = 0 then .inr none else .inr (some s⁻¹))
                  else .inl (s * t⁻¹, t⁻¹)

theorem inv_inv' (x : K) : (x⁻¹)⁻¹ = x := by grind
theorem inv_ne (x : K) (h : x ≠ 0) : x⁻¹ ≠ 0 := by grind

/-- Each chart is injective. -/
theorem chart_injective (c : Cone) (v w : K × K) (h : chart c v = chart c w) : v = w := by
  rcases v with ⟨x, y⟩; rcases w with ⟨x', y'⟩
  cases c <;> simp only [chart] at h
  · simpa using h
  · split at h <;> split at h <;> simp at h <;> grind
  · split at h <;> split at h
    all_goals (try split at h) <;> (try split at h) <;> simp at h
    all_goals grind

/-- The three charts cover P². -/
theorem chart_cover (p : P2 K) : ∃ c v, chart c v = p := by
  rcases p with ⟨a, b⟩ | _ | u
  · exact ⟨c0, (a, b), rfl⟩
  · exact ⟨c2, (0, 0), by simp [chart]⟩
  · exact ⟨c1, (0, u), by simp [chart]⟩

/-- A point lies in chart c. -/
def InChart (c : Cone) (p : P2 K) : Prop := ∃ v, chart c v = p

/-- Transition maps. Every one is a monomial map: a matrix in GL₂(ℤ) acting on
exponents, so unlike P¹ × P¹ a chart change can mix the two coordinates. -/
def trans : Cone → Cone → K × K → K × K
  | c0, c1, (a, b) => (a⁻¹, a⁻¹ * b)
  | c0, c2, (a, b) => (b⁻¹, a * b⁻¹)
  | c1, c0, (v, u) => (v⁻¹, u * v⁻¹)
  | c1, c2, (v, u) => (v * u⁻¹, u⁻¹)
  | c2, c0, (t, s) => (s * t⁻¹, t⁻¹)
  | c2, c1, (t, s) => (t * s⁻¹, s⁻¹)
  | _, _, v => v

/-- Where chart c's point v also lies in chart d: the coordinate the
transition divides by is nonzero. -/
def Dom : Cone → Cone → K × K → Prop
  | c0, c1, (a, _) => a ≠ 0
  | c0, c2, (_, b) => b ≠ 0
  | c1, c0, (v, _) => v ≠ 0
  | c1, c2, (_, u) => u ≠ 0
  | c2, c0, (t, _) => t ≠ 0
  | c2, c1, (_, s) => s ≠ 0
  | _, _, _ => True

theorem chart_trans (c d : Cone) (v : K × K) (h : Dom c d v) :
    chart d (trans c d v) = chart c v := by
  rcases v with ⟨x, y⟩
  cases c <;> cases d <;> simp only [Dom, trans, chart] at h ⊢
  all_goals (repeat' split) <;> simp_all <;> grind

/-- A point of chart c lies in chart d exactly on the transition's domain. -/
theorem inChart_iff (c d : Cone) (v : K × K) : InChart d (chart c v) ↔ Dom c d v := by
  constructor
  · rintro ⟨⟨x', y'⟩, hw⟩
    rcases v with ⟨x, y⟩
    cases c <;> cases d <;> simp only [Dom, chart] at hw ⊢
    all_goals (repeat' split at hw) <;> simp_all <;> grind
  · intro h; exact ⟨_, chart_trans c d v h⟩

/-- The cocycle condition, from injectivity: c → d → e equals c → e on the
triple overlap. -/
theorem trans_cocycle (c d e : Cone) (v : K × K) (hd : Dom c d v) (he : Dom c e v) :
    trans d e (trans c d v) = trans c e v := by
  have hde : Dom d e (trans c d v) := by
    rw [← inChart_iff, chart_trans c d v hd, inChart_iff]; exact he
  apply chart_injective e
  rw [chart_trans d e _ hde, chart_trans c d v hd, chart_trans c e v he]

/-- Each chart's origin is one of the three fixed points [0:0:1], [1:0:0],
[0:1:0], and it lies in that chart only. -/
theorem origin_in_one_chart (c d : Cone) : InChart d (chart c ((0 : K), (0 : K))) ↔ d = c := by
  rw [inChart_iff]; cases c <;> cases d <;> simp [Dom]

/-- Points of the torus (both coordinates nonzero) lie in all three charts. -/
theorem torus_in_all (c d : Cone) (x y : K) (hx : x ≠ 0) (hy : y ≠ 0) :
    InChart d (chart c (x, y)) := by
  rw [inChart_iff]; cases c <;> cases d <;> simp [Dom, hx, hy]

/-! Mult-add coordinates: a = 1 + u, b = 1 + w. -/

/-- In mult-add coordinates, the chart change c0 → c1, (a, b) ↦ (1/a, b/a),
is (u, w) ↦ (u⁻, w ◯ u⁻): the mult-add inverse, then a mult-add. -/
theorem trans01_circ (u w : K) (h : 1 + u ≠ 0) :
    trans c0 c1 (1 + u, 1 + w) = (1 + circInv u, 1 + circ w (circInv u)) := by
  simp only [trans, circ, circInv, Prod.mk.injEq]; constructor <;> grind

theorem trans02_circ (u w : K) (h : 1 + w ≠ 0) :
    trans c0 c2 (1 + u, 1 + w) = (1 + circInv w, 1 + circ u (circInv w)) := by
  simp only [trans, circ, circInv, Prod.mk.injEq]; constructor <;> grind

/-! Windings x ↦ xⁿ. -/

theorem inv_pow' (x : K) (n : Nat) : (x⁻¹) ^ n = (x ^ n)⁻¹ := by
  induction n with
  | zero => grind
  | succ n ih => rw [Semiring.pow_succ, Semiring.pow_succ, ih]; grind

theorem pow_succ_eq_zero (x : K) (n : Nat) : x ^ (n + 1) = 0 ↔ x = 0 := by
  induction n with
  | zero => grind
  | succ n ih => rw [Semiring.pow_succ]; grind

/-- The winding by n on chart coordinates. -/
def pw (n : Nat) (v : K × K) : K × K := (v.1 ^ n, v.2 ^ n)

/-- The winding [x:y:z] ↦ [xⁿ:yⁿ:zⁿ] on P² itself. -/
def wind (n : Nat) : P2 K → P2 K
  | .inl (a, b) => .inl (a ^ n, b ^ n)
  | .inr (some u) => .inr (some (u ^ n))
  | .inr none => .inr none

/-- Windings commute with every chart change: monomial maps commute with
x ↦ xⁿ. -/
theorem trans_wind (c d : Cone) (n : Nat) (v : K × K) :
    trans c d (pw n v) = pw n (trans c d v) := by
  rcases v with ⟨x, y⟩
  cases c <;> cases d <;> simp [trans, pw, inv_pow', mul_pow']

/-- So the winding is one well-defined map of P², read the same way in
every chart. -/
theorem chart_wind (c : Cone) (n : Nat) (v : K × K) :
    chart c (pw (n + 1) v) = wind (n + 1) (chart c v) := by
  rcases v with ⟨x, y⟩
  cases c
  · rfl
  · by_cases hx : x = 0 <;> simp [chart, pw, wind, hx, pow_succ_eq_zero, inv_pow', mul_pow']
  · by_cases hx : x = 0 <;> by_cases hy : y = 0 <;>
      simp [chart, pw, wind, hx, hy, pow_succ_eq_zero, inv_pow', mul_pow']

/-! Counting with a = q - 1: one torus (a²), three circles (a each), three
fixed points; in ordinary coordinates q² + q + 1. -/
theorem fan_count (a : Nat) : a ^ 2 + 3 * a + 3 = (a + 1) ^ 2 + (a + 1) + 1 := by grind

end CircleOp.Toric.P2
