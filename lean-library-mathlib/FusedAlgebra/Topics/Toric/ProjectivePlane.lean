module

public import FusedAlgebra.Basic

/-!
# The projective plane `P²` from three toric charts

Gluing at the level of points, over a field `K`, as in `FusedAlgebra.Topics.Toric.ToricCharts`.

* The fan of `P²` has three rays, `e₁`, `e₂` and `-e₁ - e₂`, and three cones (`Cone`). The cones
  are not all quadrants, so `P²` is not a product and its charts do not just invert one
  coordinate at a time.
* Every chart change (`trans`) is a monomial map. `c0 → c1` is `(a, b) ↦ (a⁻¹, a⁻¹ * b)`, which
  mixes the two coordinates. In mult-add coordinates `a = 1 + u`, `b = 1 + w` it is the mult-add
  inverse of `u`, then a mult-add onto `w` (`trans01_multAdd`, `trans02_multAdd`): the gluing of
  `P²` is mult-add division.
* The charts are injective and cover (`chart_injective`, `chart_cover`), the transitions carry
  one chart to the other on the overlap (`chart_trans`, `inChart_iff`) and satisfy the cocycle
  condition (`trans_cocycle`).
* Windings `x ↦ x ^ n` commute with every chart change (`trans_wind`), so
  `[x : y : z] ↦ [x ^ n : y ^ n : z ^ n]` is one well-defined map `wind n` of the plane
  (`chart_wind`).
* Points: the torus lies in all three charts, each fixed point in one. With `a = q - 1` the count
  is `a ^ 2 + 3 * a + 3 = q ^ 2 + q + 1` (`fan_count`, `card_plane`).

## Relation to mathlib

The point set here is `Plane K = (K × K) ⊕ Option K`, the cell decomposition
`𝔸² ⊔ 𝔸¹ ⊔ pt` with `inl (a, b) = [a : b : 1]`, `inr (some u) = [1 : u : 0]` and
`inr none = [0 : 1 : 0]`. Mathlib's point-level projective plane is
`Projectivization K (Fin 3 → K)`; the map just described is `toProjectivization`, and through it
the three charts are the standard affine charts `(a, b) ↦ [a : b : 1]`, `(v, u) ↦ [1 : u : v]`
and `(t, s) ↦ [s : 1 : t]` (`toProjectivization_chart_c0`, `_c1`, `_c2`). The map is a
bijection, but that is not proved here (see the port notes). This file does not build `P²` as a
scheme; mathlib's scheme-level object is `AlgebraicGeometry.ProjectiveSpace`.
-/

@[expose] public section

namespace MultAdd.Toric.P2

variable {K : Type*} [Field K] [DecidableEq K]

/-- A point `[x : y : z]` of `P²`. `inl (a, b)` is `[a : b : 1]`; `inr (some u)` is
`[1 : u : 0]`; `inr none` is `[0 : 1 : 0]`. -/
abbrev Plane (K : Type*) := (K × K) ⊕ Option K

/-- The three cones of the fan: `c0 = ⟨e₁, e₂⟩`, `c1 = ⟨e₂, -e₁ - e₂⟩`, `c2 = ⟨-e₁ - e₂, e₁⟩`. -/
inductive Cone | c0 | c1 | c2
  deriving DecidableEq, Repr

open Cone

/-- The three charts, one per cone of the fan of `P²`. `c0` sees `(x / z, y / z)`, `c1` sees
`(z / x, y / x)`, `c2` sees `(z / y, x / y)`. -/
def chart : Cone → K × K → Plane K
  | c0, (a, b) => .inl (a, b)
  | c1, (v, u) => if v = 0 then .inr (some u) else .inl (v⁻¹, u * v⁻¹)
  | c2, (t, s) => if t = 0 then (if s = 0 then .inr none else .inr (some s⁻¹))
                  else .inl (s * t⁻¹, t⁻¹)

/-- Each chart is injective. -/
theorem chart_injective (c : Cone) : Function.Injective (chart (K := K) c) := by
  rintro ⟨x, y⟩ ⟨x', y'⟩ h
  cases c <;> simp only [chart] at h
  · simpa using h
  · by_cases hx : x = 0 <;> by_cases hx' : x' = 0
    · simpa [hx, hx'] using h
    · simp [hx, hx'] at h
    · simp [hx, hx'] at h
    · obtain ⟨rfl, h₂⟩ : x = x' ∧ y * x⁻¹ = y' * x'⁻¹ := by simpa [hx, hx'] using h
      rw [mul_right_cancel₀ (inv_ne_zero hx) h₂]
  · by_cases hx : x = 0 <;> by_cases hx' : x' = 0
    · by_cases hy : y = 0 <;> by_cases hy' : y' = 0 <;> simp_all
    · by_cases hy : y = 0 <;> simp [hx, hx', hy] at h
    · by_cases hy' : y' = 0 <;> simp [hx, hx', hy'] at h
    · obtain ⟨h₂, rfl⟩ : y * x⁻¹ = y' * x'⁻¹ ∧ x = x' := by simpa [hx, hx'] using h
      rw [mul_right_cancel₀ (inv_ne_zero hx) h₂]

/-- The three charts cover `P²`. -/
theorem chart_cover (p : Plane K) : ∃ c v, chart c v = p := by
  rcases p with ⟨a, b⟩ | _ | u
  · exact ⟨c0, (a, b), rfl⟩
  · exact ⟨c2, (0, 0), by simp [chart]⟩
  · exact ⟨c1, (0, u), by simp [chart]⟩

/-- A point lies in chart `c`. -/
def InChart (c : Cone) (p : Plane K) : Prop := ∃ v, chart c v = p

/-- Transition maps. Every one is a monomial map: a matrix in `GL₂(ℤ)` acting on exponents, so
unlike `P¹ × P¹` a chart change can mix the two coordinates. -/
def trans : Cone → Cone → K × K → K × K
  | c0, c1, (a, b) => (a⁻¹, a⁻¹ * b)
  | c0, c2, (a, b) => (b⁻¹, a * b⁻¹)
  | c1, c0, (v, u) => (v⁻¹, u * v⁻¹)
  | c1, c2, (v, u) => (v * u⁻¹, u⁻¹)
  | c2, c0, (t, s) => (s * t⁻¹, t⁻¹)
  | c2, c1, (t, s) => (t * s⁻¹, s⁻¹)
  | _, _, v => v

/-- Where chart `c`'s point `v` also lies in chart `d`: the coordinate the transition divides by
is nonzero. -/
def Dom : Cone → Cone → K × K → Prop
  | c0, c1, (a, _) => a ≠ 0
  | c0, c2, (_, b) => b ≠ 0
  | c1, c0, (v, _) => v ≠ 0
  | c1, c2, (_, u) => u ≠ 0
  | c2, c0, (t, _) => t ≠ 0
  | c2, c1, (_, s) => s ≠ 0
  | _, _, _ => True

/-- Where two charts overlap, the transition carries one chart's coordinates to the other's. -/
theorem chart_trans (c d : Cone) (v : K × K) (h : Dom c d v) :
    chart d (trans c d v) = chart c v := by
  rcases v with ⟨x, y⟩
  cases c <;> cases d <;> simp only [Dom, trans, chart] at h ⊢
  all_goals split_ifs <;> simp_all
  all_goals field_simp

/-- A point of chart `c` lies in chart `d` exactly on the transition's domain. -/
theorem inChart_iff (c d : Cone) (v : K × K) : InChart d (chart c v) ↔ Dom c d v := by
  constructor
  · rintro ⟨⟨x', y'⟩, hw⟩
    rcases v with ⟨x, y⟩
    cases c <;> cases d <;> simp only [Dom, chart] at hw ⊢
    all_goals rintro rfl; split_ifs at hw <;> simp_all
  · intro h; exact ⟨_, chart_trans c d v h⟩

omit [DecidableEq K] in
/-- The cocycle condition, from injectivity: `c → d → e` equals `c → e` on the triple overlap. -/
theorem trans_cocycle (c d e : Cone) (v : K × K) (hd : Dom c d v) (he : Dom c e v) :
    trans d e (trans c d v) = trans c e v := by
  classical
  have hde : Dom d e (trans c d v) := by
    rw [← inChart_iff, chart_trans c d v hd, inChart_iff]; exact he
  apply chart_injective e
  rw [chart_trans d e _ hde, chart_trans c d v hd, chart_trans c e v he]

/-- Each chart's origin is one of the three fixed points `[0 : 0 : 1]`, `[1 : 0 : 0]`,
`[0 : 1 : 0]`, and it lies in that chart only. -/
theorem origin_in_one_chart (c d : Cone) : InChart d (chart c ((0 : K), (0 : K))) ↔ d = c := by
  rw [inChart_iff]; cases c <;> cases d <;> simp [Dom]

/-- Points of the torus (both coordinates nonzero) lie in all three charts. -/
theorem torus_in_all (c d : Cone) (x y : K) (hx : x ≠ 0) (hy : y ≠ 0) :
    InChart d (chart c (x, y)) := by
  rw [inChart_iff]; cases c <;> cases d <;> simp [Dom, hx, hy]

/-! ### Mult-add coordinates: `a = 1 + u`, `b = 1 + w` -/

omit [DecidableEq K] in
/-- In mult-add coordinates, the chart change `c0 → c1`, `(a, b) ↦ (a⁻¹, a⁻¹ * b)`, is the
mult-add inverse of `u`, then a mult-add onto `w`. -/
theorem trans01_multAdd (u w : K) (h : 1 + u ≠ 0) :
    trans c0 c1 (1 + u, 1 + w) = (1 + multAddInv u, 1 + (w ✳ multAddInv u)) := by
  rw [one_add_multAdd, ← inv_one_add u h, mul_comm]; rfl

omit [DecidableEq K] in
theorem trans02_multAdd (u w : K) (h : 1 + w ≠ 0) :
    trans c0 c2 (1 + u, 1 + w) = (1 + multAddInv w, 1 + (u ✳ multAddInv w)) := by
  rw [one_add_multAdd, ← inv_one_add w h]; rfl

/-! ### Windings `x ↦ x ^ n`

The winding on chart coordinates is the pointwise power `v ^ n` in the monoid `K × K`. -/

/-- The winding `[x : y : z] ↦ [x ^ n : y ^ n : z ^ n]` on `P²` itself. -/
def wind (n : ℕ) : Plane K → Plane K
  | .inl (a, b) => .inl (a ^ n, b ^ n)
  | .inr (some u) => .inr (some (u ^ n))
  | .inr none => .inr none

omit [DecidableEq K] in
/-- Windings commute with every chart change: monomial maps commute with `x ↦ x ^ n`. -/
theorem trans_wind (c d : Cone) (n : ℕ) (v : K × K) : trans c d (v ^ n) = trans c d v ^ n := by
  rcases v with ⟨x, y⟩
  cases c <;> cases d <;> simp [trans, mul_pow]

/-- So the winding is one well-defined map of `P²`, read the same way in every chart. -/
theorem chart_wind (c : Cone) {n : ℕ} (hn : n ≠ 0) (v : K × K) :
    chart c (v ^ n) = wind n (chart c v) := by
  rcases v with ⟨x, y⟩
  cases c
  · rfl
  · by_cases hx : x = 0 <;> simp [chart, wind, hx, hn, mul_pow]
  · by_cases hx : x = 0 <;> by_cases hy : y = 0 <;> simp [chart, wind, hx, hy, hn, mul_pow]

/-! ### Counting -/

/-- With `a = q - 1`: one torus (`a ^ 2`), three circles (`a` each), three fixed points; in
ordinary coordinates `q ^ 2 + q + 1`. -/
theorem fan_count {R : Type*} [CommSemiring R] (a : R) :
    a ^ 2 + 3 * a + 3 = (a + 1) ^ 2 + (a + 1) + 1 := by
  ring

omit [Field K] [DecidableEq K] in
/-- Over a finite field with `q` elements, `P²` has `q ^ 2 + q + 1` points. -/
theorem card_plane [Fintype K] :
    Fintype.card (Plane K) = Fintype.card K ^ 2 + Fintype.card K + 1 := by
  simp [sq, add_assoc]

/-! ### Homogeneous coordinates

The map to mathlib's `Projectivization K (Fin 3 → K)`, and the three charts read through it:
they are the standard affine charts `z ≠ 0`, `x ≠ 0`, `y ≠ 0`. -/

/-- The point of mathlib's projective plane with the given homogeneous coordinates. -/
def toProjectivization : Plane K → Projectivization K (Fin 3 → K)
  | .inl (a, b) => Projectivization.mk K ![a, b, 1] (by simp)
  | .inr (some u) => Projectivization.mk K ![1, u, 0] (by simp)
  | .inr none => Projectivization.mk K ![0, 1, 0] (by simp)

/-- Chart `c0` is `(a, b) ↦ [a : b : 1]`. -/
theorem toProjectivization_chart_c0 (a b : K) :
    toProjectivization (chart c0 (a, b)) = Projectivization.mk K ![a, b, 1] (by simp) :=
  rfl

/-- Chart `c1` is `(v, u) ↦ [1 : u : v]`. -/
theorem toProjectivization_chart_c1 (v u : K) :
    toProjectivization (chart c1 (v, u)) = Projectivization.mk K ![1, u, v] (by simp) := by
  by_cases hv : v = 0
  · simp [chart, toProjectivization, hv]
  · simp only [chart, hv, ite_false, toProjectivization, Projectivization.mk_eq_mk_iff']
    exact ⟨v⁻¹, by ext i; fin_cases i <;> simp [hv, mul_comm]⟩

/-- Chart `c2` is `(t, s) ↦ [s : 1 : t]`. -/
theorem toProjectivization_chart_c2 (t s : K) :
    toProjectivization (chart c2 (t, s)) = Projectivization.mk K ![s, 1, t] (by simp) := by
  by_cases ht : t = 0
  · by_cases hs : s = 0
    · simp [chart, toProjectivization, ht, hs]
    · simp only [chart, ht, hs, ite_true, ite_false, toProjectivization,
        Projectivization.mk_eq_mk_iff']
      exact ⟨s⁻¹, by ext i; fin_cases i <;> simp [hs]⟩
  · simp only [chart, ht, ite_false, toProjectivization, Projectivization.mk_eq_mk_iff']
    exact ⟨t⁻¹, by ext i; fin_cases i <;> simp [ht, mul_comm]⟩

end MultAdd.Toric.P2
