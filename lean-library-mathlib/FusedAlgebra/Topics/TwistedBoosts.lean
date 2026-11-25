module

public import FusedAlgebra.Topics.Relativity

/-!
# Twisted boosts: two boosts in different directions

A velocity in the plane is a complex number `u` (a half-velocity on the disk). The dial `q : ℝ`
is the number in the denominator `1 - q * a * b` of the one-dimensional law: `q = -1` is
relativity (`c = 1`), `q = 0` is Galilean, `q = 1` is rotations (tangent addition).

## Main definitions

* `TwistedBoosts.dot`, `TwistedBoosts.cross`: the dot and cross parts of the corner `conj u * v`.
* `TwistedBoosts.boost`: the dial's boost `!![1, u; -q * conj u, 1]`.
* `TwistedBoosts.tX`, `TwistedBoosts.tY`, `TwistedBoosts.den`: the diagonal entries
  `X ± i * Y` of a product of two boosts, and the denominator
  `den = X - i * Y = 1 - q * conj u * v`.
* `TwistedBoosts.oplus`, `TwistedBoosts.twist`: the fused sum `(u + v) / den` and the twist
  `conj den / den`.
* `TwistedBoosts.shape`, `TwistedBoosts.shapeSubmonoid`, `TwistedBoosts.sdet`: the matrices
  `!![α, β; -q * conj β, conj α]`, which are closed under products, and their determinant.

## Main results

* `corner_split`: the corner `conj u * v` splits into a dot part and a cross part.
* `boost_mul`: the product of two boosts.
* `twist_trivial_iff`: there is no twist exactly when `q * cross u v = 0`.
* `order_swap`: order matters only through the twist, `oplus q u v = twist q u v * oplus q v u`.
* `shape_mul`, `det_shape`, `sdet_mul`: one shape per dial, with a determinant that multiplies.
* `light_wall_2d`, `under_light_2d`: at `q = -1`, below light stays below light in every
  direction.

For `q = -1` this is the action of `SU(1,1)` on the disk (Möbius addition of velocities and the
Thomas–Wigner rotation), for `q = 1` it is `SU(2)`; both are classical. The core-Lean original
modelled complex numbers as pairs of rationals with hand-written operations and 2×2 matrices as
a four-field structure; here they are mathlib's `ℂ` and `Matrix (Fin 2) (Fin 2) ℂ`, with a real
dial `q`. Every original statement is the special case of rational coordinates.
-/

@[expose] public section

namespace TwistedBoosts

open MultAdd ThermoFused Relativity Matrix Complex ComplexConjugate

/-! ### The corner splits into a dot part and a cross part -/

/-- The dot product of two plane vectors, seen as complex numbers. -/
def dot (u v : ℂ) : ℝ := u.re * v.re + u.im * v.im

/-- The cross product of two plane vectors, seen as complex numbers. -/
def cross (u v : ℂ) : ℝ := u.re * v.im - u.im * v.re

/-- **The corner splits**: `conj u * v = dot u v + i * cross u v`. -/
theorem corner_split (u v : ℂ) : conj u * v = ⟨dot u v, cross u v⟩ := by
  apply Complex.ext
  · simp [dot]
  · simp [cross, sub_eq_add_neg]

/-- The other order is the conjugate. -/
theorem corner_swap (u v : ℂ) : conj v * u = conj (conj u * v) := by
  rw [map_mul, conj_conj, mul_comm]

theorem cross_anti (u v : ℂ) : cross v u = -cross u v := by
  unfold cross; ring

theorem dot_comm (u v : ℂ) : dot v u = dot u v := by
  unfold dot; ring

/-- Twice the signed area of the triangle `p`, `r`, `s` (the shoelace formula). -/
def area2 (p r s : ℂ) : ℝ := (r.re - p.re) * (s.im - p.im) - (r.im - p.im) * (s.re - p.re)

/-- The cross part is twice the signed area of the triangle `0`, `u`, `v`. -/
theorem cross_is_area (u v : ℂ) : cross u v = area2 0 u v := by
  simp [cross, area2]

/-- Collinear speeds (both on the real axis) have no cross part. -/
theorem cross_collinear (a c : ℝ) : cross a c = 0 := by
  simp [cross]

/-! ### The boost of the dial -/

/-- The dial's boost `!![1, u; -q * conj u, 1]`. -/
def boost (q : ℝ) (u : ℂ) : Matrix (Fin 2) (Fin 2) ℂ := !![1, u; -q * conj u, 1]

/-- At `q = -1` on the real axis the boost is `!![1, u; u, 1]`. -/
theorem boost_real (u : ℝ) : boost (-1) u = !![1, (u : ℂ); u, 1] := by
  simp [boost]

/-- At `q = -1` on the real axis the boost is the boost of `Relativity`, read in `ℂ`. -/
theorem boost_real_eq_map (u : ℝ) : boost (-1) u = (Relativity.boost u).map (↑) := by
  rw [boost_real]
  ext i j
  fin_cases i <;> fin_cases j <;> simp [Relativity.boost]

/-- The real part `X` of the entry that carries the twist. -/
def tX (q : ℝ) (u v : ℂ) : ℝ := 1 - q * dot u v

/-- The imaginary part `Y` of the entry that carries the twist. -/
def tY (q : ℝ) (u v : ℂ) : ℝ := q * cross u v

/-! ### The product of two boosts -/

/-- The top-left entry of a product of boosts is `1 - q * u * conj v = X + i * Y`. -/
theorem entry_tl (q : ℝ) (u v : ℂ) : 1 - q * (u * conj v) = ⟨tX q u v, tY q u v⟩ := by
  apply Complex.ext <;>
    simp only [tX, tY, dot, cross, sub_re, sub_im, one_re, one_im, mul_re, mul_im, ofReal_re,
      ofReal_im, conj_re, conj_im] <;>
    ring

/-- The bottom-right entry of a product of boosts is `1 - q * conj u * v = X - i * Y`. -/
theorem entry_br (q : ℝ) (u v : ℂ) : 1 - q * (conj u * v) = ⟨tX q u v, -tY q u v⟩ := by
  apply Complex.ext <;>
    simp only [tX, tY, dot, cross, sub_re, sub_im, one_re, one_im, mul_re, mul_im, ofReal_re,
      ofReal_im, conj_re, conj_im] <;>
    ring

/-- **Two boosts**: `B u * B v = !![X + i * Y, u + v; -q * conj (u + v), X - i * Y]`. -/
theorem boost_mul (q : ℝ) (u v : ℂ) :
    boost q u * boost q v =
      !![⟨tX q u v, tY q u v⟩, u + v; -q * conj (u + v), ⟨tX q u v, -tY q u v⟩] := by
  rw [← entry_tl, ← entry_br, boost, boost, mul_fin_two]
  ext i j
  fin_cases i <;> fin_cases j <;> simp <;> ring

/-- The diagonal entries are conjugates. -/
theorem diag_conj (q : ℝ) (u v : ℂ) :
    (boost q u * boost q v) 0 0 = conj ((boost q u * boost q v) 1 1) := by
  rw [boost_mul]
  apply Complex.ext <;> simp

/-- In one dimension at `q = -1` the product is `(1 + u * v) • B (vel u v)`, as in
`Relativity.boost_mul`. -/
theorem boost_mul_real (u v : ℝ) :
    boost (-1) u * boost (-1) v =
      (!![1 + u * v, u + v; u + v, 1 + u * v] : Matrix (Fin 2) (Fin 2) ℂ) := by
  rw [boost_mul, ← entry_tl, ← entry_br]
  ext i j
  fin_cases i <;> fin_cases j <;> simp

/-! ### The twist -/

/-- **No twist exactly when the cross part is switched off**: `X + i * Y = X - i * Y` if and only
if `q * cross u v = 0`. -/
theorem twist_trivial_iff (q : ℝ) (u v : ℂ) :
    (⟨tX q u v, tY q u v⟩ : ℂ) = ⟨tX q u v, -tY q u v⟩ ↔ q * cross u v = 0 := by
  simp [Complex.ext_iff, tY, self_eq_neg]

/-- Collinear boosts never twist. -/
theorem collinear_no_twist (q a c : ℝ) : tY q a c = 0 := by
  simp [tY, cross]

/-- The Galilean dial `q = 0` never twists, in any directions. -/
theorem galilean_no_twist (u v : ℂ) : tY 0 u v = 0 ∧ tX 0 u v = 1 := by
  simp [tY, tX]

/-- The cosine of the twist angle, as a rational function of `X` and `Y`. -/
noncomputable def cosT (q : ℝ) (u v : ℂ) : ℝ :=
  (tX q u v * tX q u v - tY q u v * tY q u v) / (tX q u v * tX q u v + tY q u v * tY q u v)

/-- The sine of the twist angle, as a rational function of `X` and `Y`. -/
noncomputable def sinT (q : ℝ) (u v : ℂ) : ℝ :=
  2 * tX q u v * tY q u v / (tX q u v * tX q u v + tY q u v * tY q u v)

/-- The twist lies on the unit circle. -/
theorem twist_on_circle (q : ℝ) (u v : ℂ)
    (h : tX q u v * tX q u v + tY q u v * tY q u v ≠ 0) :
    cosT q u v * cosT q u v + sinT q u v * sinT q u v = 1 := by
  have h' : tX q u v ^ 2 + tY q u v ^ 2 ≠ 0 := by rwa [sq, sq]
  unfold cosT sinT
  field_simp; ring

/-- The twist turns `X - i * Y` into `X + i * Y`. -/
theorem twist_turns (q : ℝ) (u v : ℂ) (h : tX q u v * tX q u v + tY q u v * tY q u v ≠ 0) :
    (⟨cosT q u v, sinT q u v⟩ : ℂ) * ⟨tX q u v, -tY q u v⟩ = ⟨tX q u v, tY q u v⟩ := by
  have h' : tX q u v ^ 2 + tY q u v ^ 2 ≠ 0 := by rwa [sq, sq]
  apply Complex.ext <;> simp only [mul_re, mul_im, cosT, sinT] <;> field_simp <;> ring

/-- **The half-angle of the twist**: `tan (ω / 2) = sin ω / (1 + cos ω) = Y / X`. -/
theorem twist_half_angle (q : ℝ) (u v : ℂ) (hX : tX q u v ≠ 0) :
    sinT q u v / (1 + cosT q u v) = tY q u v / tX q u v := by
  have h' : tX q u v ^ 2 + tY q u v ^ 2 ≠ 0 := by positivity
  have e : 1 + cosT q u v =
      2 * tX q u v * tX q u v / (tX q u v * tX q u v + tY q u v * tY q u v) := by
    unfold cosT; field_simp; ring
  rw [e, sinT]
  field_simp

/-! ### Order matters only through the twist -/

/-- The denominator of the fused sum, `X - i * Y = 1 - q * conj u * v`. -/
def den (q : ℝ) (u v : ℂ) : ℂ := ⟨tX q u v, -tY q u v⟩

theorem den_eq (q : ℝ) (u v : ℂ) : den q u v = 1 - q * (conj u * v) :=
  (entry_br q u v).symm

/-- The other order has the conjugate denominator. -/
theorem den_swap (q : ℝ) (u v : ℂ) : den q v u = conj (den q u v) := by
  apply Complex.ext <;> simp only [den, tX, tY, dot, cross, conj_re, conj_im] <;> ring

theorem normSq_den (q : ℝ) (u v : ℂ) :
    normSq (den q u v) = tX q u v * tX q u v + tY q u v * tY q u v := by
  simp [den]

theorem den_ne_zero {q : ℝ} {u v : ℂ} (h : tX q u v * tX q u v + tY q u v * tY q u v ≠ 0) :
    den q u v ≠ 0 := by
  rwa [← normSq_den, Ne, normSq_eq_zero] at h

/-- The fused sum of two plane velocities, `(u + v) / (1 - q * conj u * v)`. -/
noncomputable def oplus (q : ℝ) (u v : ℂ) : ℂ := (u + v) / den q u v

/-- The fused sum written without complex division, as in the original definition. -/
theorem oplus_eq_smul (q : ℝ) (u v : ℂ) :
    oplus q u v = (1 / normSq (den q u v)) • ((u + v) * conj (den q u v)) := by
  rw [oplus, div_eq_mul_inv, Complex.inv_def, real_smul, one_div]; ring

/-- The twist `conj D / D` of two boosts, where `D = den q u v`. -/
noncomputable def twist (q : ℝ) (u v : ℂ) : ℂ := conj (den q u v) / den q u v

/-- The twist written without complex division, as in the original definition. -/
theorem twist_eq_smul (q : ℝ) (u v : ℂ) :
    twist q u v = (1 / normSq (den q u v)) • (conj (den q u v) * conj (den q u v)) := by
  rw [twist, div_eq_mul_inv, Complex.inv_def, real_smul, one_div]; ring

/-- The twist is exactly the unit-circle point `(cos ω, sin ω)`. -/
theorem twist_eq (q : ℝ) (u v : ℂ) (h : tX q u v * tX q u v + tY q u v * tY q u v ≠ 0) :
    twist q u v = ⟨cosT q u v, sinT q u v⟩ := by
  rw [twist, div_eq_iff (den_ne_zero h), den, twist_turns q u v h]
  apply Complex.ext <;> simp

/-- The twist is a unit complex number. -/
theorem norm_twist (q : ℝ) (u v : ℂ) (h : tX q u v * tX q u v + tY q u v * tY q u v ≠ 0) :
    ‖twist q u v‖ = 1 := by
  rw [twist, norm_div, RCLike.norm_conj, div_self (norm_ne_zero_iff.2 (den_ne_zero h))]

/-- **Order swap**: `oplus q u v = twist q u v * oplus q v u`. -/
theorem order_swap (q : ℝ) (u v : ℂ) (h : tX q u v * tX q u v + tY q u v * tY q u v ≠ 0) :
    oplus q u v = twist q u v * oplus q v u := by
  have hD := den_ne_zero h
  have hD' : conj (den q u v) ≠ 0 := (map_ne_zero _).2 hD
  rw [oplus, oplus, twist, den_swap q u v, add_comm v u]
  field_simp

/-! ### One shape per dial, and the light wall in every direction -/

/-- The dial's shape `!![α, β; -q * conj β, conj α]`. -/
def shape (q : ℝ) (α β : ℂ) : Matrix (Fin 2) (Fin 2) ℂ := !![α, β; -q * conj β, conj α]

/-- Boosts have the shape. -/
theorem boost_shape (q : ℝ) (u : ℂ) : boost q u = shape q 1 u := by
  simp [boost, shape]

/-- **The shape is closed under products**: one group for each `q` (`q = -1` the Lorentz spin
group `SU(1,1)`, `q = 0` the Galilean one, `q = 1` the group `SU(2)`), up to scaling. -/
theorem shape_mul (q : ℝ) (α β γ δ : ℂ) :
    shape q α β * shape q γ δ = shape q (α * γ - q * (β * conj δ)) (α * δ + β * conj γ) := by
  rw [shape, shape, shape, mul_fin_two]
  ext i j
  fin_cases i <;> fin_cases j <;> simp <;> ring

/-- The matrices of the dial's shape form a submonoid of the 2×2 complex matrices. -/
def shapeSubmonoid (q : ℝ) : Submonoid (Matrix (Fin 2) (Fin 2) ℂ) where
  carrier := Set.range fun p : ℂ × ℂ ↦ shape q p.1 p.2
  one_mem' := ⟨(1, 0), by simp [shape, one_fin_two]⟩
  mul_mem' := by
    rintro _ _ ⟨⟨α, β⟩, rfl⟩ ⟨⟨γ, δ⟩, rfl⟩
    exact ⟨(_, _), (shape_mul q α β γ δ).symm⟩

theorem shape_mem_shapeSubmonoid (q : ℝ) (α β : ℂ) : shape q α β ∈ shapeSubmonoid q :=
  ⟨(α, β), rfl⟩

/-- Determinant of the shape: `|α| ^ 2 + q * |β| ^ 2`. -/
def sdet (q : ℝ) (α β : ℂ) : ℝ := normSq α + q * normSq β

/-- `sdet` is the determinant of the shape. -/
theorem det_shape (q : ℝ) (α β : ℂ) : (shape q α β).det = sdet q α β := by
  rw [shape, det_fin_two_of, sdet, ofReal_add, ofReal_mul, ← mul_conj, ← mul_conj]; ring

/-- **The determinant multiplies** (a two-square identity for each `q`). This is `Matrix.det_mul`
for the shape. -/
theorem sdet_mul (q : ℝ) (α β γ δ : ℂ) :
    sdet q (α * γ - q * (β * conj δ)) (α * δ + β * conj γ) = sdet q α β * sdet q γ δ := by
  apply ofReal_injective
  rw [ofReal_mul, ← det_shape, ← det_shape, ← det_shape, ← det_mul, shape_mul]

/-- For two boosts: `|1 - q * conj u * v| ^ 2 + q * |u + v| ^ 2 = (1 + q * |u| ^ 2) *
(1 + q * |v| ^ 2)`. -/
theorem boost_det (q : ℝ) (u v : ℂ) :
    normSq (den q u v) + q * normSq (u + v) = (1 + q * normSq u) * (1 + q * normSq v) := by
  simp only [den, tX, tY, dot, cross, normSq_apply, add_re, add_im]; ring

/-- **The light wall in two dimensions** (`q = -1`):
`|1 + conj u * v| ^ 2 - |u + v| ^ 2 = (1 - |u| ^ 2) * (1 - |v| ^ 2)`. -/
theorem light_wall_2d (u v : ℂ) :
    normSq (den (-1) u v) - normSq (u + v) = (1 - normSq u) * (1 - normSq v) := by
  linear_combination boost_det (-1) u v

/-- **Below light stays below light, in any directions**: if `|u|, |v| < 1` then
`|u + v| < |1 + conj u * v|`, so `|oplus (-1) u v| < 1`. -/
theorem under_light_2d (u v : ℂ) (hu : normSq u < 1) (hv : normSq v < 1) :
    normSq (u + v) < normSq (den (-1) u v) := by
  have hp := mul_pos (sub_pos.2 hu) (sub_pos.2 hv)
  linarith [light_wall_2d u v]

/-- At the sphere end `q = 1` there is no wall: the determinant of a boost is at least `1`. -/
theorem sphere_no_wall (u : ℂ) : 1 ≤ sdet 1 1 u := by
  simpa [sdet] using normSq_nonneg u

/-! ### Half a kick, and the tie to mult-add -/

/-- The half-velocity `t` composed with itself gives the velocity `2 * t / (1 + t ^ 2)`. -/
theorem half_kick {K : Type*} [Field K] (t : K) : vel t t = 2 * t / (1 + t * t) := by
  rw [vel, two_mul]

/-- On the mult-add ruler, the velocity is the half-velocity mult-added to itself: the
half-velocity is half a kick. -/
theorem half_kick_ruler {K : Type*} [Field K] (t : K) (h : 1 - t ≠ 0) (h2 : 1 + t * t ≠ 0) :
    dop (vel t t) = dop t ✳ dop t :=
  dop_vel t t h h h2

/-- Both parts of the twist's tangent `Y / X` come from the one corner `conj u * v`. -/
theorem twist_from_corner (q : ℝ) (u v : ℂ) :
    tX q u v = 1 - q * (conj u * v).re ∧ tY q u v = q * (conj u * v).im := by
  rw [corner_split]; exact ⟨rfl, rfl⟩

end TwistedBoosts
