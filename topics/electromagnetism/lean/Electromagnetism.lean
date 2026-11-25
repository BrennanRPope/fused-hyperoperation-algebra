/-
  Electromagnetism through the circle operation (called mult-add before 1.0): the elementary parts
  =======================================================

  Companion checks: ../checks.py (sympy, numpy, scipy). Core Lean 4.37.0-rc1, no mathlib,
  no `sorry`. Imports Base, Topics.Thermodynamics (velocity addition `vel`) and
  Topics.Relativity (`dop`, the Doppler-squared-minus-one variable) and Topics.TwistedBoosts
  (complex numbers as pairs of rationals).

  What is proved here (the algebra; the calculus, Clifford products and numerics are in the checks):
  1. The dial is set by the wave equation: a dial boost x' = g(x − vt), t' = g(t + qvx) leaves
     the wave operator's cross term at zero exactly when c²q = −1, and then the operator keeps
     its wave speed.
  2. Boosting E and B along x, written with the Doppler factor D (all rational): the light-cone
     parts E_y ± B_z, E_z ∓ B_y are scaled by D⁻¹ and D; boosts compose by multiplying D; the
     speed of D₁D₂ is Einstein's sum; D² − 1 composes by the circle operation; E² − B² and E·B are invariant;
     null fields stay null; a light wave's amplitude scales by 1/D; a pure E field gains B.
  3. The gauge phase: a complex a is a unit phase exactly when a ◯ ā = 0; unit phases are closed
     under ◯. The split-complex (boost) version is the same rule with j² = +1.
  4. Lattice gauge: a plaquette is the ◯ of its four links (two run backwards, as circle
     inverses); two plaquettes ◯ into the big loop (discrete Stokes).
  5. Charge conservation in momentum space: k_μ k_ν F^{μν} = 0 for any antisymmetric F.
  6. Born–Infeld's field cap: E² = D²/(b² + D²) · b² stays below b², with the exact inverse.
-/
import FusedAlgebra.Base
import FusedAlgebra.Topics.Thermodynamics
import FusedAlgebra.Topics.Relativity
import FusedAlgebra.Topics.TwistedBoosts

namespace Electromagnetism
open CircleOp ThermoFused Relativity

open scoped CircleOp

/-! ## 1. The wave equation picks the dial -/

/-- Cross-term coefficient of ∂_X∂_T in ∂_t² − c²∂_x² after the dial boost (divided by −2g²). -/
theorem wave_cross_zero_iff (v c q : Rat) (hv : v ≠ 0) :
    v * (1 + c * c * q) = 0 ↔ c * c * q = -1 := by
  constructor
  · intro h
    have h' : 1 + c * c * q = 0 := by
      rcases Rat.mul_eq_zero.mp h with h1 | h1
      · exact absurd h1 hv
      · exact h1
    grind
  · intro h; grind

/-- Galilean boosts (q = 0) leave a cross term −2g²v: the wave equation is not Galilean. -/
theorem galilean_cross (v c : Rat) (hv : v ≠ 0) : v * (1 + c * c * 0) ≠ 0 := by
  grind

/-- At c²q = −1 the T-T coefficient is −1/c² times the X-X coefficient: same wave speed. -/
theorem wave_speed_kept (v c q : Rat) (hc : c ≠ 0) (h : c * c * q = -1) :
    (1 - c * c * q * q * v * v) * (c * c) = -(v * v - c * c) := by
  have hq : q * (c * c) = -1 := by grind
  have : c * c * q * q * v * v * (c * c) = (q * (c * c)) * (q * (c * c)) * v * v := by grind
  grind

/-! ## 2. The fields ride the boosts -/

/-- The six field components (c = 1). -/
structure EB where
  ex : Rat
  ey : Rat
  ez : Rat
  bx : Rat
  by_ : Rat
  bz : Rat

theorem EB.ext' {f g : EB} (h1 : f.ex = g.ex) (h2 : f.ey = g.ey) (h3 : f.ez = g.ez)
    (h4 : f.bx = g.bx) (h5 : f.by_ = g.by_) (h6 : f.bz = g.bz) : f = g := by
  cases f; cases g; simp_all

/-- γ and γv as rational functions of the Doppler factor D. -/
def gam (D : Rat) : Rat := (D * D + 1) / (2 * D)
def gv (D : Rat) : Rat := (D * D - 1) / (2 * D)
/-- The speed of the boost with Doppler factor D. -/
def vD (D : Rat) : Rat := (D * D - 1) / (D * D + 1)

/-- The boost of E and B along x. -/
def boost (D : Rat) (f : EB) : EB :=
  ⟨f.ex, gam D * f.ey - gv D * f.bz, gam D * f.ez + gv D * f.by_,
   f.bx, gam D * f.by_ + gv D * f.ez, gam D * f.bz - gv D * f.ey⟩

theorem gam_sub_gv (D : Rat) (hD : D ≠ 0) : gam D - gv D = 1 / D := by
  unfold gam gv; grind

theorem gam_add_gv (D : Rat) (hD : D ≠ 0) : gam D + gv D = D := by
  unfold gam gv; grind

/-- γ² − (γv)² = 1. -/
theorem gam_sq (D : Rat) (hD : D ≠ 0) : gam D * gam D - gv D * gv D = 1 := by
  have h1 := gam_sub_gv D hD
  have h2 := gam_add_gv D hD
  have : gam D * gam D - gv D * gv D = (gam D - gv D) * (gam D + gv D) := by grind
  rw [this, h1, h2]; grind

/-- **Light-cone part 1**: E_y + B_z is scaled by 1/D. -/
theorem lc_plus (D : Rat) (f : EB) (hD : D ≠ 0) :
    (boost D f).ey + (boost D f).bz = (f.ey + f.bz) / D := by
  have h := gam_sub_gv D hD
  have e : (boost D f).ey + (boost D f).bz = (gam D - gv D) * (f.ey + f.bz) := by
    unfold boost; grind
  rw [e, h]; grind

/-- **Light-cone part 2**: E_y − B_z is scaled by D. -/
theorem lc_minus (D : Rat) (f : EB) (hD : D ≠ 0) :
    (boost D f).ey - (boost D f).bz = D * (f.ey - f.bz) := by
  have h := gam_add_gv D hD
  have e : (boost D f).ey - (boost D f).bz = (gam D + gv D) * (f.ey - f.bz) := by
    unfold boost; grind
  rw [e, h]

/-- The other pair: E_z − B_y scaled by 1/D, E_z + B_y by D. -/
theorem lc_pair2 (D : Rat) (f : EB) (hD : D ≠ 0) :
    (boost D f).ez - (boost D f).by_ = (f.ez - f.by_) / D ∧
    (boost D f).ez + (boost D f).by_ = D * (f.ez + f.by_) := by
  have h1 := gam_sub_gv D hD
  have h2 := gam_add_gv D hD
  have e1 : (boost D f).ez - (boost D f).by_ = (gam D - gv D) * (f.ez - f.by_) := by
    unfold boost; grind
  have e2 : (boost D f).ez + (boost D f).by_ = (gam D + gv D) * (f.ez + f.by_) := by
    unfold boost; grind
  rw [e1, e2, h1, h2]; constructor <;> grind

/-- The ratio of the two light-cone parts is multiplied by D², the circle-operation lift. -/
theorem lc_ratio (D : Rat) (f : EB) (hD : D ≠ 0) :
    ((boost D f).ey - (boost D f).bz) * (f.ey + f.bz)
      = (D * D) * ((f.ey - f.bz) * ((boost D f).ey + (boost D f).bz)) := by
  rw [lc_minus D f hD, lc_plus D f hD]
  have : D * D * ((f.ey - f.bz) * ((f.ey + f.bz) / D)) = D * (f.ey - f.bz) * (f.ey + f.bz) * (D / D) := by
    grind
  rw [this]; grind

/-- D² is the lift of the speed: D² = (1 + v)/(1 − v), and dop v = D² − 1. -/
theorem dop_vD (D : Rat) (h : D * D + 1 ≠ 0) : dop (vD D) = D * D - 1 := by
  unfold dop vD
  have e : 1 - (D * D - 1) / (D * D + 1) = 2 / (D * D + 1) := by grind
  rw [e]; grind

/-- **D² − 1 composes by the circle operation**: (D₁D₂)² − 1 = (D₁² − 1) ◯ (D₂² − 1). -/
theorem doppler_circ (D1 D2 : Rat) :
    (D1 * D2) * (D1 * D2) - 1 = (D1 * D1 - 1) ◯ (D2 * D2 - 1) := by
  unfold circ; grind

/-- So the field's light-cone ratio, under two boosts, rides the circle operation of the two d's. -/
theorem dop_compose (D1 D2 : Rat) (h1 : D1 * D1 + 1 ≠ 0) (h2 : D2 * D2 + 1 ≠ 0)
    (h12 : (D1 * D2) * (D1 * D2) + 1 ≠ 0) :
    dop (vD (D1 * D2)) = dop (vD D1) ◯ dop (vD D2) := by
  rw [dop_vD _ h12, dop_vD _ h1, dop_vD _ h2]; exact doppler_circ D1 D2

/-- **Speeds of composed boosts add by Einstein's rule**: v(D₁D₂) = v(D₁) ⊕ v(D₂). -/
theorem vD_mul (D1 D2 : Rat) (h1 : D1 * D1 + 1 ≠ 0) (h2 : D2 * D2 + 1 ≠ 0)
    (h12 : (D1 * D1 + 1) * (D2 * D2 + 1) + (D1 * D1 - 1) * (D2 * D2 - 1) ≠ 0) :
    vD (D1 * D2) = vel (vD D1) (vD D2) := by
  unfold vD vel
  have den : 1 + (D1 * D1 - 1) / (D1 * D1 + 1) * ((D2 * D2 - 1) / (D2 * D2 + 1))
      = ((D1 * D1 + 1) * (D2 * D2 + 1) + (D1 * D1 - 1) * (D2 * D2 - 1)) / ((D1 * D1 + 1) * (D2 * D2 + 1)) := by
    grind
  have num : (D1 * D1 - 1) / (D1 * D1 + 1) + (D2 * D2 - 1) / (D2 * D2 + 1)
      = ((D1 * D1 - 1) * (D2 * D2 + 1) + (D2 * D2 - 1) * (D1 * D1 + 1)) / ((D1 * D1 + 1) * (D2 * D2 + 1)) := by
    grind
  rw [den, num]
  have hp : (D1 * D1 + 1) * (D2 * D2 + 1) ≠ 0 := by grind
  have k1 : (D1 * D1 + 1) * (D2 * D2 + 1) + (D1 * D1 - 1) * (D2 * D2 - 1) = 2 * ((D1 * D2) * (D1 * D2) + 1) := by grind
  have k2 : (D1 * D1 - 1) * (D2 * D2 + 1) + (D2 * D2 - 1) * (D1 * D1 + 1) = 2 * ((D1 * D2) * (D1 * D2) - 1) := by grind
  rw [k1, k2]
  have h3 : (D1 * D2) * (D1 * D2) + 1 ≠ 0 := by grind
  generalize (D1 * D2) * (D1 * D2) = P at *
  generalize (D1 * D1 + 1) * (D2 * D2 + 1) = Q at *
  grind

/-- The invariant E² − B². -/
def inv1 (f : EB) : Rat := f.ex * f.ex + f.ey * f.ey + f.ez * f.ez - f.bx * f.bx - f.by_ * f.by_ - f.bz * f.bz
/-- The invariant E·B. -/
def inv2 (f : EB) : Rat := f.ex * f.bx + f.ey * f.by_ + f.ez * f.bz

/-- **E² − B² is the same in every frame.** -/
theorem inv1_boost (D : Rat) (f : EB) (hD : D ≠ 0) : inv1 (boost D f) = inv1 f := by
  have h := gam_sq D hD
  unfold inv1 boost; simp only
  have e : (gam D * f.ey - gv D * f.bz) * (gam D * f.ey - gv D * f.bz)
      + (gam D * f.ez + gv D * f.by_) * (gam D * f.ez + gv D * f.by_)
      - (gam D * f.by_ + gv D * f.ez) * (gam D * f.by_ + gv D * f.ez)
      - (gam D * f.bz - gv D * f.ey) * (gam D * f.bz - gv D * f.ey)
      = (gam D * gam D - gv D * gv D) * (f.ey * f.ey + f.ez * f.ez - f.by_ * f.by_ - f.bz * f.bz) := by grind
  grind

/-- **E·B is the same in every frame.** -/
theorem inv2_boost (D : Rat) (f : EB) (hD : D ≠ 0) : inv2 (boost D f) = inv2 f := by
  have h := gam_sq D hD
  unfold inv2 boost; simp only
  have e : (gam D * f.ey - gv D * f.bz) * (gam D * f.by_ + gv D * f.ez)
      + (gam D * f.ez + gv D * f.by_) * (gam D * f.bz - gv D * f.ey)
      = (gam D * gam D - gv D * gv D) * (f.ey * f.by_ + f.ez * f.bz) := by grind
  grind

/-- Light stays light: a null field (E² = B², E ⊥ B) is null in every frame. -/
theorem null_stays_null (D : Rat) (f : EB) (hD : D ≠ 0) (h1 : inv1 f = 0) (h2 : inv2 f = 0) :
    inv1 (boost D f) = 0 ∧ inv2 (boost D f) = 0 := by
  rw [inv1_boost D f hD, inv2_boost D f hD]; exact ⟨h1, h2⟩

/-- A light wave moving along +x: E_y = B_z = a. -/
def wave (a : Rat) : EB := ⟨0, a, 0, 0, 0, a⟩

theorem wave_null (a : Rat) : inv1 (wave a) = 0 ∧ inv2 (wave a) = 0 := by
  unfold inv1 inv2 wave; constructor <;> grind

/-- **The wave's amplitude is Doppler-shifted**: after the boost it is the wave of amplitude a/D. -/
theorem wave_boost (D a : Rat) (hD : D ≠ 0) : boost D (wave a) = wave (a / D) := by
  have h := gam_sub_gv D hD
  have key : gam D * a - gv D * a = a / D := by
    rw [show gam D * a - gv D * a = (gam D - gv D) * a by grind, h]; grind
  apply EB.ext' <;> unfold boost wave <;> simp only <;> grind

/-- Its energy density (E² + B²)/2 = a² is divided by D², the inverse circle-operation lift. -/
theorem wave_energy (D a : Rat) (hD : D ≠ 0) :
    (a / D) * (a / D) * (D * D) = a * a := by
  grind

/-- **Magnetism from electricity**: a pure E_y, boosted, gets B_z = −γv E_y. -/
theorem purcell (D e : Rat) : (boost D ⟨0, e, 0, 0, 0, 0⟩).bz = -(gv D * e) := by
  unfold boost; simp only; grind

/-! ## 3. The gauge phase as a dial -/

open TwistedBoosts in
/-- The circle operation on complex numbers (pairs of rationals). -/
def maC (a b : C) : C := C.add (C.add a b) (C.mul a b)

open TwistedBoosts in
/-- **a ◯ ā = |1 + a|² − 1** (real): so a is a unit phase exactly when a ◯ ā = 0. -/
theorem maC_conj (a : C) : maC a (C.conj a) = ⟨C.nsq (C.add C.one a) - 1, 0⟩ := by
  unfold maC C.add C.mul C.conj C.nsq C.one; simp only; congr 1 <;> grind

open TwistedBoosts in
/-- The lift of a ◯ b is the product of lifts (complex multiplication). -/
theorem maC_lift (a b : C) : C.add C.one (maC a b) = C.mul (C.add C.one a) (C.add C.one b) := by
  unfold maC C.add C.mul C.one; simp only; congr 1 <;> grind

open TwistedBoosts in
/-- |xy|² = |x|²|y|² for complex pairs. -/
theorem nsq_mul (u w : C) : C.nsq (C.mul u w) = C.nsq u * C.nsq w := by
  unfold C.nsq C.mul; simp only; grind

open TwistedBoosts in
/-- **Unit phases are closed under the circle operation**: if a ◯ ā = 0 and b ◯ b̄ = 0, the same holds for a ◯ b. -/
theorem phase_closed (a b : C)
    (ha : C.nsq (C.add C.one a) = 1) (hb : C.nsq (C.add C.one b) = 1) :
    C.nsq (C.add C.one (maC a b)) = 1 := by
  rw [maC_lift, nsq_mul, ha, hb]; grind

open TwistedBoosts in
/-- The wall a = −1 is never a phase: |1 + (−1)|² = 0. -/
theorem wall_not_phase : C.nsq (C.add C.one ⟨-1, 0⟩) = 0 := by
  unfold C.nsq C.add C.one; simp only; grind

/-- Split-complex numbers x + jy with j² = +1 (the boost plane). -/
structure S where
  re : Rat
  im : Rat

def S.mul (x y : S) : S := ⟨x.re * y.re + x.im * y.im, x.re * y.im + x.im * y.re⟩
def S.add (x y : S) : S := ⟨x.re + y.re, x.im + y.im⟩
def S.conj (x : S) : S := ⟨x.re, -x.im⟩
def S.one : S := ⟨1, 0⟩
/-- The Lorentz "norm" x² − y². -/
def S.nsq (x : S) : Rat := x.re * x.re - x.im * x.im
def maS (a b : S) : S := S.add (S.add a b) (S.mul a b)

/-- **The boost version of the same rule**: a ◯ ā = (|1 + a|²_Lorentz − 1). -/
theorem maS_conj (a : S) : maS a (S.conj a) = ⟨S.nsq (S.add S.one a) - 1, 0⟩ := by
  unfold maS S.add S.mul S.conj S.nsq S.one; simp only; congr 1 <;> grind

theorem maS_lift (a b : S) : S.add S.one (maS a b) = S.mul (S.add S.one a) (S.add S.one b) := by
  unfold maS S.add S.mul S.one; simp only; congr 1 <;> grind

theorem S_nsq_mul (u w : S) : S.nsq (S.mul u w) = S.nsq u * S.nsq w := by
  unfold S.nsq S.mul; simp only; grind

/-- Boosts (a ◯ ā = 0 in the split plane) are closed under the circle operation too. -/
theorem boost_closed (a b : S)
    (ha : S.nsq (S.add S.one a) = 1) (hb : S.nsq (S.add S.one b) = 1) :
    S.nsq (S.add S.one (maS a b)) = 1 := by
  rw [maS_lift, S_nsq_mul, ha, hb]; grind

/-- Half-angle charts: (1 + i t₁)(1 + i t₂) = (1 − t₁t₂) + i(t₁ + t₂) and
    (1 + j s₁)(1 + j s₂) = (1 + s₁s₂) + j(s₁ + s₂): the landscape law at q = +1 and q = −1. -/
theorem half_charts (t1 t2 : Rat) :
    TwistedBoosts.C.mul ⟨1, t1⟩ ⟨1, t2⟩ = ⟨1 - t1 * t2, t1 + t2⟩ ∧
    S.mul ⟨1, t1⟩ ⟨1, t2⟩ = ⟨1 + t1 * t2, t1 + t2⟩ := by
  unfold TwistedBoosts.C.mul S.mul; constructor <;> (simp only; congr 1 <;> grind)

/-! ## 4. Lattice gauge: plaquettes compose by the circle operation -/

/-- **A plaquette is the circle operation of its links**, the two backward links entering as
    circle inverses: 1 + (a₁ ◯ a₂ ◯ a₃⁻ ◯ a₄⁻) = (1 + a₁)(1 + a₂)/((1 + a₃)(1 + a₄)). -/
theorem plaquette (a1 a2 a3 a4 : Rat) (h3 : 1 + a3 ≠ 0) (h4 : 1 + a4 ≠ 0) :
    1 + (a1 ◯ a2 ◯ circInv a3 ◯ circInv a4) = (1 + a1) * (1 + a2) / ((1 + a3) * (1 + a4)) := by
  rw [circ_lift, circ_lift, circ_lift, ← inv_one_add a3 h3, ← inv_one_add a4 h4]
  grind

/-- **Discrete Stokes**: with link lifts b₁, s, t₁, l (left square) and b₂, r, t₂, s (right square,
    sharing s), the two plaquettes multiply to the big rectangle and the shared link cancels. -/
theorem stokes (b1 s t1 l b2 r t2 : Rat) (hs : s ≠ 0) (ht1 : t1 ≠ 0) (hl : l ≠ 0) (ht2 : t2 ≠ 0) :
    (b1 * s / (t1 * l)) * (b2 * r / (t2 * s)) = b1 * b2 * r / (t2 * t1 * l) := by
  grind

/-- In circle-operation form: the big loop's value is the ◯ of the two plaquette values. -/
theorem stokes_circ (p1 p2 big : Rat) (h : (1 + p1) * (1 + p2) = 1 + big) : big = p1 ◯ p2 := by
  unfold circ; grind

/-! ## 5. Charge conservation is forced -/

/-- Sum over four indices. -/
def sum4 (f : Fin 4 → Rat) : Rat := f 0 + f 1 + f 2 + f 3

/-- **k_μ k_ν F^{μν} = 0 for any antisymmetric F**: in momentum space, k·J = 0 follows from
    J = k·F alone. -/
theorem charge_conserved (k : Fin 4 → Rat) (F : Fin 4 → Fin 4 → Rat)
    (anti : ∀ i j, F i j = -F j i) :
    sum4 (fun i => sum4 (fun j => k i * k j * F i j)) = 0 := by
  unfold sum4
  have d0 := anti 0 0; have d1 := anti 1 1; have d2 := anti 2 2; have d3 := anti 3 3
  have a01 := anti 0 1; have a02 := anti 0 2; have a03 := anti 0 3
  have a12 := anti 1 2; have a13 := anti 1 3; have a23 := anti 2 3
  grind

/-! ## 6. Born–Infeld's field cap -/

theorem sq_nn (s : Rat) : 0 ≤ s * s := by
  rcases (@Rat.le_total 0 s) with h | h
  · exact Rat.mul_nonneg h h
  · have h1 : 0 ≤ -s := by grind
    have := Rat.mul_nonneg h1 h1
    grind

theorem sq_pos' (b : Rat) (hb : b ≠ 0) : 0 < b * b := by
  have h := sq_nn b
  have hne : b * b ≠ 0 := by
    intro h0; rcases Rat.mul_eq_zero.mp h0 with h1 | h1 <;> exact hb h1
  grind

/-- With E² = b²D²/(b² + D²), the field stays below the wall b² (stated with the positive
    denominator cleared): b²·D² < b²·(b² + D²). -/
theorem bi_cap (b D : Rat) (hb : b ≠ 0) : b * b * (D * D) < b * b * (b * b + D * D) := by
  have h1 := sq_pos' b hb
  have h2 : 0 < (b * b) * (b * b) := Rat.mul_pos h1 h1
  grind

/-- The denominator b² + D² is positive. -/
theorem bi_den_pos (b D : Rat) (hb : b ≠ 0) : 0 < b * b + D * D := by
  have h1 := sq_pos' b hb
  have h2 := sq_nn D
  grind

/-- And the exact inverse: with e = E²/b² = D²/(b² + D²), D²/b² = e/(1 − e). -/
theorem bi_inverse (b D : Rat) (hb : b ≠ 0) :
    (D * D / (b * b + D * D)) / (1 - D * D / (b * b + D * D)) = D * D / (b * b) := by
  have hb2 : b * b ≠ 0 := by have := sq_pos' b hb; grind
  have hden : b * b + D * D ≠ 0 := by have := bi_den_pos b D hb; grind
  have e1 : 1 - D * D / (b * b + D * D) = b * b / (b * b + D * D) := by grind
  rw [e1]
  generalize b * b + D * D = S at *
  generalize b * b = B at *
  grind

end Electromagnetism
