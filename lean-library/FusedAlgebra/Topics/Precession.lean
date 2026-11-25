/-
  Precession through the circle operation  (DRAFT, Nov 25 2026)
  ==============================================================

  Companion: ../README.md and ../checks.py. Core Lean 4.37.0-rc1 only (no mathlib),
  no `sorry`. Imports Base only. Quaternions are 4-tuples of rationals.

  The circle operation on any associative algebra: a ∘ₖ b = a + b + k·ab.
  It keeps its lift, 1 + k(a ∘ₖ b) = (1 + ka)(1 + kb), even when ab ≠ ba.

    1. Swap = turn. Swapping the inputs gives the same result turned by the first
       input:  (1 + ka)·(1 + k(b ∘ₖ a)) = (1 + k(a ∘ₖ b))·(1 + ka).
    2. The swap's size is k times the commutator, and for two pure quaternions
       (two spatial directions) that is 2k·(u × v): a direction at right angles
       to both inputs. This is the "90 degrees" of every precession below.
    3. Swapping never changes the size (the quaternion norm is multiplicative).
    4. Inputs along one line commute, so there is no turn: one-dimensional
       versions of the theory have no precession.
    5. Rodrigues: (1 + g)(1 + h) = (1 − g·h) + (g + h + g × h) for pure g, h.
       The dot part of the corner is the plane dial d = +1 (the tan law's
       denominator); the cross part is the tilt.
    6. Gyroscope: torque is horizontal and at 90° to the axle; steady precession
       of a horizontal axle needs exactly Ωₚ·I₃ω₃ = m·g·l.
    7. Larmor: the turn B × S keeps the spin's length and its angle to B.
    8. Spinning BTZ: drag speed u = (ρ − 1)/(ρ + 1) with ρ = T_R/T_L; its
       relativity-page coordinate 2u/(1 − u) is ρ − 1, so stacking two drags
       multiplies the temperature ratios.

  What Lean does not check (real numbers, integrals): twist = d × area, Foucault,
  Montgomery's rigid-body phase, the heavy-top simulation, Berry's half angle,
  the BTZ gyroscope rate. Those are in ../checks.py.
-/
import FusedAlgebra.Base

namespace Precession

/-! ## 0. Quaternions over the rationals -/

structure Q where
  w : Rat
  x : Rat
  y : Rat
  z : Rat

@[ext] theorem Q.ext' {a b : Q} (h0 : a.w = b.w) (h1 : a.x = b.x) (h2 : a.y = b.y)
    (h3 : a.z = b.z) : a = b := by
  cases a; cases b; simp_all

def Q.one : Q := ⟨1, 0, 0, 0⟩
def Q.add (a b : Q) : Q := ⟨a.w + b.w, a.x + b.x, a.y + b.y, a.z + b.z⟩
def Q.sub (a b : Q) : Q := ⟨a.w - b.w, a.x - b.x, a.y - b.y, a.z - b.z⟩
def Q.smul (r : Rat) (a : Q) : Q := ⟨r * a.w, r * a.x, r * a.y, r * a.z⟩
def Q.mul (a b : Q) : Q :=
  ⟨a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z,
   a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y,
   a.w * b.y - a.x * b.z + a.y * b.w + a.z * b.x,
   a.w * b.z + a.x * b.y - a.y * b.x + a.z * b.w⟩
def Q.nsq (a : Q) : Rat := a.w * a.w + a.x * a.x + a.y * a.y + a.z * a.z

/-- The circle operation on quaternions: a ∘ₖ b = a + b + k·ab. -/
def circ (k : Rat) (a b : Q) : Q := Q.add (Q.add a b) (Q.smul k (Q.mul a b))

/-- The lift 1 + k·a. -/
def lift (k : Rat) (a : Q) : Q := Q.add Q.one (Q.smul k a)

theorem Q.mul_assoc (a b c : Q) : Q.mul (Q.mul a b) c = Q.mul a (Q.mul b c) := by
  ext <;> simp [Q.mul] <;> grind

/-- The lift survives non-commutativity. -/
theorem lift_circ (k : Rat) (a b : Q) :
    lift k (circ k a b) = Q.mul (lift k a) (lift k b) := by
  ext <;> simp [lift, circ, Q.add, Q.smul, Q.mul, Q.one] <;> grind

theorem circ_assoc (k : Rat) (a b c : Q) :
    circ k (circ k a b) c = circ k a (circ k b c) := by
  ext <;> simp [circ, Q.add, Q.smul, Q.mul] <;> grind

theorem circ_zero (k : Rat) (a : Q) : circ k a ⟨0, 0, 0, 0⟩ = a := by
  ext <;> simp [circ, Q.add, Q.smul, Q.mul] <;> grind

/-! ## 1. Swap = turn -/

/-- **Swap = turn**: the swapped result is the original turned (conjugated) by the
    first input's lift. With L = 1 + ka:  L · lift(b ∘ a) = lift(a ∘ b) · L. -/
theorem swap_is_turn (k : Rat) (a b : Q) :
    Q.mul (lift k a) (lift k (circ k b a)) = Q.mul (lift k (circ k a b)) (lift k a) := by
  rw [lift_circ, lift_circ, ← Q.mul_assoc]

/-- The swap's size is k times the commutator. -/
theorem swap_diff (k : Rat) (a b : Q) :
    Q.sub (circ k a b) (circ k b a) = Q.smul k (Q.sub (Q.mul a b) (Q.mul b a)) := by
  ext <;> simp [circ, Q.add, Q.sub, Q.smul, Q.mul] <;> grind

/-! ## 2. The 90 degrees: the commutator of two directions is twice their cross product -/

def dot3 (u v : Q) : Rat := u.x * v.x + u.y * v.y + u.z * v.z
def cross3 (u v : Q) : Q :=
  ⟨0, u.y * v.z - u.z * v.y, u.z * v.x - u.x * v.z, u.x * v.y - u.y * v.x⟩
/-- A pure quaternion: a spatial direction (real part 0). -/
def pure (a b c : Rat) : Q := ⟨0, a, b, c⟩

theorem pure_commutator (a b c d e f : Rat) :
    Q.sub (Q.mul (pure a b c) (pure d e f)) (Q.mul (pure d e f) (pure a b c))
      = Q.smul 2 (cross3 (pure a b c) (pure d e f)) := by
  ext <;> simp [Q.sub, Q.mul, Q.smul, cross3, pure] <;> grind

theorem swap_pure (k a b c d e f : Rat) :
    Q.sub (circ k (pure a b c) (pure d e f)) (circ k (pure d e f) (pure a b c))
      = Q.smul (2 * k) (cross3 (pure a b c) (pure d e f)) := by
  ext <;> simp [circ, Q.add, Q.sub, Q.smul, Q.mul, cross3, pure] <;> grind

theorem cross_perp_left (u v : Q) : dot3 u (cross3 u v) = 0 := by
  simp [dot3, cross3]; grind

theorem cross_perp_right (u v : Q) : dot3 v (cross3 u v) = 0 := by
  simp [dot3, cross3]; grind

/-! ## 3. Swapping keeps the size -/

/-- Euler's four-square identity: the quaternion norm multiplies. -/
theorem nsq_mul (a b : Q) : Q.nsq (Q.mul a b) = Q.nsq a * Q.nsq b := by
  simp [Q.nsq, Q.mul]; grind

theorem swap_keeps_size (k : Rat) (a b : Q) :
    Q.nsq (lift k (circ k a b)) = Q.nsq (lift k (circ k b a)) := by
  rw [lift_circ, lift_circ, nsq_mul, nsq_mul]; grind

/-! ## 4. One line, no turn -/

theorem collinear_commute (k r s : Rat) (c : Q) :
    circ k (Q.smul r c) (Q.smul s c) = circ k (Q.smul s c) (Q.smul r c) := by
  ext <;> simp [circ, Q.add, Q.smul, Q.mul] <;> grind

/-! ## 5. Rodrigues: the corner splits into the plane dial (dot) and the tilt (cross) -/

theorem rodrigues (a b c d e f : Rat) :
    Q.mul (Q.add Q.one (pure a b c)) (Q.add Q.one (pure d e f))
      = Q.add ⟨1 - dot3 (pure a b c) (pure d e f), 0, 0, 0⟩
          (Q.add (Q.add (pure a b c) (pure d e f)) (cross3 (pure a b c) (pure d e f))) := by
  ext <;> simp [Q.mul, Q.add, Q.one, pure, dot3, cross3] <;> grind

/-- Same axis: the cross part vanishes and the tan law (a + b)/(1 − ab) is left. -/
theorem rodrigues_same_axis (a b : Rat) :
    Q.mul (Q.add Q.one (pure a 0 0)) (Q.add Q.one (pure b 0 0)) = ⟨1 - a * b, a + b, 0, 0⟩ := by
  ext <;> simp [Q.mul, Q.add, Q.one, pure] <;> grind

/-! ## 6. The gyroscope -/

/-- Torque of gravity on an axle a with the centre of mass at distance l. -/
def torque (m g l : Rat) (a : Q) : Q := cross3 (Q.smul l a) (pure 0 0 (-(m * g)))

theorem torque_perp_axle (m g l : Rat) (a : Q) : dot3 a (torque m g l a) = 0 := by
  simp [torque, dot3, cross3, Q.smul, pure]; grind

theorem torque_horizontal (m g l : Rat) (a : Q) : (torque m g l a).z = 0 := by
  simp [torque, cross3, Q.smul, pure]; grind

/-- Steady precession of a horizontal axle along x: the turn of the angular momentum
    by the precession, Ωₚ ẑ × L with L = I₃ω₃ x̂ + I₁Ωₚ ẑ, equals the torque exactly
    when Ωₚ·I₃ω₃ = m·g·l. -/
theorem gyroscope_steady (m g l I1 I3 w3 Wp : Rat) (h : Wp * (I3 * w3) = m * g * l) :
    cross3 (pure 0 0 Wp) (pure (I3 * w3) 0 (I1 * Wp)) = torque m g l (pure 1 0 0) := by
  ext <;> simp [cross3, torque, pure, Q.smul] <;> grind

/-! ## 7. Larmor: precession keeps length and tilt -/

theorem larmor_keeps_length (B S : Q) : dot3 S (cross3 B S) = 0 := by
  simp [dot3, cross3]; grind

theorem larmor_keeps_tilt (B S : Q) : dot3 B (cross3 B S) = 0 := by
  simp [dot3, cross3]; grind

/-! ## 8. Spinning BTZ: drag speed from the two temperatures -/

/-- Drag speed (horizon angular velocity times ℓ) from ρ = T_R / T_L. -/
def drag (ρ : Rat) : Rat := (ρ - 1) / (ρ + 1)

/-- Its relativity-page coordinate 2u/(1 − u) (Doppler factor squared minus one)
    is ρ − 1, so the ratio of temperatures is the lift. -/
theorem drag_coordinate (ρ : Rat) (hρ : ρ + 1 ≠ 0) :
    2 * drag ρ / (1 - drag ρ) = ρ - 1 := by
  unfold drag
  have e : (1 : Rat) - (ρ - 1) / (ρ + 1) = 2 / (ρ + 1) := by grind
  rw [e]
  grind

/-- Two drags stacked by Einstein addition multiply the temperature ratios. -/
theorem drag_stack (ρ σ : Rat) (hρ : 0 < ρ) (hσ : 0 < σ) :
    (drag ρ + drag σ) / (1 + drag ρ * drag σ) = drag (ρ * σ) := by
  unfold drag
  have h1 : ρ + 1 ≠ 0 := by intro h; grind
  have h2 : σ + 1 ≠ 0 := by intro h; grind
  have h3 : ρ * σ + 1 ≠ 0 := by
    have : 0 < ρ * σ := Rat.mul_pos hρ hσ
    intro h; grind
  have hn : (ρ - 1) / (ρ + 1) + (σ - 1) / (σ + 1) = 2 * (ρ * σ - 1) / ((ρ + 1) * (σ + 1)) := by
    grind
  have hd : 1 + (ρ - 1) / (ρ + 1) * ((σ - 1) / (σ + 1)) = 2 * (ρ * σ + 1) / ((ρ + 1) * (σ + 1)) := by
    grind
  rw [hn, hd]
  have h12 : (ρ + 1) * (σ + 1) ≠ 0 := by
    intro h0; rcases Rat.mul_eq_zero.mp h0 with h | h <;> contradiction
  grind

/-- The parallel thread's spinning bath (left movers 4× hotter than right movers):
    the horizon drags at 3/5 of light speed, against the hotter side. -/
theorem spinning_bath_drag : drag (1 / 4) = -3 / 5 := by
  unfold drag; grind

end Precession
