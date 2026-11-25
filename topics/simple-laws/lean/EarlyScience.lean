/-
  Fused Hyperoperation Algebra: early science, era by era, on the landscape.
  Core Lean 4.37 only (no mathlib), no `sorry`. Imports the library.

  Two kinds of law:
    * combination laws (two things of one kind make one thing of that kind) get a point (k, q);
    * product laws (y = c·x^m·z^n: Ohm, Newton's F = ma, Boyle, Kepler III, gravity) are not
      points. They are the circle operation (called mult-add before 1.0) itself, read on percent changes: every product law makes
      percent changes combine by the circle operation, and its exponents are circle powers.
-/
import FusedAlgebra.Base

namespace SimpleLaws
open CircleOp

/-! ## Counting -/

/-- non-empty choices: a set of m things has 2^m − 1 non-empty subsets; joining two separate
    sets combines those counts by the circle operation (pick from either, or both, but not nothing) -/
theorem nonempty_subsets (m n : Nat) :
    circ ((2 : Int) ^ m - 1) ((2 : Int) ^ n - 1) = (2 : Int) ^ (m + n) - 1 := by
  unfold circ; rw [pow_add']; grind

/-! ## Measurement -/

/-- relative errors in a product combine by the circle operation, exactly (the textbook rule "they add"
    is the k = 0 shadow that drops the corner e₁e₂) -/
theorem product_error (x y e1 e2 : Int) :
    (x * (1 + e1)) * (y * (1 + e2)) = (x * y) * (1 + circ e1 e2) := by
  unfold circ; grind

/-- the corner that the first-order rule drops -/
theorem product_error_corner (e1 e2 : Int) : circ e1 e2 - (e1 + e2) = e1 * e2 := by
  unfold circ; grind

/-! ## Archimedes' lever and Boyle's law: products held constant -/

/-- if one factor of a constant product grows by a, the other must change by b with a ◯ b = 0:
    the lever (weight × arm) and Boyle (pressure × volume) are the same picture -/
theorem constant_product (P V a b : Int) (hPV : P * V ≠ 0) :
    (P * (1 + a)) * (V * (1 + b)) = P * V ↔ circ a b = 0 := by
  have e : (P * (1 + a)) * (V * (1 + b)) = P * V + P * V * circ a b := by unfold circ; grind
  rw [e]
  constructor
  · intro h
    have h0 : P * V * circ a b = 0 := by omega
    rcases Int.mul_eq_zero.mp h0 with h1 | h1
    · exact absurd h1 hPV
    · exact h1
  · intro h; rw [h]; grind

/-- on a field the partner change is the circle inverse −a/(1 + a) -/
theorem boyle_partner {K : Type} [Lean.Grind.Field K] (a b : K) (ha : 1 + a ≠ 0)
    (h : circ a b = 0) : b = circInv a := by
  have h1 : b * (1 + a) = -a := by unfold circ at h; grind
  unfold circInv
  grind

/-- pressure gains combine at k = +1 while the matching volume losses combine at k = −1:
    if (1 + a)(1 − d₁) = 1 and (1 + b)(1 − d₂) = 1 then the totals match the same way -/
theorem boyle_gain_loss (a b d1 d2 : Int) (h1 : (1 + a) * (1 - d1) = 1) (h2 : (1 + b) * (1 - d2) = 1) :
    (1 + circ a b) * (1 - circK (-1) d1 d2) = 1 := by
  rw [circ_lift, circK_either]
  have : (1 + a) * (1 + b) * (1 - (1 - (1 - d1) * (1 - d2))) = ((1 + a) * (1 - d1)) * ((1 + b) * (1 - d2)) := by grind
  rw [this, h1, h2]; grind

/-! ## Kepler's third law and Newton's inverse square: circle powers -/

/-- T² ∝ a³ on percent changes: (1 + t)² = (1 + α)³ is t ◯ t = α ◯ α ◯ α
    (the square's gnomon equals the cube's layers) -/
theorem kepler_third (t α : Int) :
    (1 + t) * (1 + t) = (1 + α) * (1 + α) * (1 + α) ↔ circ t t = circ (circ α α) α := by
  rw [← circ_lift, ← circ_lift, ← circ_lift]; omega

/-- four times as far out, eight times as long a year -/
example : (1 + 7) * (1 + 7) = (1 + 3) * (1 + 3) * (1 + 3) := by decide

/-- gravity: when the distance grows by ρ the force changes by the circle inverse of ρ ◯ ρ -/
theorem inverse_square {K : Type} [Lean.Grind.Field K] (ρ : K) (h : 1 + circ ρ ρ ≠ 0) :
    ((1 + ρ) * (1 + ρ))⁻¹ = 1 + circInv (circ ρ ρ) := by
  rw [← circ_lift]; exact inv_one_add _ h

/-! ## Newton's law of cooling: the gap left behind compounds at k = −1 -/

/-- each step closes a fraction a of the temperature gap; after n steps the part closed is
    1 − (1 − a)ⁿ, the union-of-failures power -/
theorem cooling (a : Int) (n : Nat) : wpow (-1) a n = 1 - (1 - a) ^ n := wpow_neg_one a n

/-! ## Optics: lenses in contact add; lenses apart are the circle operation with k = −d -/

/-- a 2×2 ray matrix [[A, B], [C, D]] -/
structure Ray where
  A : Int
  B : Int
  C : Int
  D : Int

def Ray.mul (M N : Ray) : Ray :=
  ⟨M.A * N.A + M.B * N.C, M.A * N.B + M.B * N.D, M.C * N.A + M.D * N.C, M.C * N.B + M.D * N.D⟩

/-- a thin lens of power P, and a gap of length d -/
def lens (P : Int) : Ray := ⟨1, 0, -P, 1⟩
def gap (d : Int) : Ray := ⟨1, d, 0, 1⟩

/-- two lenses in contact: powers add (the point (0, 0)) -/
theorem lenses_in_contact (P1 P2 : Int) : ((lens P2).mul (lens P1)).C = -(P1 + P2) := by
  simp only [Ray.mul, lens]; grind

/-- two lenses a distance d apart: total power P₁ + P₂ − d·P₁P₂, the circle operation at k = −d -/
theorem lenses_apart (P1 d P2 : Int) :
    ((lens P2).mul ((gap d).mul (lens P1))).C = -(circK (-d) P1 P2) := by
  simp only [Ray.mul, lens, gap, circK]; grind

/-- each lens's squeeze factor 1 − d·P (how much a parallel beam narrows over the gap) multiplies -/
theorem squeeze_multiplies (P1 d P2 : Int) :
    1 - d * circK (-d) P1 P2 = (1 - d * P1) * (1 - d * P2) := by
  unfold circK; grind

/-- the wall: a first lens that focuses exactly on the second (d·P₁ = 1) makes the second lens
    change nothing: the field lens of telescope eyepieces -/
theorem field_lens (P1 d P2 : Int) (h : d * P1 = 1) : circK (-d) P1 P2 = P1 := by
  unfold circK
  have : -d * P1 * P2 = -(d * P1) * P2 := by grind
  rw [this, h]; grind

/-! ## Circuits and springs: series and parallel -/

/-- parallel resistors (or springs in series, capacitors in series) are addition read through
    the inversion lens x ↦ 1/x: (RS/(R + S))·(1/R + 1/S) = 1 -/
theorem parallel_is_inverted_addition {K : Type} [Lean.Grind.Field K] (R S : K)
    (hR : R ≠ 0) (hS : S ≠ 0) (hRS : R + S ≠ 0) :
    (R * S / (R + S)) * (R⁻¹ + S⁻¹) = 1 := by
  grind

/-- the resistor ladder (1 Ω in series, then 1 Ω across, repeated): one section sends a load
    Z = n/d to 1 + Z/(1 + Z) = (2n + d)/(n + d) -/
def ladderStep (p : Int × Int) : Int × Int := (2 * p.1 + p.2, p.1 + p.2)

def ladder : Nat → Int × Int
  | 0 => (0, 1)
  | n + 1 => ladderStep (ladder n)

def fib : Nat → Int
  | 0 => 0
  | 1 => 1
  | n + 2 => fib (n + 1) + fib n

/-- n sections ending in a short circuit give F(2n)/F(2n − 1): 1, 3/2, 8/5, 21/13, ... -/
theorem ladder_fib (n : Nat) : ladder (n + 1) = (fib (2 * n + 2), fib (2 * n + 1)) := by
  induction n with
  | zero => decide
  | succ n ih =>
    have e1 : 2 * (n + 1) + 2 = (2 * n + 2) + 2 := by omega
    have e2 : 2 * (n + 1) + 1 = (2 * n + 1) + 2 := by omega
    rw [e1, e2]
    simp only [ladder]
    rw [show ladderStep (ladder n) = ladder (n + 1) from rfl, ih]
    simp only [ladderStep, fib]
    simp only [Prod.mk.injEq]
    constructor <;> grind

/-- the walls of one section are its fixed points: Z = (1 + 2Z)/(1 + Z) ⇔ Z² = Z + 1,
    the golden ratio and −1/golden ratio (two real walls) -/
theorem ladder_walls (n d : Int) : (2 * n + d) * d = n * (n + d) ↔ n * n = n * d + d * d := by
  constructor <;> intro h <;> grind

/-- every ladder value sits on the hyperbola n² − nd − d² = −1, which hugs the golden wall:
    so Z² − Z − 1 = −1/d² and the ladder walks into the wall -/
theorem ladder_hyperbola (k : Nat) : (ladder k).1 * (ladder k).1 - (ladder k).1 * (ladder k).2
    - (ladder k).2 * (ladder k).2 = -1 := by
  induction k with
  | zero => decide
  | succ k ih => simp only [ladder, ladderStep]; grind

end SimpleLaws
