module

public import FusedAlgebra.Basic

/-!
# Conservation of mass and charge

**Mass.** The invariant mass of two parts is a law of cosines,
`M ^ 2 = m₁ ^ 2 + m₂ ^ 2 + 2 * m₁ * m₂ * c` (`mass2`). `c = 1` is Lavoisier (`M = m₁ + m₂`
exactly: the border, a flat triangle, `lavoisier`); `c < 1` is a bound state (`bound_lighter`);
`c > 1` is a free pair in motion (`free_heavier`).

**Charge.** Charges add, and the fugacity `z ^ Q` multiplies (`pow_add`), so its excess
`z ^ Q - 1` mult-adds (`fugacity_multAdd`). Every additive conservation law is the log of a
mult-add lift.

`circle_closed` is the Brahmagupta–Fibonacci identity on the unit circle: in mathlib,
multiplicativity of `Complex.normSq` (`map_mul`), or closure of `Circle` under multiplication.

The core-Lean original worked over `ℚ`; here the identities hold in any commutative ring and the
inequalities in any linearly ordered commutative ring.
-/

@[expose] public section

namespace ChemConservation

open MultAdd

variable {R : Type*}

section CommRing

variable [CommRing R]

/-- Invariant mass squared of two parts with "cosine" `c`. -/
def mass2 (m1 m2 c : R) : R := m1 ^ 2 + m2 ^ 2 + 2 * m1 * m2 * c

/-- Lavoisier's border: `c = 1` gives exactly `(m₁ + m₂) ^ 2`. -/
theorem lavoisier (m1 m2 : R) : mass2 m1 m2 1 = (m1 + m2) ^ 2 := by
  unfold mass2; ring

/-- The mass deficit squared is `2 * m₁ * m₂ * (1 - c)`. -/
theorem deficit (m1 m2 c : R) : (m1 + m2) ^ 2 - mass2 m1 m2 c = 2 * m1 * m2 * (1 - c) := by
  unfold mass2; ring

/-- With a binding energy `B`, `M = m₁ + m₂ - B`:
`B * (2 * (m₁ + m₂) - B) = 2 * m₁ * m₂ * (1 - c)`, exactly. -/
theorem binding_exact (m1 m2 c B : R) (h : (m1 + m2 - B) ^ 2 = mass2 m1 m2 c) :
    B * (2 * (m1 + m2) - B) = 2 * m1 * m2 * (1 - c) := by
  unfold mass2 at h
  linear_combination -h

/-- The cosine law reads as a weighted mult-add in the masses' squares: with `a = m₁ ^ 2`,
`b = m₂ ^ 2` and the corner taken through `√(a * b) = m₁ * m₂`, the extra piece beyond plain
addition is `2 * c * (m₁ * m₂)`. At `c = 1` it is the perfect square. -/
theorem cosine_extra (m1 m2 c : R) : mass2 m1 m2 c - (m1 ^ 2 + m2 ^ 2) = 2 * c * (m1 * m2) := by
  unfold mass2; ring

/-- Small-angle form: if `c = 1 - t / 2` (`t = θ ^ 2`), the squared deficit is `m₁ * m₂ * t`.
Stated without division, with `s = t / 2`. -/
theorem small_angle (m1 m2 s : R) : (m1 + m2) ^ 2 - mass2 m1 m2 (1 - s) = m1 * m2 * (2 * s) := by
  unfold mass2; ring

/-! ### Charge: additive charges, multiplicative fugacities -/

/-- The excess `z ^ Q - 1` of the fugacity mult-adds as charges add. -/
theorem fugacity_multAdd (z : R) (Q1 Q2 : ℕ) :
    z ^ (Q1 + Q2) - 1 = (z ^ Q1 - 1) ✳ (z ^ Q2 - 1) := by
  rw [pow_add]; unfold multAdd; ring

/-- A neutral system (charge `0`) has fugacity `1`: excess `0`, mult-add's identity. -/
theorem neutral_identity (z : R) : z ^ 0 - 1 = 0 := by simp

/-- Rotation (the phase `e ^ (i * Q * φ)`) in coordinates: if both factors lie on the unit
circle, so does their product. -/
theorem circle_closed (a b c d : R) (h1 : a ^ 2 + b ^ 2 = 1) (h2 : c ^ 2 + d ^ 2 = 1) :
    (a * c - b * d) ^ 2 + (a * d + b * c) ^ 2 = 1 := by
  linear_combination (c ^ 2 + d ^ 2) * h1 + h2

end CommRing

section Ordered

variable [CommRing R] [LinearOrder R] [IsStrictOrderedRing R]

/-- Bound (`c < 1`) means the whole weighs less than its parts, for positive masses. -/
theorem bound_lighter (m1 m2 c : R) (h1 : 0 < m1) (h2 : 0 < m2) (hc : c < 1) :
    mass2 m1 m2 c < (m1 + m2) ^ 2 := by
  rw [← sub_pos, deficit]
  exact mul_pos (mul_pos (mul_pos two_pos h1) h2) (sub_pos.2 hc)

/-- Free pairs in motion (`c > 1`) weigh more than their rest masses. -/
theorem free_heavier (m1 m2 c : R) (h1 : 0 < m1) (h2 : 0 < m2) (hc : 1 < c) :
    (m1 + m2) ^ 2 < mass2 m1 m2 c := by
  rw [← sub_neg, deficit]
  exact mul_neg_of_pos_of_neg (mul_pos (mul_pos two_pos h1) h2) (sub_neg.2 hc)

end Ordered

end ChemConservation
