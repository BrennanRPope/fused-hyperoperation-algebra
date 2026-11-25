/-
  Several waves: where the circle operation (called mult-add before 1.0) sits, and the new directions each pair of waves adds
  ===================================================================================
  Companion checks: ../waves/ (rank.py, kp.py, mirror_check.py, relation4c.py).
  Core Lean 4 only (no mathlib), no `sorry`. Imports only FusedAlgebra.Base.
  Intended path: FusedAlgebra/Topics/SeveralWaves.lean.
  Not compiled in the session that wrote it (no Lean there); written in the library's
  `unfold …; grind` style for its first compile. Every identity was checked with sympy.

  Picture. One wave is the circle operation: its tau function is 1 + x, and (1 + a)(1 + b) = 1 + a ◯ b
  (Base's circ_lift). Two and three waves are SolitonsRungs' tau2 and tau3.
  With g waves the tau function is a sum over the subsets S of the waves of
    ∏_{i ∈ S} x_i · ∏_{i < j in S} A_ij,
  so each pair of waves carries one number, the pair weight A_ij. A_ij = 1 for every pair is
  g separate copies of the circle operation; A_ij = 0 is exclusion (the two waves cannot both be on).

  What the file proves:
    A. Four waves with every pair weight 1 is four separate circle operations; a wave whose weights
       to the others are all 1 splits off as its own factor.
    B. KP (the 2+1-dimensional wave equation) fixes each wave by two points p, q on a line:
       the dispersion relation holds for every p, q, and Hirota's two-wave condition forces
       the pair weight to be the cross ratio (p₁ − p₂)(q₁ − q₂) / ((p₁ − q₂)(q₁ − p₂)).
    C. That cross ratio is 1 exactly when one of the waves is trivial (p = q): two genuine
       waves always interact. It is 0 when the waves share a point (resonance). It does not
       change under any affine change p ↦ λp + μ, nor under p ↦ 1/p, so it only depends on
       the 2g points up to the 3-dimensional group of the line.
    D. Counting. g waves have g(g − 1)/2 pair weights; 2g points up to that group give 2g − 3.
       The excess is (g − 2)(g − 3)/2, the same number as the excess of Riemann's period
       matrices (g(g + 1)/2) over curves (3g − 3). Both are zero up to three waves and
       positive from four. The symmetric KdV case (q = −p) has g − 1 numbers and an excess
       (g − 1)(g − 2)/2, positive from three.
    E. The one rule four waves obey: six weights, five numbers, one polynomial `rule4`,
       checked on two configurations of points; an array of all 2s breaks it.
  LIBRARY NOTE: first compiled in the shared library, unchanged, including the four rule4
  decide +kernel checks over Rat. No exact twins.
-/
import FusedAlgebra.Base

namespace SeveralWaves
open CircleOp

open scoped CircleOp

/-! ## A. Waves with pair weights -/

/-- the four-wave tau function: each subset of waves carries the product of its pair weights -/
def tau4 (A12 A13 A14 A23 A24 A34 a b c d : Int) : Int :=
  1 + a + b + c + d
  + A12 * a * b + A13 * a * c + A14 * a * d + A23 * b * c + A24 * b * d + A34 * c * d
  + A12 * A13 * A23 * a * b * c + A12 * A14 * A24 * a * b * d
  + A13 * A14 * A34 * a * c * d + A23 * A24 * A34 * b * c * d
  + A12 * A13 * A14 * A23 * A24 * A34 * a * b * c * d

/-- every pair weight 1: four separate circle operations -/
theorem tau4_free (a b c d : Int) :
    tau4 1 1 1 1 1 1 a b c d = (1 + a) * (1 + b) * (1 + c) * (1 + d) := by
  unfold tau4; grind

/-- a fourth wave with weight 1 to the others splits off as its own circle-operation factor -/
theorem tau4_splits (A12 A13 A23 a b c d : Int) :
    tau4 A12 A13 1 A23 1 1 a b c d
      = (1 + a + b + c + A12 * a * b + A13 * a * c + A23 * b * c
          + A12 * A13 * A23 * a * b * c) * (1 + d) := by
  unfold tau4; grind

/-- weight 0 between waves 1 and 2: no term has both of them on -/
theorem tau4_exclusion (A13 A14 A23 A24 A34 a b c d : Int) :
    tau4 0 A13 A14 A23 A24 A34 a b c d
      = 1 + a + b + c + d + A13 * a * c + A14 * a * d + A23 * b * c + A24 * b * d
        + A34 * c * d + A13 * A14 * A34 * a * c * d + A23 * A24 * A34 * b * c * d := by
  unfold tau4; grind

/-! ## B. KP fixes the pair weight -/

/-- the KP polynomial of a wave number (kx, ky, kt): Hirota's D_x⁴ − 4 D_x D_t + 3 D_y² -/
def kpP (kx ky kt : Int) : Int := kx * kx * kx * kx - 4 * kx * kt + 3 * ky * ky

/-- every pair of points p, q gives a wave: (p − q, p² − q², p³ − q³) solves the dispersion -/
theorem kp_dispersion (p q : Int) :
    kpP (p - q) (p * p - q * q) (p * p * p - q * q * q) = 0 := by
  unfold kpP; grind

/-- Hirota's two-wave condition: P(k₁ + k₂)·A + P(k₁ − k₂) = 0, with A the cross ratio
    (written with denominators cleared) -/
theorem kp_pair_weight (p1 q1 p2 q2 : Int) :
    kpP ((p1 - q1) + (p2 - q2)) ((p1 * p1 - q1 * q1) + (p2 * p2 - q2 * q2))
        ((p1 * p1 * p1 - q1 * q1 * q1) + (p2 * p2 * p2 - q2 * q2 * q2))
        * ((p1 - p2) * (q1 - q2))
      + kpP ((p1 - q1) - (p2 - q2)) ((p1 * p1 - q1 * q1) - (p2 * p2 - q2 * q2))
        ((p1 * p1 * p1 - q1 * q1 * q1) - (p2 * p2 * p2 - q2 * q2 * q2))
        * ((p1 - q2) * (q1 - p2)) = 0 := by
  unfold kpP; grind

/-- the KdV case q = −p: the cross ratio becomes ((p₁ − p₂)/(p₁ + p₂))² -/
theorem kdv_cross_ratio (p1 p2 : Int) :
    (p1 - p2) * (-p1 - -p2) = -((p1 - p2) * (p1 - p2)) ∧
    (p1 - -p2) * (-p1 - p2) = -((p1 + p2) * (p1 + p2)) := by
  constructor <;> grind

/-! ## C. What the cross ratio sees -/

/-- numerator minus denominator of the cross ratio: A − 1 carries (a_i − b_i)(a_j − b_j) -/
theorem cross_ratio_minus_one (ai bi aj bj : Int) :
    (ai - aj) * (bi - bj) - (ai - bj) * (bi - aj) = (ai - bi) * (aj - bj) := by
  grind

/-- the pair weight is 1 exactly when one of the two waves is trivial -/
theorem free_iff_trivial (ai bi aj bj : Int) :
    (ai - aj) * (bi - bj) = (ai - bj) * (bi - aj) ↔ ai = bi ∨ aj = bj := by
  have h := cross_ratio_minus_one ai bi aj bj
  constructor
  · intro e
    have h0 : (ai - bi) * (aj - bj) = 0 := by grind
    rcases Int.mul_eq_zero.mp h0 with h1 | h1
    · left; omega
    · right; omega
  · intro e
    rcases e with e | e
    · subst e; grind
    · subst e; grind

/-- the pair weight is 0 exactly when the waves share a point (resonance) -/
theorem resonance_iff (ai bi aj bj : Int) :
    (ai - aj) * (bi - bj) = 0 ↔ ai = aj ∨ bi = bj := by
  constructor
  · intro h
    rcases Int.mul_eq_zero.mp h with h1 | h1
    · left; omega
    · right; omega
  · intro e
    rcases e with e | e
    · subst e; grind
    · subst e; grind

/-- an affine change p ↦ λp + μ scales numerator and denominator by λ² alike -/
theorem cross_ratio_affine (l m ai bi aj bj : Int) :
    ((l * ai + m) - (l * aj + m)) * ((l * bi + m) - (l * bj + m)) = l * l * ((ai - aj) * (bi - bj)) ∧
    ((l * ai + m) - (l * bj + m)) * ((l * bi + m) - (l * aj + m)) = l * l * ((ai - bj) * (bi - aj)) := by
  constructor <;> grind

/-- p ↦ 1/p keeps the cross ratio (cross-multiplied, with p ↦ 1/p written as p ↦ P/p
    on the common denominator P = a_i b_i a_j b_j) -/
theorem cross_ratio_inversion (ai bi aj bj : Int) :
    ((aj - ai) * (bj - bi)) * ((ai - bj) * (bi - aj))
      = ((ai - aj) * (bi - bj)) * ((bj - ai) * (aj - bi)) := by
  grind

/-! ## D. Counting the new directions -/

/-- pair weights minus the numbers 2g points carry, doubled -/
theorem excess_points (g : Int) : g * (g - 1) - 2 * (2 * g - 3) = (g - 2) * (g - 3) := by
  grind

/-- Riemann's period numbers minus curves, doubled: the same excess -/
theorem excess_curves (g : Int) : g * (g + 1) - 2 * (3 * g - 3) = (g - 2) * (g - 3) := by
  grind

/-- the two counts agree for every number of waves -/
theorem excess_agree (g : Int) :
    g * (g - 1) - 2 * (2 * g - 3) = g * (g + 1) - 2 * (3 * g - 3) := by
  grind

/-- the symmetric KdV waves: g − 1 numbers, excess (g − 1)(g − 2)/2 -/
theorem excess_kdv (g : Int) : g * (g - 1) - 2 * (g - 1) = (g - 1) * (g - 2) := by
  grind

/-- the excess is zero up to three waves and positive from four -/
theorem excess_from_four :
    ∀ g : Nat, g < 40 → 2 ≤ g → (0 < (g - 2) * (g - 3) ↔ 4 ≤ g) := by decide

/-- for KdV it is positive already from three waves -/
theorem kdv_excess_from_three :
    ∀ g : Nat, g < 40 → 2 ≤ g → (0 < (g - 1) * (g - 2) ↔ 3 ≤ g) := by decide

/-! ## E. The rule four waves obey

  Four waves have six pair weights but their eight points, up to the group of the line, carry
  only five numbers. So the six weights obey one polynomial rule. `rule4` is that rule: found
  by exact interpolation, irreducible, symmetric in the four waves, degree 2 in each weight
  (393 terms). sympy shows it vanishes identically on the cross ratios with the first wave at
  (0, ∞) and a₁ = 1 (waves/relcheck2.py); invariance under the line's group (section C) carries
  that to every configuration. Five weights fix the sixth up to two values.
  The library checks it on concrete configurations; the symbolic identity is too large for
  `grind` and stays a sympy check. -/

/-- the pair weight of waves (a_i, b_i) and (a_j, b_j): the cross ratio -/
def crossR (ai bi aj bj : Rat) : Rat := (ai - aj) * (bi - bj) / ((ai - bj) * (bi - aj))

/-- the four-wave rule: zero on every array of weights that comes from points on a line -/
def rule4 (w01 w02 w03 w12 w13 w23 : Rat) : Rat :=
  -w01 * w01 * w02 * w02 * w03 * w03 * w12 * w12 * w13 * w13 * w23 * w23
    + w01 * w01 * w02 * w02 * w03 * w03 * w12 * w12 * w13 * w13 * w23
    + w01 * w01 * w02 * w02 * w03 * w03 * w12 * w12 * w13 * w23 * w23
    - w01 * w01 * w02 * w02 * w03 * w03 * w12 * w12 * w13 * w23
    + w01 * w01 * w02 * w02 * w03 * w03 * w12 * w13 * w13 * w23 * w23
    - w01 * w01 * w02 * w02 * w03 * w03 * w12 * w13 * w13 * w23
    - w01 * w01 * w02 * w02 * w03 * w03 * w12 * w13 * w23 * w23
    + w01 * w01 * w02 * w02 * w03 * w03 * w12 * w13 * w23
    + w01 * w01 * w02 * w02 * w03 * w12 * w12 * w13 * w13 * w23 * w23
    - w01 * w01 * w02 * w02 * w03 * w12 * w12 * w13 * w13 * w23
    - w01 * w01 * w02 * w02 * w03 * w12 * w12 * w13 * w23 * w23
    + 2 * w01 * w01 * w02 * w02 * w03 * w12 * w12 * w13 * w23
    - w01 * w01 * w02 * w02 * w03 * w12 * w12 * w13
    - w01 * w01 * w02 * w02 * w03 * w12 * w12 * w23 + w01 * w01 * w02 * w02 * w03 * w12 * w12
    - w01 * w01 * w02 * w02 * w03 * w12 * w13 * w13 * w23
    + w01 * w01 * w02 * w02 * w03 * w12 * w13 * w13
    - w01 * w01 * w02 * w02 * w03 * w12 * w13 * w23 * w23
    + 2 * w01 * w01 * w02 * w02 * w03 * w12 * w13 * w23 - w01 * w01 * w02 * w02 * w03 * w12 * w13
    + w01 * w01 * w02 * w02 * w03 * w12 * w23 * w23 - w01 * w01 * w02 * w02 * w03 * w12 * w23
    - w01 * w01 * w02 * w02 * w12 * w12 * w13 * w23 + w01 * w01 * w02 * w02 * w12 * w12 * w13
    + w01 * w01 * w02 * w02 * w12 * w12 * w23 - w01 * w01 * w02 * w02 * w12 * w12
    + w01 * w01 * w02 * w02 * w12 * w13 * w23 - w01 * w01 * w02 * w02 * w12 * w13
    - w01 * w01 * w02 * w02 * w12 * w23 + w01 * w01 * w02 * w02 * w12
    + w01 * w01 * w02 * w03 * w03 * w12 * w12 * w13 * w13 * w23 * w23
    - w01 * w01 * w02 * w03 * w03 * w12 * w12 * w13 * w13 * w23
    - w01 * w01 * w02 * w03 * w03 * w12 * w12 * w13 * w23
    + w01 * w01 * w02 * w03 * w03 * w12 * w12 * w13
    - w01 * w01 * w02 * w03 * w03 * w12 * w13 * w13 * w23 * w23
    + 2 * w01 * w01 * w02 * w03 * w03 * w12 * w13 * w13 * w23
    - w01 * w01 * w02 * w03 * w03 * w12 * w13 * w13
    - w01 * w01 * w02 * w03 * w03 * w12 * w13 * w23 * w23
    + 2 * w01 * w01 * w02 * w03 * w03 * w12 * w13 * w23 - w01 * w01 * w02 * w03 * w03 * w12 * w13
    - w01 * w01 * w02 * w03 * w03 * w13 * w13 * w23 + w01 * w01 * w02 * w03 * w03 * w13 * w13
    + w01 * w01 * w02 * w03 * w03 * w13 * w23 * w23 - w01 * w01 * w02 * w03 * w03 * w13 * w23
    - w01 * w01 * w02 * w03 * w12 * w12 * w13 * w13 * w23 * w23
    + w01 * w01 * w02 * w03 * w12 * w12 * w13 * w13 * w23
    - w01 * w01 * w02 * w03 * w12 * w12 * w13 * w23 * w23
    + 2 * w01 * w01 * w02 * w03 * w12 * w12 * w13 * w23 - w01 * w01 * w02 * w03 * w12 * w12 * w13
    + w01 * w01 * w02 * w03 * w12 * w12 * w23 - w01 * w01 * w02 * w03 * w12 * w12
    - w01 * w01 * w02 * w03 * w12 * w13 * w13 * w23 * w23
    + 2 * w01 * w01 * w02 * w03 * w12 * w13 * w13 * w23 - w01 * w01 * w02 * w03 * w12 * w13 * w13
    + 6 * w01 * w01 * w02 * w03 * w12 * w13 * w23 * w23
    - 12 * w01 * w01 * w02 * w03 * w12 * w13 * w23 + 6 * w01 * w01 * w02 * w03 * w12 * w13
    - w01 * w01 * w02 * w03 * w12 * w23 * w23 + 2 * w01 * w01 * w02 * w03 * w12 * w23
    - w01 * w01 * w02 * w03 * w12 + w01 * w01 * w02 * w03 * w13 * w13 * w23
    - w01 * w01 * w02 * w03 * w13 * w13 - w01 * w01 * w02 * w03 * w13 * w23 * w23
    + 2 * w01 * w01 * w02 * w03 * w13 * w23 - w01 * w01 * w02 * w03 * w13
    - w01 * w01 * w02 * w03 * w23 * w23 + w01 * w01 * w02 * w03 * w23
    + w01 * w01 * w02 * w12 * w12 * w13 * w23 * w23 - w01 * w01 * w02 * w12 * w12 * w13 * w23
    - w01 * w01 * w02 * w12 * w12 * w23 + w01 * w01 * w02 * w12 * w12
    - w01 * w01 * w02 * w12 * w13 * w23 * w23 + 2 * w01 * w01 * w02 * w12 * w13 * w23
    - w01 * w01 * w02 * w12 * w13 - w01 * w01 * w02 * w12 * w23 * w23
    + 2 * w01 * w01 * w02 * w12 * w23 - w01 * w01 * w02 * w12 - w01 * w01 * w02 * w13 * w23
    + w01 * w01 * w02 * w13 + w01 * w01 * w02 * w23 * w23 - w01 * w01 * w02 * w23
    - w01 * w01 * w03 * w03 * w12 * w13 * w13 * w23 + w01 * w01 * w03 * w03 * w12 * w13 * w13
    + w01 * w01 * w03 * w03 * w12 * w13 * w23 - w01 * w01 * w03 * w03 * w12 * w13
    + w01 * w01 * w03 * w03 * w13 * w13 * w23 - w01 * w01 * w03 * w03 * w13 * w13
    - w01 * w01 * w03 * w03 * w13 * w23 + w01 * w01 * w03 * w03 * w13
    + w01 * w01 * w03 * w12 * w13 * w13 * w23 * w23 - w01 * w01 * w03 * w12 * w13 * w13 * w23
    - w01 * w01 * w03 * w12 * w13 * w23 * w23 + 2 * w01 * w01 * w03 * w12 * w13 * w23
    - w01 * w01 * w03 * w12 * w13 - w01 * w01 * w03 * w12 * w23 + w01 * w01 * w03 * w12
    - w01 * w01 * w03 * w13 * w13 * w23 + w01 * w01 * w03 * w13 * w13
    - w01 * w01 * w03 * w13 * w23 * w23 + 2 * w01 * w01 * w03 * w13 * w23 - w01 * w01 * w03 * w13
    + w01 * w01 * w03 * w23 * w23 - w01 * w01 * w03 * w23 - w01 * w01 * w12 * w13 * w23 * w23
    + w01 * w01 * w12 * w13 * w23 + w01 * w01 * w12 * w23 * w23 - w01 * w01 * w12 * w23
    + w01 * w01 * w13 * w23 * w23 - w01 * w01 * w13 * w23 - w01 * w01 * w23 * w23
    + w01 * w01 * w23 + w01 * w02 * w02 * w03 * w03 * w12 * w12 * w13 * w13 * w23 * w23
    - w01 * w02 * w02 * w03 * w03 * w12 * w12 * w13 * w23 * w23
    - w01 * w02 * w02 * w03 * w03 * w12 * w12 * w13 * w23
    + w01 * w02 * w02 * w03 * w03 * w12 * w12 * w23
    - w01 * w02 * w02 * w03 * w03 * w12 * w13 * w13 * w23 * w23
    - w01 * w02 * w02 * w03 * w03 * w12 * w13 * w13 * w23
    + 2 * w01 * w02 * w02 * w03 * w03 * w12 * w13 * w23 * w23
    + 2 * w01 * w02 * w02 * w03 * w03 * w12 * w13 * w23
    - w01 * w02 * w02 * w03 * w03 * w12 * w23 * w23 - w01 * w02 * w02 * w03 * w03 * w12 * w23
    + w01 * w02 * w02 * w03 * w03 * w13 * w13 * w23
    - w01 * w02 * w02 * w03 * w03 * w13 * w23 * w23 - w01 * w02 * w02 * w03 * w03 * w13 * w23
    + w01 * w02 * w02 * w03 * w03 * w23 * w23
    - w01 * w02 * w02 * w03 * w12 * w12 * w13 * w13 * w23 * w23
    - w01 * w02 * w02 * w03 * w12 * w12 * w13 * w13 * w23
    + w01 * w02 * w02 * w03 * w12 * w12 * w13 * w23 * w23
    + 2 * w01 * w02 * w02 * w03 * w12 * w12 * w13 * w23 + w01 * w02 * w02 * w03 * w12 * w12 * w13
    - w01 * w02 * w02 * w03 * w12 * w12 * w23 - w01 * w02 * w02 * w03 * w12 * w12
    - w01 * w02 * w02 * w03 * w12 * w13 * w13 * w23 * w23
    + 6 * w01 * w02 * w02 * w03 * w12 * w13 * w13 * w23 - w01 * w02 * w02 * w03 * w12 * w13 * w13
    + 2 * w01 * w02 * w02 * w03 * w12 * w13 * w23 * w23
    - 12 * w01 * w02 * w02 * w03 * w12 * w13 * w23 + 2 * w01 * w02 * w02 * w03 * w12 * w13
    - w01 * w02 * w02 * w03 * w12 * w23 * w23 + 6 * w01 * w02 * w02 * w03 * w12 * w23
    - w01 * w02 * w02 * w03 * w12 - w01 * w02 * w02 * w03 * w13 * w13 * w23
    - w01 * w02 * w02 * w03 * w13 * w13 + w01 * w02 * w02 * w03 * w13 * w23 * w23
    + 2 * w01 * w02 * w02 * w03 * w13 * w23 + w01 * w02 * w02 * w03 * w13
    - w01 * w02 * w02 * w03 * w23 * w23 - w01 * w02 * w02 * w03 * w23
    + w01 * w02 * w02 * w12 * w12 * w13 * w13 * w23 - w01 * w02 * w02 * w12 * w12 * w13 * w23
    - w01 * w02 * w02 * w12 * w12 * w13 + w01 * w02 * w02 * w12 * w12
    - w01 * w02 * w02 * w12 * w13 * w13 * w23 - w01 * w02 * w02 * w12 * w13 * w13
    + 2 * w01 * w02 * w02 * w12 * w13 * w23 + 2 * w01 * w02 * w02 * w12 * w13
    - w01 * w02 * w02 * w12 * w23 - w01 * w02 * w02 * w12 + w01 * w02 * w02 * w13 * w13
    - w01 * w02 * w02 * w13 * w23 - w01 * w02 * w02 * w13 + w01 * w02 * w02 * w23
    - w01 * w02 * w03 * w03 * w12 * w12 * w13 * w13 * w23 * w23
    - w01 * w02 * w03 * w03 * w12 * w12 * w13 * w13 * w23
    - w01 * w02 * w03 * w03 * w12 * w12 * w13 * w23 * w23
    + 6 * w01 * w02 * w03 * w03 * w12 * w12 * w13 * w23 - w01 * w02 * w03 * w03 * w12 * w12 * w13
    - w01 * w02 * w03 * w03 * w12 * w12 * w23 - w01 * w02 * w03 * w03 * w12 * w12
    + w01 * w02 * w03 * w03 * w12 * w13 * w13 * w23 * w23
    + 2 * w01 * w02 * w03 * w03 * w12 * w13 * w13 * w23 + w01 * w02 * w03 * w03 * w12 * w13 * w13
    + 2 * w01 * w02 * w03 * w03 * w12 * w13 * w23 * w23
    - 12 * w01 * w02 * w03 * w03 * w12 * w13 * w23 + 2 * w01 * w02 * w03 * w03 * w12 * w13
    + w01 * w02 * w03 * w03 * w12 * w23 * w23 + 2 * w01 * w02 * w03 * w03 * w12 * w23
    + w01 * w02 * w03 * w03 * w12 - w01 * w02 * w03 * w03 * w13 * w13 * w23
    - w01 * w02 * w03 * w03 * w13 * w13 - w01 * w02 * w03 * w03 * w13 * w23 * w23
    + 6 * w01 * w02 * w03 * w03 * w13 * w23 - w01 * w02 * w03 * w03 * w13
    - w01 * w02 * w03 * w03 * w23 * w23 - w01 * w02 * w03 * w03 * w23
    + 2 * w01 * w02 * w03 * w12 * w12 * w13 * w13 * w23 * w23
    + 2 * w01 * w02 * w03 * w12 * w12 * w13 * w13 * w23
    + 2 * w01 * w02 * w03 * w12 * w12 * w13 * w23 * w23
    - 12 * w01 * w02 * w03 * w12 * w12 * w13 * w23 + 2 * w01 * w02 * w03 * w12 * w12 * w13
    + 2 * w01 * w02 * w03 * w12 * w12 * w23 + 2 * w01 * w02 * w03 * w12 * w12
    + 2 * w01 * w02 * w03 * w12 * w13 * w13 * w23 * w23
    - 12 * w01 * w02 * w03 * w12 * w13 * w13 * w23 + 2 * w01 * w02 * w03 * w12 * w13 * w13
    - 12 * w01 * w02 * w03 * w12 * w13 * w23 * w23 + 40 * w01 * w02 * w03 * w12 * w13 * w23
    - 12 * w01 * w02 * w03 * w12 * w13 + 2 * w01 * w02 * w03 * w12 * w23 * w23
    - 12 * w01 * w02 * w03 * w12 * w23 + 2 * w01 * w02 * w03 * w12
    + 2 * w01 * w02 * w03 * w13 * w13 * w23 + 2 * w01 * w02 * w03 * w13 * w13
    + 2 * w01 * w02 * w03 * w13 * w23 * w23 - 12 * w01 * w02 * w03 * w13 * w23
    + 2 * w01 * w02 * w03 * w13 + 2 * w01 * w02 * w03 * w23 * w23 + 2 * w01 * w02 * w03 * w23
    - w01 * w02 * w12 * w12 * w13 * w13 * w23 * w23 - w01 * w02 * w12 * w12 * w13 * w13 * w23
    - w01 * w02 * w12 * w12 * w13 * w23 * w23 + 6 * w01 * w02 * w12 * w12 * w13 * w23
    - w01 * w02 * w12 * w12 * w13 - w01 * w02 * w12 * w12 * w23 - w01 * w02 * w12 * w12
    + w01 * w02 * w12 * w13 * w13 * w23 * w23 + 2 * w01 * w02 * w12 * w13 * w13 * w23
    + w01 * w02 * w12 * w13 * w13 + 2 * w01 * w02 * w12 * w13 * w23 * w23
    - 12 * w01 * w02 * w12 * w13 * w23 + 2 * w01 * w02 * w12 * w13 + w01 * w02 * w12 * w23 * w23
    + 2 * w01 * w02 * w12 * w23 + w01 * w02 * w12 - w01 * w02 * w13 * w13 * w23
    - w01 * w02 * w13 * w13 - w01 * w02 * w13 * w23 * w23 + 6 * w01 * w02 * w13 * w23
    - w01 * w02 * w13 - w01 * w02 * w23 * w23 - w01 * w02 * w23
    + w01 * w03 * w03 * w12 * w12 * w13 * w13 * w23 - w01 * w03 * w03 * w12 * w12 * w13 * w23
    - w01 * w03 * w03 * w12 * w12 * w13 + w01 * w03 * w03 * w12 * w12
    - w01 * w03 * w03 * w12 * w13 * w13 * w23 - w01 * w03 * w03 * w12 * w13 * w13
    + 2 * w01 * w03 * w03 * w12 * w13 * w23 + 2 * w01 * w03 * w03 * w12 * w13
    - w01 * w03 * w03 * w12 * w23 - w01 * w03 * w03 * w12 + w01 * w03 * w03 * w13 * w13
    - w01 * w03 * w03 * w13 * w23 - w01 * w03 * w03 * w13 + w01 * w03 * w03 * w23
    - w01 * w03 * w12 * w12 * w13 * w13 * w23 * w23 - w01 * w03 * w12 * w12 * w13 * w13 * w23
    + w01 * w03 * w12 * w12 * w13 * w23 * w23 + 2 * w01 * w03 * w12 * w12 * w13 * w23
    + w01 * w03 * w12 * w12 * w13 - w01 * w03 * w12 * w12 * w23 - w01 * w03 * w12 * w12
    - w01 * w03 * w12 * w13 * w13 * w23 * w23 + 6 * w01 * w03 * w12 * w13 * w13 * w23
    - w01 * w03 * w12 * w13 * w13 + 2 * w01 * w03 * w12 * w13 * w23 * w23
    - 12 * w01 * w03 * w12 * w13 * w23 + 2 * w01 * w03 * w12 * w13 - w01 * w03 * w12 * w23 * w23
    + 6 * w01 * w03 * w12 * w23 - w01 * w03 * w12 - w01 * w03 * w13 * w13 * w23
    - w01 * w03 * w13 * w13 + w01 * w03 * w13 * w23 * w23 + 2 * w01 * w03 * w13 * w23
    + w01 * w03 * w13 - w01 * w03 * w23 * w23 - w01 * w03 * w23
    + w01 * w12 * w12 * w13 * w13 * w23 * w23 - w01 * w12 * w12 * w13 * w23 * w23
    - w01 * w12 * w12 * w13 * w23 + w01 * w12 * w12 * w23 - w01 * w12 * w13 * w13 * w23 * w23
    - w01 * w12 * w13 * w13 * w23 + 2 * w01 * w12 * w13 * w23 * w23 + 2 * w01 * w12 * w13 * w23
    - w01 * w12 * w23 * w23 - w01 * w12 * w23 + w01 * w13 * w13 * w23 - w01 * w13 * w23 * w23
    - w01 * w13 * w23 + w01 * w23 * w23 - w02 * w02 * w03 * w03 * w12 * w13 * w23 * w23
    + w02 * w02 * w03 * w03 * w12 * w13 * w23 + w02 * w02 * w03 * w03 * w12 * w23 * w23
    - w02 * w02 * w03 * w03 * w12 * w23 + w02 * w02 * w03 * w03 * w13 * w23 * w23
    - w02 * w02 * w03 * w03 * w13 * w23 - w02 * w02 * w03 * w03 * w23 * w23
    + w02 * w02 * w03 * w03 * w23 + w02 * w02 * w03 * w12 * w13 * w13 * w23 * w23
    - w02 * w02 * w03 * w12 * w13 * w13 * w23 - w02 * w02 * w03 * w12 * w13 * w23 * w23
    + 2 * w02 * w02 * w03 * w12 * w13 * w23 - w02 * w02 * w03 * w12 * w13
    - w02 * w02 * w03 * w12 * w23 + w02 * w02 * w03 * w12 - w02 * w02 * w03 * w13 * w13 * w23
    + w02 * w02 * w03 * w13 * w13 - w02 * w02 * w03 * w13 * w23 * w23
    + 2 * w02 * w02 * w03 * w13 * w23 - w02 * w02 * w03 * w13 + w02 * w02 * w03 * w23 * w23
    - w02 * w02 * w03 * w23 - w02 * w02 * w12 * w13 * w13 * w23 + w02 * w02 * w12 * w13 * w13
    + w02 * w02 * w12 * w13 * w23 - w02 * w02 * w12 * w13 + w02 * w02 * w13 * w13 * w23
    - w02 * w02 * w13 * w13 - w02 * w02 * w13 * w23 + w02 * w02 * w13
    + w02 * w03 * w03 * w12 * w12 * w13 * w23 * w23 - w02 * w03 * w03 * w12 * w12 * w13 * w23
    - w02 * w03 * w03 * w12 * w12 * w23 + w02 * w03 * w03 * w12 * w12
    - w02 * w03 * w03 * w12 * w13 * w23 * w23 + 2 * w02 * w03 * w03 * w12 * w13 * w23
    - w02 * w03 * w03 * w12 * w13 - w02 * w03 * w03 * w12 * w23 * w23
    + 2 * w02 * w03 * w03 * w12 * w23 - w02 * w03 * w03 * w12 - w02 * w03 * w03 * w13 * w23
    + w02 * w03 * w03 * w13 + w02 * w03 * w03 * w23 * w23 - w02 * w03 * w03 * w23
    - w02 * w03 * w12 * w12 * w13 * w13 * w23 * w23 + w02 * w03 * w12 * w12 * w13 * w13 * w23
    - w02 * w03 * w12 * w12 * w13 * w23 * w23 + 2 * w02 * w03 * w12 * w12 * w13 * w23
    - w02 * w03 * w12 * w12 * w13 + w02 * w03 * w12 * w12 * w23 - w02 * w03 * w12 * w12
    - w02 * w03 * w12 * w13 * w13 * w23 * w23 + 2 * w02 * w03 * w12 * w13 * w13 * w23
    - w02 * w03 * w12 * w13 * w13 + 6 * w02 * w03 * w12 * w13 * w23 * w23
    - 12 * w02 * w03 * w12 * w13 * w23 + 6 * w02 * w03 * w12 * w13 - w02 * w03 * w12 * w23 * w23
    + 2 * w02 * w03 * w12 * w23 - w02 * w03 * w12 + w02 * w03 * w13 * w13 * w23
    - w02 * w03 * w13 * w13 - w02 * w03 * w13 * w23 * w23 + 2 * w02 * w03 * w13 * w23
    - w02 * w03 * w13 - w02 * w03 * w23 * w23 + w02 * w03 * w23
    + w02 * w12 * w12 * w13 * w13 * w23 * w23 - w02 * w12 * w12 * w13 * w13 * w23
    - w02 * w12 * w12 * w13 * w23 + w02 * w12 * w12 * w13 - w02 * w12 * w13 * w13 * w23 * w23
    + 2 * w02 * w12 * w13 * w13 * w23 - w02 * w12 * w13 * w13 - w02 * w12 * w13 * w23 * w23
    + 2 * w02 * w12 * w13 * w23 - w02 * w12 * w13 - w02 * w13 * w13 * w23 + w02 * w13 * w13
    + w02 * w13 * w23 * w23 - w02 * w13 * w23 - w03 * w03 * w12 * w12 * w13 * w23
    + w03 * w03 * w12 * w12 * w13 + w03 * w03 * w12 * w12 * w23 - w03 * w03 * w12 * w12
    + w03 * w03 * w12 * w13 * w23 - w03 * w03 * w12 * w13 - w03 * w03 * w12 * w23
    + w03 * w03 * w12 + w03 * w12 * w12 * w13 * w13 * w23 * w23
    - w03 * w12 * w12 * w13 * w13 * w23 - w03 * w12 * w12 * w13 * w23 * w23
    + 2 * w03 * w12 * w12 * w13 * w23 - w03 * w12 * w12 * w13 - w03 * w12 * w12 * w23
    + w03 * w12 * w12 - w03 * w12 * w13 * w13 * w23 + w03 * w12 * w13 * w13
    - w03 * w12 * w13 * w23 * w23 + 2 * w03 * w12 * w13 * w23 - w03 * w12 * w13
    + w03 * w12 * w23 * w23 - w03 * w12 * w23 - w12 * w12 * w13 * w13 * w23 * w23
    + w12 * w12 * w13 * w13 * w23 + w12 * w12 * w13 * w23 * w23 - w12 * w12 * w13 * w23
    + w12 * w13 * w13 * w23 * w23 - w12 * w13 * w13 * w23 - w12 * w13 * w23 * w23
    + w12 * w13 * w23

/-- waves (0, 2), (1, 5), (3, −1), (7, 4): the weights -/
theorem config_one_weights :
    crossR 0 2 1 5 = -3/5 ∧ crossR 0 2 3 (-1) = 9 ∧ crossR 0 2 7 4 = 7/10 ∧
    crossR 1 5 3 (-1) = -3 ∧ crossR 1 5 7 4 = -1 ∧ crossR 3 (-1) 7 4 = 5/2 := by
  unfold crossR; refine ⟨?_, ?_, ?_, ?_, ?_, ?_⟩ <;> decide +kernel

/-- … and they obey the rule -/
theorem config_one_rule : rule4 (-3/5) 9 (7/10) (-3) (-1) (5/2) = 0 := by
  unfold rule4; decide +kernel

/-- waves (1, −3), (2, 6), (4, 0), (9, 5): the weights -/
theorem config_two_weights :
    crossR 1 (-3) 2 6 = 9/25 ∧ crossR 1 (-3) 4 0 = -9/7 ∧ crossR 1 (-3) 9 5 = 4/3 ∧
    crossR 2 6 4 0 = -3 ∧ crossR 2 6 9 5 = -7/9 ∧ crossR 4 0 9 5 = 25/9 := by
  unfold crossR; refine ⟨?_, ?_, ?_, ?_, ?_, ?_⟩ <;> decide +kernel

theorem config_two_rule : rule4 (9/25) (-9/7) (4/3) (-3) (-7/9) (25/9) = 0 := by
  unfold rule4; decide +kernel

/-- the rule is a real constraint: an array of all 2s breaks it, so no four waves on a line
    have every pair weight 2 -/
theorem rule4_all_two : rule4 2 2 2 2 2 2 = -32 := by
  unfold rule4; decide +kernel

/-- five weights fix the sixth up to two values: with the first five of configuration one,
    both 5/2 (the actual weight) and 92/77 satisfy the rule -/
theorem rule4_second_value : rule4 (-3/5) 9 (7/10) (-3) (-1) (92/77) = 0 := by
  unfold rule4; decide +kernel

end SeveralWaves
