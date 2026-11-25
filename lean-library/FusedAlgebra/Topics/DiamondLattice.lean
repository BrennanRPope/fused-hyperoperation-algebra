/-
  The diamond hierarchical lattice on the landscape
  =================================================

  Companion page: ../index.html.  Numbers: ../checks.py.
  Core Lean 4.37 only (no mathlib), no `sorry`, standard axioms only.
  Intended library path: FusedAlgebra/Topics/DiamondLattice.lean
  Imports FusedAlgebra.Topics.Statistics for the landscape law `lens k q a b`
  = (a + b + k·ab)/(1 − q·ab).

  In the transmissivity t = v/(Q + v) of a Potts bond (Q colours, v = e^K − 1),
  series squares t and parallel is the landscape law F_{Q−2, 1−Q}. One diamond step is
  T(t) = F_{Q−2,1−Q}(t², t²) = (2t² + (Q−2)t⁴)/(1 + (Q−1)t⁴).

  1. Walls and lift: the parallel law's walls are t = 1 (perfect bond) and t = 1/(1 − Q)
     (the colouring bond v = −1); the lift φ(t) = (1 − t)/(1 + (Q−1)t) = 1/(1 + v) turns it
     into multiplication; p² − 4q = Q².
  2. Fixed points: T(t) − t = −t(t − 1)((Q−1)t³ + t² + t − 1)/D. Writing Q = (y−1)²(y+1),
     the point t = 1/y is fixed, and the slope there is 4/(1 + t)² = (2y/(y+1))².
     Closed forms at Q = (m−1)²(m+1): Q = 3, 16, 45 give t = 1/2, 1/3, 1/4.
  3. The fold: 32/27 − (y−1)²(y+1) = (y + 1/3)²(5/3 − y); slope −1 at y = −1/3.
  4. Special Q: Q = 0 parabolic, Q = 1 polynomial, Q = 2 inversion symmetry, Q = 3 mirror
     (T(−1/2) = 1/2, and y = 2 is the only way the colouring wall mirrors the critical point).

  Sources: Berker–Ostlund 1979; Derrida–De Seze–Itzykson 1983; Tsallis–Levy 1981
  (transmissivity); Chang–Roeder–Shrock 2020 (fixed cubic, 32/27, q = 3 crossing).
  Ours: the landscape reading (walls, lift, Q² discriminant), the parametrisation
  Q = (y−1)²(y+1) with the slope 4/(1+t)², and the Q = 3 mirror explanation.
  LIBRARY NOTE: in the shared library at FusedAlgebra/Topics/DiamondLattice.lean; merged unchanged; no exact twins.
-/
import FusedAlgebra.Base
import FusedAlgebra.Topics.Statistics

namespace DiamondLattice
open MultAdd
open StatFused (lens)

local infixl:65 " ✳ " => multAdd

set_option linter.unusedSectionVars false

/-! # 0. The map -/

/-- One diamond step: series squares t, then the parallel landscape law F_{Q−2,1−Q}. -/
def dmap (Q t : Rat) : Rat := lens (Q - 2) (1 - Q) (t * t) (t * t)

def dnum (Q t : Rat) : Rat := 2 * t ^ 2 + (Q - 2) * t ^ 4
def dden (Q t : Rat) : Rat := 1 + (Q - 1) * t ^ 4

/-- The quotient-rule derivative of dnum/dden (checked against sympy in checks.py). -/
def dder (Q t : Rat) : Rat := 4 * t * (1 - t ^ 2) * (1 + (Q - 1) * t ^ 2) / (dden Q t) ^ 2

theorem dmap_formula (Q t : Rat) : dmap Q t = dnum Q t / dden Q t := by
  unfold dmap lens dnum dden
  have h1 : t * t + t * t + (Q - 2) * (t * t) * (t * t) = 2 * t ^ 2 + (Q - 2) * t ^ 4 := by grind
  have h2 : 1 - (1 - Q) * (t * t) * (t * t) = 1 + (Q - 1) * t ^ 4 := by grind
  rw [h1, h2]

/-- Series squares, so the step is even. -/
theorem dmap_even (Q t : Rat) : dmap Q (-t) = dmap Q t := by
  unfold dmap; rw [show -t * -t = t * t by grind]

/-- The quotient rule: N'D − ND' = 4t(1 − t²)(1 + (Q−1)t²). -/
theorem quotient_rule (Q t : Rat) :
    (4 * t + 4 * (Q - 2) * t ^ 3) * dden Q t - dnum Q t * (4 * (Q - 1) * t ^ 3)
      = 4 * t * (1 - t ^ 2) * (1 + (Q - 1) * t ^ 2) := by
  unfold dden dnum; grind

/-- Critical points besides 0, ±1: where t² = 1/(1−Q), i.e. 1 + (Q−1)t² = 0. -/
theorem crit_colour (Q t : Rat) (h : 1 + (Q - 1) * t ^ 2 = 0) :
    4 * t * (1 - t ^ 2) * (1 + (Q - 1) * t ^ 2) = 0 := by rw [h]; grind

/-! # 1. Walls and the lift -/

/-- The perfect bond t = 1 absorbs under the parallel law. -/
theorem wall_one (Q b : Rat) (h : 1 - (1 - Q) * 1 * b ≠ 0) : lens (Q - 2) (1 - Q) 1 b = 1 := by
  unfold lens at *; grind

/-- The colouring bond v = −1 sits at t = 1/(1 − Q). -/
theorem colour_bond (Q : Rat) (h : 1 - Q ≠ 0) : (-1) / (Q + -1) = 1 / (1 - Q) := by
  grind

/-- The second wall: t = 1/(1 − Q) absorbs under the parallel law. -/
theorem wall_colour (Q b : Rat) (hQ : 1 - Q ≠ 0) (h : 1 - b ≠ 0) :
    lens (Q - 2) (1 - Q) (1 / (1 - Q)) b = 1 / (1 - Q) := by
  unfold lens
  have hd : 1 - (1 - Q) * (1 / (1 - Q)) * b = 1 - b := by grind
  rw [hd]; grind

/-- The parallel law doubles its own colouring wall to itself: the critical value of T. -/
theorem colour_wall_doubles (Q : Rat) (hQ : 1 - Q ≠ 0) (h : Q ≠ 0) :
    lens (Q - 2) (1 - Q) (1 / (1 - Q)) (1 / (1 - Q)) = 1 / (1 - Q) :=
  wall_colour Q (1 / (1 - Q)) hQ (by grind)

theorem quot_helper (x1 x2 y1 y2 D : Rat) (hD : D ≠ 0) (h1 : y1 ≠ 0) (h2 : y2 ≠ 0) :
    x1 * x2 / D / (y1 * y2 / D) = x1 / y1 * (x2 / y2) := by grind

theorem slope_helper (u w : Rat) (hu : u ≠ 0) (hw : w ≠ 0) :
    4 * w ^ 2 * u ^ 2 / (u ^ 2 * w) ^ 2 = 4 / u ^ 2 := by grind

/-- The lift of the parallel law. -/
def lift (Q t : Rat) : Rat := (1 - t) / (1 + (Q - 1) * t)

/-- **The lift turns parallel into multiplication.** -/
theorem lift_mul (Q a b : Rat) (h : 1 - (1 - Q) * a * b ≠ 0) (ha : 1 + (Q - 1) * a ≠ 0)
    (hb : 1 + (Q - 1) * b ≠ 0) :
    lift Q (lens (Q - 2) (1 - Q) a b) = lift Q a * lift Q b := by
  unfold lift lens
  have hnum : 1 - (a + b + (Q - 2) * a * b) / (1 - (1 - Q) * a * b)
      = (1 - a) * (1 - b) / (1 - (1 - Q) * a * b) := by grind
  have hden : 1 + (Q - 1) * ((a + b + (Q - 2) * a * b) / (1 - (1 - Q) * a * b))
      = (1 + (Q - 1) * a) * (1 + (Q - 1) * b) / (1 - (1 - Q) * a * b) := by grind
  rw [hnum, hden]; exact quot_helper _ _ _ _ _ h ha hb

/-- The lift of a bond's transmissivity is 1/(1 + v) = e^{−K}. -/
theorem lift_weight (Q v : Rat) (h : Q + v ≠ 0) (h1 : 1 + v ≠ 0) (hQ : Q ≠ 0) :
    lift Q (v / (Q + v)) = 1 / (1 + v) := by
  unfold lift
  have hn : 1 - v / (Q + v) = Q / (Q + v) := by grind
  have hd : 1 + (Q - 1) * (v / (Q + v)) = Q * (1 + v) / (Q + v) := by grind
  rw [hn, hd]; grind

/-- So one diamond step squares in t, then squares in the lift. -/
theorem lift_dmap (Q t : Rat) (h : 1 - (1 - Q) * (t * t) * (t * t) ≠ 0)
    (ha : 1 + (Q - 1) * (t * t) ≠ 0) : lift Q (dmap Q t) = lift Q (t * t) * lift Q (t * t) :=
  lift_mul Q _ _ h ha ha

/-- **The bend's discriminant is a square**: p² − 4q = Q² on the "or" line. -/
theorem disc_square (Q : Rat) : (Q - 2) ^ 2 - 4 * (1 - Q) = Q ^ 2 := by grind

/-- The two walls 1 and 1/(1−Q) differ by −Q/(1−Q): they merge only at Q = 0. -/
theorem wall_gap (Q : Rat) (h : 1 - Q ≠ 0) : 1 - 1 / (1 - Q) = -Q / (1 - Q) := by grind

/-- The bend is the reciprocal of the second wall. -/
theorem bend_recip (Q : Rat) : 1 / (1 / (1 - Q)) = 1 - Q := by grind

/-! # 2. Fixed points: one cubic, one parametrisation, one slope -/

/-- **T(t) − t factors**: the nontrivial fixed points solve (Q−1)t³ + t² + t − 1 = 0. -/
theorem fixed_factor (Q t : Rat) :
    dnum Q t - t * dden Q t = -t * (t - 1) * ((Q - 1) * t ^ 3 + t ^ 2 + t - 1) := by
  unfold dnum dden; grind

/-- Writing Q = (y−1)²(y+1) solves the cubic in y. -/
theorem param_cubic (y Q : Rat) (hQ : Q = (y - 1) ^ 2 * (y + 1)) :
    (Q - 1) + y + y ^ 2 - y ^ 3 = 0 := by subst hQ; grind

/-- ... and t = 1/y solves the fixed-point cubic. -/
theorem param_fixed (y Q : Rat) (hy : y ≠ 0) (hQ : Q = (y - 1) ^ 2 * (y + 1)) :
    (Q - 1) * (1 / y) ^ 3 + (1 / y) ^ 2 + 1 / y - 1 = 0 := by
  have hc := param_cubic y Q hQ
  have e : (Q - 1) * (1 / y) ^ 3 + (1 / y) ^ 2 + 1 / y - 1 = ((Q - 1) + y + y ^ 2 - y ^ 3) / y ^ 3 := by
    grind
  rw [e, hc]; grind

/-- **Every nontrivial fixed point is t = 1/y with Q = (y−1)²(y+1).** -/
theorem dmap_fixed_param (y Q : Rat) (hy : y ≠ 0) (hQ : Q = (y - 1) ^ 2 * (y + 1))
    (hD : dden Q (1 / y) ≠ 0) : dmap Q (1 / y) = 1 / y := by
  rw [dmap_formula]
  have hf := fixed_factor Q (1 / y)
  rw [param_fixed y Q hy hQ] at hf
  have : dnum Q (1 / y) = 1 / y * dden Q (1 / y) := by grind
  rw [this]; grind

/-- The fixed point in the Chang–Roeder–Shrock form q² + 2qv − v³ = 0, v = y² − 1. -/
theorem crs_cubic (y : Rat) :
    ((y - 1) ^ 2 * (y + 1)) ^ 2 + 2 * ((y - 1) ^ 2 * (y + 1)) * (y ^ 2 - 1) - (y ^ 2 - 1) ^ 3 = 0 := by
  grind

/-- At a nontrivial fixed point the denominator factors: D = (1 + t)²(1 − t). -/
theorem dden_at_fixed (Q t : Rat) (h : (Q - 1) * t ^ 3 + t ^ 2 + t - 1 = 0) :
    dden Q t = (1 + t) ^ 2 * (1 - t) := by
  unfold dden
  have : (Q - 1) * t ^ 4 = t * ((Q - 1) * t ^ 3) := by grind
  rw [this, show (Q - 1) * t ^ 3 = 1 - t - t ^ 2 by grind]; grind

theorem slope_num_at_fixed (Q t : Rat) (h : (Q - 1) * t ^ 3 + t ^ 2 + t - 1 = 0) :
    t * (4 * t * (1 - t ^ 2) * (1 + (Q - 1) * t ^ 2)) = t * (4 * (1 - t) ^ 2 * (1 + t) ^ 2) := by
  have e : t * (1 + (Q - 1) * t ^ 2) = 1 - t ^ 2 := by
    have : t * (1 + (Q - 1) * t ^ 2) = t + (Q - 1) * t ^ 3 := by grind
    rw [this, show (Q - 1) * t ^ 3 = 1 - t - t ^ 2 by grind]; grind
  have : t * (4 * t * (1 - t ^ 2) * (1 + (Q - 1) * t ^ 2)) = 4 * t * (1 - t ^ 2) * (t * (1 + (Q - 1) * t ^ 2)) := by
    grind
  rw [this, e]; grind

/-- **The slope at every nontrivial fixed point is 4/(1 + t)².** -/
theorem slope_at_fixed (Q t : Rat) (h : (Q - 1) * t ^ 3 + t ^ 2 + t - 1 = 0)
    (ht : t ≠ 0) (h1 : 1 + t ≠ 0) (h2 : 1 - t ≠ 0) : dder Q t = 4 / (1 + t) ^ 2 := by
  unfold dder
  have hn := slope_num_at_fixed Q t h
  have hn' : 4 * t * (1 - t ^ 2) * (1 + (Q - 1) * t ^ 2) = 4 * (1 - t) ^ 2 * (1 + t) ^ 2 := by
    grind
  rw [hn', dden_at_fixed Q t h]
  exact slope_helper _ _ h1 h2

/-- In y: 4/(1 + 1/y)² = (2y/(y + 1))². -/
theorem slope_in_y (y : Rat) (hy : y ≠ 0) (h1 : y + 1 ≠ 0) :
    4 / (1 + 1 / y) ^ 2 = (2 * y / (y + 1)) ^ 2 := by
  have e : 1 + 1 / y = (y + 1) / y := by grind
  rw [e]; grind

/-! # 3. Closed forms: Q = (m−1)²(m+1) gives t* = 1/m -/

theorem three_fixed : dmap 3 (1 / 2) = 1 / 2 := by unfold dmap lens; grind
theorem three_slope : dder 3 (1 / 2) = 16 / 9 := by unfold dder dden; grind
theorem sixteen_fixed : dmap 16 (1 / 3) = 1 / 3 := by unfold dmap lens; grind
theorem sixteen_slope : dder 16 (1 / 3) = 9 / 4 := by unfold dder dden; grind
theorem fortyfive_fixed : dmap 45 (1 / 4) = 1 / 4 := by unfold dmap lens; grind
theorem fortyfive_slope : dder 45 (1 / 4) = 64 / 25 := by unfold dder dden; grind

/-- Percolation (Q = 1): at the golden point t² = 1 − t the slope is 4 − 4t. -/
theorem golden_slope (t : Rat) (h : t * t = 1 - t) (ht : t ≠ 0) (h1 : 1 + t ≠ 0) (h2 : 1 - t ≠ 0) :
    dder 1 t = 4 - 4 * t := by
  rw [slope_at_fixed 1 t (by grind) ht h1 h2]
  have : (1 + t) ^ 2 * (4 - 4 * t) = 4 := by grind
  grind

/-- Ising (Q = 2, tribonacci point t³ + t² + t = 1): slope 4/(1 + t)². -/
theorem ising_slope (t : Rat) (h : t ^ 3 + t ^ 2 + t = 1) (ht : t ≠ 0) (h1 : 1 + t ≠ 0)
    (h2 : 1 - t ≠ 0) : dder 2 t = 4 / (1 + t) ^ 2 :=
  slope_at_fixed 2 t (by grind) ht h1 h2

/-! # 4. The fold at 32/27 -/

/-- **The fold**: Q(y) = (y−1)²(y+1) touches 32/27 at y = −1/3. -/
theorem fold_identity (y : Rat) :
    32 / 27 - (y - 1) ^ 2 * (y + 1) = (y + 1 / 3) ^ 2 * (5 / 3 - y) := by grind

theorem fold_value : ((-1 / 3 : Rat) - 1) ^ 2 * ((-1 / 3) + 1) = 32 / 27 := by grind

/-- At Q = 32/27 the ferromagnetic point is t* = 3/5 (y = 5/3), slope 25/16. -/
theorem fold_ferro_fixed : dmap (32 / 27) (3 / 5) = 3 / 5 := by unfold dmap lens; grind
theorem fold_ferro_slope : dder (32 / 27) (3 / 5) = 25 / 16 := by unfold dder dden; grind

/-- At the fold the slope factor 2y/(y + 1) is −1: neutral. -/
theorem fold_neutral : 2 * (-1 / 3 : Rat) / ((-1 / 3) + 1) = -1 := by grind

/-- Attracting exactly when (1 − y)(1 + 3y) > 0: (y + 1)² − (2y)² = (1 − y)(1 + 3y). -/
theorem attract_identity (y : Rat) : (y + 1) ^ 2 - (2 * y) ^ 2 = (1 - y) * (1 + 3 * y) := by grind

/-- The fold point in t and in x = e^K: t = 1/y = −3 and x = y² = 1/9. -/
theorem fold_point : 1 / (-1 / 3 : Rat) = -3 ∧ (-1 / 3 : Rat) ^ 2 = 1 / 9 := by
  constructor <;> grind

/-! # 5. Special Q along the "or" line -/

/-- Q = 0, where the two walls merge (p² = 4q): T = 2t²/(1 + t²) ... -/
theorem q_zero_map (t : Rat) (h : 1 - t ^ 2 ≠ 0) (h' : 1 + t ^ 2 ≠ 0) :
    dmap 0 t = 2 * t ^ 2 / (1 + t ^ 2) := by
  rw [dmap_formula]; unfold dnum dden
  have e1 : 2 * t ^ 2 + (0 - 2) * t ^ 4 = (2 * t ^ 2) * (1 - t ^ 2) := by grind
  have e2 : 1 + (0 - 1) * t ^ 4 = (1 + t ^ 2) * (1 - t ^ 2) := by grind
  rw [e1, e2]; grind

/-- ... whose fixed point t = 1 is double: T(t) − t = −t(t − 1)²/(1 + t²) ... -/
theorem q_zero_double (t : Rat) (h' : 1 + t ^ 2 ≠ 0) :
    2 * t ^ 2 / (1 + t ^ 2) - t = -t * (t - 1) ^ 2 / (1 + t ^ 2) := by grind

/-- ... and parabolic: the derivative 4t/(1 + t²)² of 2t²/(1 + t²) is 1 at t = 1. -/
theorem q_zero_parabolic : 4 * (1 : Rat) / (1 + 1 ^ 2) ^ 2 = 1 := by grind

/-- Q = 1: the second wall moves to infinity, the series wall, and T is a polynomial. -/
theorem q_one_poly (t : Rat) : dmap 1 t = 2 * t ^ 2 - t ^ 4 := by
  rw [dmap_formula]; unfold dnum dden; grind

/-- Q = 2: walls 0, ∞, 1, −1 are harmonic and T(1/t) = T(t). -/
theorem q_two_inversion (t : Rat) (ht : t ≠ 0) (h : 1 + t ^ 4 ≠ 0) : dmap 2 (1 / t) = dmap 2 t := by
  rw [dmap_formula, dmap_formula]; unfold dnum dden
  have h4 : t ^ 4 ≠ 0 := by grind
  have e1 : 2 * (1 / t) ^ 2 + (2 - 2) * (1 / t) ^ 4 = (2 * t ^ 2) / t ^ 4 := by grind
  have e2 : 1 + (2 - 1) * (1 / t) ^ 4 = (1 + t ^ 4) / t ^ 4 := by grind
  rw [e1, e2]; grind

/-- Q = 3: the colouring wall −1/2 is the mirror of the critical point 1/2, so it lands on it. -/
theorem q_three_mirror : dmap 3 (-1 / 2) = 1 / 2 := by unfold dmap lens; grind

/-- The colouring wall mirrors the critical point (1/(1−Q) = −1/y) only for y = 2, i.e. Q = 3. -/
theorem mirror_only_three (y : Rat) (h : (y - 1) ^ 2 * (y + 1) = y + 1) :
    y = -1 ∨ y = 0 ∨ y = 2 := by
  have : (y + 1) * (y * (y - 2)) = 0 := by grind
  grind

/-- The mult-add tie: parallel joins bond weights by mult-add, so the lift 1 + v multiplies. -/
theorem parallel_weights (v w : Rat) : 1 + (v ✳ w) = (1 + v) * (1 + w) := multAdd_lift v w

end DiamondLattice
