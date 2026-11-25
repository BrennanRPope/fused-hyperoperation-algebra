/-
  Fused Hyperoperation Algebra: knot theory (tangle sums are the dial)
  ====================================================================
  Core Lean 4 only (no mathlib), no `sorry`. Builds on FusedAlgebra.Topics.QAcross.Knots
  (Kauffman bracket from PD codes) and QAcross.QDial (q-integers as the dial).

  A 2-tangle T (a knot diagram cut open at four ends) has Kauffman bracket
  ⟨T⟩ = f·⟨[0]⟩ + g·⟨[∞]⟩ ([0] = two horizontal arcs, [∞] = two vertical arcs), stored as (f, g).
  Standard (Kauffman 1987, Goldman–Kauffman 1997): the sum of tangles side by side is
  (f₁f₂, f₁g₂ + g₁f₂ + d·g₁g₂), d = −A² − A⁻² the loop value.

  What this file shows, through the project's words:
    * in x = g/f the sum is the dial x₁ ◯_d x₂ with weight d (ratio_dial);
    * the lift 1 + d·x is closing the tangle (D joins left ends and right ends), and the lift
      law is ⟨D(S + T)⟩ = ⟨D(S)⟩·⟨D(T)⟩: connected sum (closeD_tsum);
    * in y = A²·g/f one crossing is the q-dial step y ↦ q·y + 1 with q = −A⁻⁴ = −t, so n twists
      give y = [n]_q (twist_formula, twist_step_y);
    * the twist's two walls merge exactly when A⁴ = −1, where d = 0, q = 1 and y is Conway's
      fraction, which adds (conway_twist, conway_add).
  LIBRARY NOTE: in the shared library at FusedAlgebra/Topics/KnotTheory.lean; one exact twin of QDial removed (see the LIBRARY comment).
-/
import FusedAlgebra.Base
import FusedAlgebra.Topics.QAcross.QDial
import FusedAlgebra.Topics.QAcross.Knots

namespace KnotTheory
open CircleOp
open Lean.Grind

open scoped CircleOp

set_option linter.unusedSectionVars false

/-! ## 1. The tangle algebra on any commutative ring -/

section Ring
variable {R : Type} [CommRing R]

/-- Side-by-side sum of tangles (f₁, g₁) + (f₂, g₂). -/
def tsum (d : R) (S T : R × R) : R × R := (S.1 * T.1, S.1 * T.2 + S.2 * T.1 + d * S.2 * T.2)
/-- Stacked sum (S above T). -/
def vsum (d : R) (S T : R × R) : R × R := (S.1 * T.2 + S.2 * T.1 + d * S.1 * T.1, S.2 * T.2)
/-- Denominator closure ⟨D(T)⟩: join the two left ends and the two right ends. -/
def closeD (d : R) (T : R × R) : R := T.1 + d * T.2
/-- Numerator closure ⟨N(T)⟩: join the two top ends and the two bottom ends. -/
def closeN (d : R) (T : R × R) : R := d * T.1 + T.2
def zeroT : R × R := (1, 0)
def infT : R × R := (0, 1)
/-- Rotating a tangle a quarter turn swaps [0] and [∞]. -/
def swapT (T : R × R) : R × R := (T.2, T.1)

theorem tsum_comm (d : R) (S T : R × R) : tsum d S T = tsum d T S := by
  unfold tsum; congr 1 <;> grind

theorem tsum_assoc (d : R) (S T U : R × R) : tsum d (tsum d S T) U = tsum d S (tsum d T U) := by
  unfold tsum; congr 1 <;> grind

/-- [0] is the identity of tangle addition. -/
theorem tsum_zeroT (d : R) (S : R × R) : tsum d S zeroT = S := by
  unfold tsum zeroT; cases S; congr 1 <;> grind

/-- The lift law is connected sum: closing a sum gives the product of the closures. -/
theorem closeD_tsum (d : R) (S T : R × R) : closeD d (tsum d S T) = closeD d S * closeD d T := by
  unfold closeD tsum; grind

/-- The homogeneous form of the dial: the [∞]-part of a sum, times d, plus its [0]-part, is the
    product (f₁ + d g₁)(f₂ + d g₂): the weighted lift 1 + d·x with x = g/f, cleared of fractions. -/
theorem tsum_lift (d : R) (S T : R × R) :
    (tsum d S T).1 + d * (tsum d S T).2 = (S.1 + d * S.2) * (T.1 + d * T.2) := by
  unfold tsum; grind

/-- The "mult-add value" P − 1 of a closed diagram turns connected sum into ◯ exactly.
    (A relabelling: it is the lift of P(K # L) = P(K)·P(L).) -/
theorem closeD_circ (d : R) (S T : R × R) :
    closeD d (tsum d S T) - 1 = (closeD d S - 1) ◯ (closeD d T - 1) := by
  rw [closeD_tsum]; unfold circ; grind

/-- The same for any multiplicative knot polynomial (Jones, Alexander, HOMFLY). -/
theorem poly_circ (P Q : R) : P * Q - 1 = (P - 1) ◯ (Q - 1) := by
  unfold circ; grind

/-- The wall: a tangle whose closure has bracket 0 absorbs everything added to it. -/
theorem closeD_wall (d : R) (S T : R × R) (h : closeD d S = 0) : closeD d (tsum d S T) = 0 := by
  rw [closeD_tsum, h]; grind

/-- The second wall, x = ∞: adding [∞] leaves only [∞], carrying the closure ⟨D(S)⟩. -/
theorem tsum_infT (d : R) (S : R × R) : tsum d S infT = (0, closeD d S) := by
  unfold tsum infT closeD; congr 1 <;> grind

/-- Stacking is adding, turned a quarter: vsum = swap ∘ tsum ∘ swap. -/
theorem vsum_eq_swap (d : R) (S T : R × R) : vsum d S T = swapT (tsum d (swapT S) (swapT T)) := by
  unfold vsum tsum swapT; congr 1 <;> grind

/-- So the numerator closure is the lift for stacking. -/
theorem closeN_vsum (d : R) (S T : R × R) : closeN d (vsum d S T) = closeN d S * closeN d T := by
  unfold closeN vsum; grind

/-- A loop on its own: [0] stacked on [0] is d·[0]. -/
theorem stack_zero_zero (d : R) : vsum d (zeroT : R × R) zeroT = (d, 0) := by
  unfold vsum zeroT; congr 1 <;> grind

/-! ## 2. One crossing, n crossings: the q-dial (A·Ai = 1, Ai = A⁻¹) -/

/-- The loop value d = −A² − A⁻². -/
def loopv (A Ai : R) : R := -(A * A) - Ai * Ai
/-- The tangle q: q = −A⁻⁴ (= −t, t the Jones variable). -/
def tq (Ai : R) : R := -(Ai * Ai * Ai * Ai)

/-- n crossings side by side: the twist tangle [n]. -/
def twist (A Ai : R) : Nat → R × R
  | 0 => zeroT
  | n + 1 => tsum (loopv A Ai) (twist A Ai n) (A, Ai)

/-- In y = A²·x the dial weight d becomes d·A⁻² = q − 1. -/
theorem weight_is_q_minus_one (A Ai : R) (h : A * Ai = 1) : loopv A Ai * (Ai * Ai) = tq Ai - 1 := by
  unfold loopv tq; grind

/-- Rescaling by A² turns the d-dial into the (q − 1)-dial. -/
theorem y_scale (A Ai x y : R) (h : A * Ai = 1) :
    A * A * circK (loopv A Ai) x y = circK (tq Ai - 1) (A * A * x) (A * A * y) := by
  unfold circK loopv tq; grind

/-- n twists: [n] = Aⁿ·[0] + Aⁿ·A⁻²·[n]_q·[∞]. -/
theorem twist_formula (A Ai : R) (h : A * Ai = 1) (n : Nat) :
    twist A Ai n = (A ^ n, A ^ n * (Ai * Ai) * QDial.qnum (tq Ai) n) := by
  induction n with
  | zero => show ((1 : R), (0 : R)) = _; rw [Semiring.pow_zero, QDial.qnum_zero]; congr 1; grind
  | succ n ih =>
    rw [twist, ih, QDial.qnum_succ, Semiring.pow_succ]
    unfold tsum loopv tq; congr 1 <;> grind

/-- One crossing closed up is a kink: ⟨D(X)⟩ = −A⁻³. -/
theorem kink (A Ai : R) (h : A * Ai = 1) : closeD (loopv A Ai) (A, Ai) = -(Ai * Ai * Ai) := by
  unfold closeD loopv; grind

/-- Closing [n] gives n kinks: (−A⁻³)ⁿ. -/
theorem twist_closeD (A Ai : R) (h : A * Ai = 1) (n : Nat) :
    closeD (loopv A Ai) (twist A Ai n) = (-(Ai * Ai * Ai)) ^ n := by
  induction n with
  | zero => rw [Semiring.pow_zero]; show (1 : R) + _ * 0 = 1; grind
  | succ n ih => rw [twist, closeD_tsum, ih, kink A Ai h, Semiring.pow_succ]

/-- The twist matrix [[A, 0], [A⁻¹, −A⁻³]] has two eigenvalues, A and −A⁻³ (the two walls);
    they coincide exactly when A⁴ = −1. -/
theorem walls_merge (A Ai : R) (h : A * Ai = 1) : A = -(Ai * Ai * Ai) ↔ A * A * A * A = -1 := by
  constructor
  · intro e; have : A * A * A * A = A * A * A * (-(Ai * Ai * Ai)) := by rw [← e]
    grind
  · intro e; have : -(Ai * Ai * Ai) = -(Ai * Ai * Ai) * (-(A * A * A * A)) := by rw [e]; grind
    grind

/-- Where the walls merge the loop value is 0 ... -/
theorem loop_zero (A Ai : R) (h : A * Ai = 1) (e : A * A * A * A = -1) : loopv A Ai = 0 := by
  unfold loopv
  have : Ai * Ai = Ai * Ai * (-(A * A * A * A)) := by rw [e]; grind
  grind

/-- ... and the tangle q is 1 (the counting point of the dial). -/
theorem q_one (A Ai : R) (h : A * Ai = 1) (e : A * A * A * A = -1) : tq Ai = 1 := by
  unfold tq
  have : Ai * Ai * Ai * Ai = Ai * Ai * Ai * Ai * (-(A * A * A * A)) := by rw [e]; grind
  grind

/-- Mirror image A ↦ A⁻¹ sends q to 1/q: q(A)·q(A⁻¹) = 1. -/
theorem mirror_q (A Ai : R) (h : A * Ai = 1) : tq Ai * tq A = 1 := by
  unfold tq; grind

/-- The Potts number of the loop: Q = d² = t + 2 + t⁻¹ with t = A⁻⁴. -/
theorem potts_Q (A Ai : R) (h : A * Ai = 1) :
    loopv A Ai * loopv A Ai = Ai * Ai * Ai * Ai + 2 + A * A * A * A := by
  unfold loopv; grind

/-! ## 3. Torus knots: Alexander polynomials as q-numbers -/

/-- [pq]_t = [p]_t·[q]_(t^p) (the q-number behind Δ(T(p,q)) = [q]_(t^p)/[q]_t). -/
theorem qnum_mul (t : R) (p n : Nat) :
    QDial.qnum t (p * n) = QDial.qnum t p * QDial.qnum (t ^ p) n := by
  induction n with
  | zero => rw [Nat.mul_zero, QDial.qnum_zero, QDial.qnum_zero]; grind
  | succ n ih =>
    rw [Nat.mul_succ, QDial.qnum_add_std, ih, QDial.qnum_add_std (t ^ p) n 1]
    have e : (t ^ p) ^ n = t ^ (p * n) := by
      rw [← pw_eq_pow, ← pw_eq_pow, ← pw_eq_pow, pw_mul]
    rw [e]
    show _ = QDial.qnum t p * (QDial.qnum (t ^ p) n + t ^ (p * n) * (1 + t ^ p * 0))
    grind

/-- T(2, n): (1 + t)·[n]_(−t) = 1 − (−t)ⁿ, so for odd n, Δ(T(2,n)) = [n]_(−t). -/
theorem alexander_T2 (t : R) (n : Nat) : (1 + t) * QDial.qnum (-t) n = 1 - (-t) ^ n := by
  have := QDial.qnum_lift (-t) n
  grind

-- LIBRARY: qnum_at_one (at q = 1 a q-number is n) was an exact twin of QDial.qnum_at_one and is removed.

end Ring

/-! ## 4. Ratios: the dial itself, Conway's fraction (fields) -/

section Field
variable {K : Type} [Field K]

/-- In x = g/f, adding tangles is the weighted mult-add with weight d. -/
theorem ratio_dial (d : K) (S T : K × K) (hS : S.1 ≠ 0) (hT : T.1 ≠ 0) :
    (tsum d S T).2 / (tsum d S T).1 = circK d (S.2 / S.1) (T.2 / T.1) := by
  unfold tsum circK; grind

/-- The lift 1 + d·x is the closure, divided by the [0]-coefficient. -/
theorem lift_is_closing (d : K) (T : K × K) (hT : T.1 ≠ 0) :
    1 + d * (T.2 / T.1) = closeD d T / T.1 := by
  unfold closeD; grind

/-- Stacking is the same dial in z = f/g. -/
theorem stack_dial (d : K) (S T : K × K) (hS : S.2 ≠ 0) (hT : T.2 ≠ 0) :
    (vsum d S T).1 / (vsum d S T).2 = circK d (S.1 / S.2) (T.1 / T.2) := by
  unfold vsum circK; grind

/-- One crossing added: y ↦ q·y + 1, with y = A²·g/f. -/
theorem twist_step_y (A Ai : K) (h : A * Ai = 1) (S : K × K) (hS : S.1 ≠ 0) :
    A * A * (tsum (loopv A Ai) S (A, Ai)).2 / (tsum (loopv A Ai) S (A, Ai)).1 =
      tq Ai * (A * A * S.2 / S.1) + 1 := by
  have hA : A ≠ 0 := by intro e; rw [e] at h; grind
  unfold tsum loopv tq; grind

/-- One crossing stacked: 1/y ↦ q⁻¹·(1/y) + 1 (the same step at 1/q, seen through 1/y). -/
theorem stack_step_y (A Ai : K) (h : A * Ai = 1) (S : K × K) (hS : S.2 ≠ 0) :
    (vsum (loopv A Ai) S (A, Ai)).1 / (A * A * (vsum (loopv A Ai) S (A, Ai)).2) =
      tq A * (S.1 / (A * A * S.2)) + 1 := by
  have hA : A ≠ 0 := by intro e; rw [e] at h; grind
  have hAi : Ai ≠ 0 := by intro e; rw [e] at h; grind
  unfold vsum loopv tq; grind

/-- The q-step is a rotation (|q| = 1) or scaling about its finite wall c = 1/(1 − q). -/
theorem step_about_wall (q y : K) (hq : 1 - q ≠ 0) : q * y + 1 - 1 / (1 - q) = q * (y - 1 / (1 - q)) := by
  grind

/-- So the twists sit on a circle around the wall: [n]_q = c − qⁿ·c. -/
theorem twists_on_circle (q : K) (hq : 1 - q ≠ 0) (n : Nat) :
    QDial.qnum q n = 1 / (1 - q) - q ^ n * (1 / (1 - q)) := by
  have := QDial.qnum_lift q n
  grind

/-- At a root of unity qᵐ = 1 the polygon closes: m twists return y to 0. -/
theorem twists_close (q : K) (hq : 1 - q ≠ 0) (m : Nat) (hm : q ^ m = 1) : QDial.qnum q m = 0 := by
  rw [twists_on_circle q hq m, hm]; grind

/-- Conway's fraction F = A²·g/f at A⁴ = −1: a crossing adds 1. -/
theorem conway_twist (A Ai : K) (h : A * Ai = 1) (e : A * A * A * A = -1) (S : K × K) (hS : S.1 ≠ 0) :
    A * A * (tsum (loopv A Ai) S (A, Ai)).2 / (tsum (loopv A Ai) S (A, Ai)).1 = A * A * S.2 / S.1 + 1 := by
  rw [twist_step_y A Ai h S hS, q_one A Ai h e]; grind

/-- Conway's fraction is additive under tangle sum (the dial at weight d = 0). -/
theorem conway_add (A Ai : K) (h : A * Ai = 1) (e : A * A * A * A = -1) (S T : K × K)
    (hS : S.1 ≠ 0) (hT : T.1 ≠ 0) :
    A * A * (tsum (loopv A Ai) S T).2 / (tsum (loopv A Ai) S T).1 =
      A * A * S.2 / S.1 + A * A * T.2 / T.1 := by
  have h0 := loop_zero A Ai h e
  rw [h0]; unfold tsum; grind

/-- ... and 1/F is additive under stacking. -/
theorem conway_stack (A Ai : K) (h : A * Ai = 1) (e : A * A * A * A = -1) (S T : K × K)
    (hS : S.2 ≠ 0) (hT : T.2 ≠ 0) :
    (vsum (loopv A Ai) S T).1 / (A * A * (vsum (loopv A Ai) S T).2) =
      S.1 / (A * A * S.2) + T.1 / (A * A * T.2) := by
  have hA : A ≠ 0 := by intro e'; rw [e'] at h; grind
  have h0 := loop_zero A Ai h e
  rw [h0]; unfold vsum; grind

end Field

/-! ## 5. Integer fractions: continued fractions and SL(2, ℤ) -/

/-- Horizontal twists: p/q ↦ (p + a q)/q. -/
def hstep (a : Int) (pq : Int × Int) : Int × Int := (pq.1 + a * pq.2, pq.2)
/-- Vertical twists: 1/F ↦ 1/F + a, i.e. p/q ↦ p/(q + a p). -/
def vstep (a : Int) (pq : Int × Int) : Int × Int := (pq.1, pq.2 + a * pq.1)

/-- h^a₁ v^a₂ h^a₃ … applied to [0] = 0/1. -/
def cfracAux : Bool → List Int → Int × Int → Int × Int
  | _, [], pq => pq
  | true, a :: as, pq => cfracAux false as (hstep a pq)
  | false, a :: as, pq => cfracAux true as (vstep a pq)
def cfrac (l : List Int) : Int × Int := cfracAux true l (0, 1)

/-- The seven rational knots of the page: 3_1, 4_1, 5_1, 5_2, 6_1, 6_2, 6_3. -/
theorem rational_fractions :
    cfrac [3] = (3, 1) ∧ cfrac [1, 1, 2] = (5, 2) ∧ cfrac [5] = (5, 1) ∧ cfrac [1, 1, 3] = (7, 2) ∧
    cfrac [1, 1, 4] = (9, 2) ∧ cfrac [2, 1, 3] = (11, 3) ∧ cfrac [1, 1, 1, 1, 2] = (13, 5) := by decide

/-- Both steps keep the determinant p·s − q·r of any two fractions (they are in SL(2, ℤ)). -/
theorem hstep_det (a : Int) (P Q : Int × Int) :
    (hstep a P).1 * (hstep a Q).2 - (hstep a P).2 * (hstep a Q).1 = P.1 * Q.2 - P.2 * Q.1 := by
  unfold hstep; grind
theorem vstep_det (a : Int) (P Q : Int × Int) :
    (vstep a P).1 * (vstep a Q).2 - (vstep a P).2 * (vstep a Q).1 = P.1 * Q.2 - P.2 * Q.1 := by
  unfold vstep; grind

/-- Twist forward, back: the steps are invertible. -/
theorem hstep_inv (a : Int) (P : Int × Int) : hstep (-a) (hstep a P) = P := by
  unfold hstep; cases P; congr 1 <;> grind
theorem vstep_inv (a : Int) (P : Int × Int) : vstep (-a) (vstep a P) = P := by
  unfold vstep; cases P; congr 1 <;> grind

/-! ## 6. Laurent brackets: the tangle algebra reproduces the PD-code brackets of QAcross.Knots -/

open QKnots

/-- Laurent polynomial in A: (lowest power, coefficients). -/
abbrev LP := Int × Poly

def ladd (a b : LP) : LP :=
  let lo := min a.1 b.1
  (lo, padd (pzeros (a.1 - lo).toNat a.2) (pzeros (b.1 - lo).toNat b.2))
def lmul (a b : LP) : LP := (a.1 + b.1, pmul a.2 b.2)
def lA : LP := (1, [1])
def lAi : LP := (-1, [1])
def lone : LP := (0, [1])
def lzero : LP := (0, [])
def ld : LP := (-2, [-1, 0, 0, 0, -1])

def ltsum (S T : LP × LP) : LP × LP :=
  (lmul S.1 T.1, ladd (ladd (lmul S.1 T.2) (lmul S.2 T.1)) (lmul ld (lmul S.2 T.2)))
def lvsum (S T : LP × LP) : LP × LP :=
  (ladd (ladd (lmul S.1 T.2) (lmul S.2 T.1)) (lmul ld (lmul S.1 T.1)), lmul S.2 T.2)
def lcloseN (T : LP × LP) : LP := normLP (ladd (lmul ld T.1) T.2)
def lcloseD (T : LP × LP) : LP := normLP (ladd T.1 (lmul ld T.2))

/-- Build a tangle from a word: true = add a crossing to the right, false = stack one below. -/
def lbuild : List Bool → LP × LP
  | [] => (lone, lzero)
  | b :: bs => if b then ltsum (lbuild bs) (lA, lAi) else lvsum (lbuild bs) (lA, lAi)

/-- Three crossings side by side, closed on top: the right-handed trefoil, the mirror of the
    left trefoil of QAcross.Knots (bracket A⁻⁷ − A⁻³ − A⁵). -/
theorem tangle_trefoil :
    lcloseN (lbuild [true, true, true]) = mirrorLP (normLP (-9, bracketP 3 6 trefoilL)) := by
  decide +kernel

/-- h v h h closed on top (fraction 5/2) is the figure-eight of QAcross.Knots. -/
theorem tangle_figEight :
    lcloseN (lbuild [true, true, false, true]) = normLP (-12, bracketP 4 8 figEight) := by
  decide +kernel

/-- h h v v closed at the sides (fraction 2/5) is the figure-eight again. -/
theorem tangle_figEight_D :
    lcloseD (lbuild [false, false, true, true]) = normLP (-12, bracketP 4 8 figEight) := by
  decide +kernel

/-- Five crossings closed on top: the cinquefoil (mirror of QAcross.Knots' 5₁). -/
theorem tangle_cinquefoil :
    lcloseN (lbuild [true, true, true, true, true]) = mirrorLP (normLP (-15, bracketP 5 10 cinquefoilL)) := by
  decide +kernel

end KnotTheory
