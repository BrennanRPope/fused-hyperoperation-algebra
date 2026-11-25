/-
  Fused Hyperoperation Algebra: knots through the Kauffman bracket (the knot q)
  ============================================================================
  Core Lean 4 only (no mathlib), no `sorry`.

  LIBRARY NOTE: copied from q-across/lean/Knots.lean (Nov 25) and compiled here for the
  first time: no sorry, standard axioms only; one fix: the section's `[CommRing R]` became `[Lean.Grind.CommRing R]`
  (core Lean's ring class; the bare name doesn't resolve without `open Lean.Grind`). No statement changed.

  What is standard: the Kauffman bracket (1987) and Jones polynomial (1984). The bracket of a
  diagram with n crossings is a sum over the 2ⁿ corners of the n-cube (each crossing smoothed
  one of two ways), each corner weighted by A^(#A − #B)·d^(loops − 1) with loop value
  d = −A² − A⁻². This is the same cube as the circle-operation corner (called mult-add corner before 1.0) expansion; Khovanov homology
  (1999) is built on exactly this cube of states.

  The q-link: d = −A⁻²·[2] with [2] = 1 + A⁴ the q-integer of QDial at q = A⁴, and the Jones
  variable is t = A⁻⁴. So the knot q and the counting q are one q.

  Representation: a diagram is a PD code (Knot Atlas convention X[a,b,c,d], a = incoming under
  strand, then counterclockwise). A polynomial is a list of Int coefficients, index = power.
  The bracket is returned as P(A) with ⟨D⟩ = A^(−3n)·P(A). The normalized invariant
  f = (−A³)^(−w)·⟨D⟩ (w = writhe, read off the oriented diagram) is returned as
  (lowest power, coefficients). The Jones polynomial is f with A = t^(−1/4).
-/
import FusedAlgebra.Base
import FusedAlgebra.Cube
import FusedAlgebra.Topics.QAcross.QDial

namespace QKnots
open CircleOp

/-! ## Polynomials as coefficient lists -/

abbrev Poly := List Int

def padd : Poly → Poly → Poly
  | [], q => q
  | a :: p, [] => a :: p
  | a :: p, b :: q => (a + b) :: padd p q

def pscale (c : Int) (p : Poly) : Poly := p.map (fun x => c * x)

def pmul : Poly → Poly → Poly
  | [], _ => []
  | a :: p, q => padd (pscale a q) (0 :: pmul p q)

def ppow (p : Poly) : Nat → Poly
  | 0 => [1]
  | n + 1 => pmul p (ppow p n)

/-- Multiply by A^k. -/
def pzeros : Nat → Poly → Poly
  | 0, p => p
  | k + 1, p => 0 :: pzeros k p

/-- Drop zero coefficients at the top. -/
def trim (p : Poly) : Poly := (p.reverse.dropWhile (fun x => x == 0)).reverse

def negOnePow : Nat → Int
  | 0 => 1
  | n + 1 => -negOnePow n

/-! ## Diagrams, smoothings, loops -/

abbrev Crossing := Nat × Nat × Nat × Nat

/-- The A-smoothing (false) joins a–b and c–d; the B-smoothing (true) joins a–d and b–c. -/
def smooth : Crossing → Bool → List (Nat × Nat)
  | (a, b, c, d), false => [(a, b), (c, d)]
  | (a, b, c, d), true => [(a, d), (b, c)]

/-- Merge the classes of p.1 and p.2 (comp sends each arc label to its class's name). -/
def relabel (comp : List Nat) (p : Nat × Nat) : List Nat :=
  if comp.getD p.1 0 == comp.getD p.2 0 then comp
  else comp.map (fun c => if c == comp.getD p.2 0 then comp.getD p.1 0 else c)

/-- Number of loops a smoothing makes: arc labels 1..m, joined by the pairs. -/
def countLoops (m : Nat) (pairs : List (Nat × Nat)) : Nat :=
  (((List.range (m + 1)).zip (pairs.foldl relabel (List.range (m + 1)))).filter
    (fun ic => decide (1 ≤ ic.1) && ic.2 == ic.1)).length

/-- All 2ⁿ states: the corners of the n-cube. -/
def allStates : Nat → List (List Bool)
  | 0 => [[]]
  | n + 1 => (allStates n).map (fun s => false :: s) ++ (allStates n).map (fun s => true :: s)

def statePairs : List Crossing → List Bool → List (Nat × Nat)
  | x :: xs, b :: bs => smooth x b ++ statePairs xs bs
  | _, _ => []

def numA (s : List Bool) : Nat := (s.filter (fun b => !b)).length
def numB (s : List Bool) : Nat := (s.filter (fun b => b)).length

/-- A corner's term times A^(3n): A^(#A − #B)·d^(L−1) with d = −A⁻²(1 + A⁴), shifted so every
    power is ≥ 0. Its power of A is 3n + #A − #B − 2(L − 1) = 2n + 2·#A + 2 − 2L. -/
def stateTerm (n m : Nat) (pd : List Crossing) (s : List Bool) : Poly :=
  pzeros (2 * n + 2 * numA s + 2 - 2 * countLoops m (statePairs pd s))
    (pscale (negOnePow (countLoops m (statePairs pd s) - 1))
      (ppow [1, 0, 0, 0, 1] (countLoops m (statePairs pd s) - 1)))

/-- The Kauffman bracket times A^(3n). -/
def bracketP (n m : Nat) (pd : List Crossing) : Poly :=
  trim ((allStates n).foldl (fun acc s => padd acc (stateTerm n m pd s)) [])

/-- (lowest power, coefficients), with zeros stripped at both ends. -/
def dropLead : Int → Poly → Int × Poly
  | _, [] => (0, [])
  | lo, c :: cs => if c == 0 then dropLead (lo + 1) cs else (lo, c :: cs)

def normLP (lp : Int × Poly) : Int × Poly := dropLead lp.1 (trim lp.2)

/-- f = (−A³)^(−w)·⟨D⟩: the bracket made invariant under all three Reidemeister moves. -/
def kauffmanF (n m : Nat) (pd : List Crossing) (w : Int) : Int × Poly :=
  normLP (-3 * w - 3 * (n : Int), pscale (if w % 2 == 0 then 1 else -1) (bracketP n m pd))

/-- The mirror image: A ↦ A⁻¹. -/
def mirrorLP (lp : Int × Poly) : Int × Poly :=
  (-(lp.1 + (lp.2.length : Int) - 1), lp.2.reverse)

/-- Value at A = 1 (so t = 1). -/
def evalOne (lp : Int × Poly) : Int := lp.2.foldl (fun x y => x + y) 0

/-- Value at t = −1 (A⁴ = −1), for invariants whose powers are all multiples of 4 (knots). -/
def evalTnegAux : Int → Poly → Int
  | _, [] => 0
  | e, c :: cs => (if e % 8 == 0 then c else if e % 8 == 4 then -c else 0) + evalTnegAux (e + 1) cs

def evalTneg (lp : Int × Poly) : Int := evalTnegAux lp.1 lp.2

/-! ## The diagrams (Knot Atlas PD codes) -/

/-- The unknot drawn with one twist (Reidemeister I kink), writhe +1. -/
def kinkUnknot : List Crossing := [(1, 1, 2, 2)]
/-- Left-handed trefoil 3₁, writhe −3. -/
def trefoilL : List Crossing := [(1, 4, 2, 5), (3, 6, 4, 1), (5, 2, 6, 3)]
/-- Right-handed trefoil: the mirror (each crossing's labels rotated by one), writhe +3. -/
def trefoilR : List Crossing := [(4, 2, 5, 1), (6, 4, 1, 3), (2, 6, 3, 5)]
/-- Figure-eight knot 4₁, writhe 0. -/
def figEight : List Crossing := [(4, 2, 5, 1), (8, 6, 1, 5), (6, 3, 7, 4), (2, 7, 3, 8)]
/-- Hopf link L2a1, oriented with writhe +2. -/
def hopf : List Crossing := [(4, 1, 3, 2), (2, 3, 1, 4)]
/-- Cinquefoil 5₁ (torus knot T(2,5)), writhe −5. -/
def cinquefoilL : List Crossing :=
  [(1, 6, 2, 7), (3, 8, 4, 9), (5, 10, 6, 1), (7, 2, 8, 3), (9, 4, 10, 5)]

/-! ## Brackets and invariants (finite computations, checked by the kernel) -/

/-- ⟨kink⟩ = −A³ (stored times A³ as −A⁶). -/
theorem kink_bracket : bracketP 1 2 kinkUnknot = [0, 0, 0, 0, 0, 0, -1] := by decide +kernel

/-- The writhe correction removes the kink: the unknot's invariant is 1. -/
theorem unknot_f : kauffmanF 1 2 kinkUnknot 1 = (0, [1]) := by decide +kernel

/-- ⟨left trefoil⟩ = A⁷ − A³ − A⁻⁵ (stored times A⁹). -/
theorem trefoilL_bracket :
    bracketP 3 6 trefoilL = [0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 1] := by
  decide +kernel

/-- Left trefoil: f = A⁴ + A¹² − A¹⁶, so V(t) = t⁻¹ + t⁻³ − t⁻⁴. -/
theorem trefoilL_f :
    kauffmanF 3 6 trefoilL (-3) = (4, [1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -1]) := by
  decide +kernel

/-- Right trefoil: f = −A⁻¹⁶ + A⁻¹² + A⁻⁴, so V(t) = t + t³ − t⁴. -/
theorem trefoilR_f :
    kauffmanF 3 6 trefoilR 3 = (-16, [-1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1]) := by
  decide +kernel

/-- The two trefoils are mirror images of each other ... -/
theorem trefoil_mirror : mirrorLP (kauffmanF 3 6 trefoilL (-3)) = kauffmanF 3 6 trefoilR 3 := by
  decide +kernel

/-- ... and the invariant tells them apart: the trefoil is chiral. -/
theorem trefoil_chiral : kauffmanF 3 6 trefoilL (-3) ≠ kauffmanF 3 6 trefoilR 3 := by
  decide +kernel

/-- Figure-eight: f = A⁻⁸ − A⁻⁴ + 1 − A⁴ + A⁸, so V(t) = t² − t + 1 − t⁻¹ + t⁻². -/
theorem figEight_f :
    kauffmanF 4 8 figEight 0 = (-8, [1, 0, 0, 0, -1, 0, 0, 0, 1, 0, 0, 0, -1, 0, 0, 0, 1]) := by
  decide +kernel

/-- The figure-eight equals its own mirror image (it is amphichiral), unlike the trefoil. -/
theorem figEight_amphichiral : mirrorLP (kauffmanF 4 8 figEight 0) = kauffmanF 4 8 figEight 0 := by
  decide +kernel

/-- Hopf link: f = −A⁻¹⁰ − A⁻², so V(t) = −t^(1/2) − t^(5/2). -/
theorem hopf_f : kauffmanF 2 4 hopf 2 = (-10, [-1, 0, 0, 0, 0, 0, 0, 0, -1]) := by
  decide +kernel

/-- Cinquefoil: V(t) = t⁻² + t⁻⁴ − t⁻⁵ + t⁻⁶ − t⁻⁷. -/
theorem cinquefoilL_f :
    kauffmanF 5 10 cinquefoilL (-5) =
      (8, [1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -1, 0, 0, 0, 1, 0, 0, 0, -1]) := by
  decide +kernel

/-! ## The dial points of the knot q -/

/-- At t = 1 (the counting point, k = 0) the Jones polynomial only counts components:
    1 for every knot, −2 for the two-component Hopf link. -/
theorem jones_at_one :
    evalOne (kauffmanF 3 6 trefoilL (-3)) = 1 ∧ evalOne (kauffmanF 3 6 trefoilR 3) = 1 ∧
    evalOne (kauffmanF 4 8 figEight 0) = 1 ∧ evalOne (kauffmanF 5 10 cinquefoilL (-5)) = 1 ∧
    evalOne (kauffmanF 2 4 hopf 2) = -2 := by
  decide +kernel

/-- At t = −1 (the parity point, k = −2) the Jones polynomial gives ± the knot's determinant:
    3 for the trefoil (its three 3-colorings beyond the trivial), 5 for the figure-eight and
    for the cinquefoil. -/
theorem jones_at_neg_one :
    evalTneg (kauffmanF 3 6 trefoilL (-3)) = -3 ∧ evalTneg (kauffmanF 3 6 trefoilR 3) = -3 ∧
    evalTneg (kauffmanF 4 8 figEight 0) = 5 ∧ evalTneg (kauffmanF 5 10 cinquefoilL (-5)) = 5 := by
  decide +kernel

/-! ## The state cube is the circle-operation cube -/

/-- The trefoil's eight states as (number of B-smoothings, loops). The B-count runs over the
    layers of the 3-cube (1, 3, 3, 1). -/
theorem trefoilL_cube :
    (allStates 3).map (fun s => (numB s, countLoops 6 (statePairs trefoilL s))) =
      [(0, 3), (1, 2), (1, 2), (2, 1), (1, 2), (2, 1), (2, 1), (3, 2)] := by
  decide +kernel

/-- The states with k B-smoothings are counted by the cube's layers (Pascal), up to 5 crossings. -/
theorem state_layers :
    ∀ n : Fin 6, ∀ k : Fin 6,
      ((allStates n.val).filter (fun s => numB s == k.val)).length = CircleOp.Cube.layer n.val k.val := by
  decide +kernel

/-! ## The local rules, for any commutative ring (A·Ai = 1, so Ai = A⁻¹) -/

section Local
variable {R : Type} [Lean.Grind.CommRing R]

/-- The loop value is minus a q-integer: d = −A⁻²·[2] with [2] = 1 + A⁴ = QDial.qnum (A⁴) 2. -/
theorem loop_value_q (A Ai : R) (h : A * Ai = 1) :
    -(A * A) - Ai * Ai = -(Ai * Ai) * QDial.qnum (A * A * A * A) 2 := by
  show -(A * A) - Ai * Ai = -(Ai * Ai) * (1 + (A * A * A * A) * (1 + (A * A * A * A) * 0))
  grind

/-- Reidemeister I: a kink multiplies the bracket by A·d + A⁻¹ = −A³. -/
theorem reidemeister_one (A Ai : R) (h : A * Ai = 1) :
    A * (-(A * A) - Ai * Ai) + Ai = -(A * A * A) := by grind

/-- Reidemeister II: in the bracket of two crossings, the cup-cap term has coefficient
    A² + A·A⁻¹·d + A⁻² = 0, and the parallel term has coefficient A·A⁻¹ = 1. -/
theorem reidemeister_two (A Ai : R) (h : A * Ai = 1) :
    A * A + A * Ai * (-(A * A) - Ai * Ai) + Ai * Ai = 0 := by grind

/-- The writhe factor undoes the kink: (−A³)·(−A⁻³) = 1. -/
theorem writhe_cancels (A Ai : R) (h : A * Ai = 1) :
    (-(A * A * A)) * (-(Ai * Ai * Ai)) = 1 := by grind

end Local

end QKnots
