module

public import FusedAlgebra.Cube
public import FusedAlgebra.Topics.QAcross.QDial

/-!
# Knots through the Kauffman bracket (the knot `q`)

The Kauffman bracket (1987) of a diagram with `n` crossings is a sum over the `2 ^ n` corners of
the `n`-cube (each crossing smoothed one of two ways), each corner weighted by
`A ^ (#A - #B) * d ^ (loops - 1)` with loop value `d = -A ^ 2 - A⁻¹ ^ 2`. This is the same cube
as the mult-add corner expansion; Khovanov homology is built on exactly this cube of states. The
Jones polynomial (1984) is the writhe-normalised bracket at `A = t ^ (-1 / 4)`. All of this is
standard; none of it is in mathlib.

The `q`-link: `d = -A⁻¹ ^ 2 * [2]` with `[2] = 1 + A ^ 4` the `q`-integer of `QDial` at
`q = A ^ 4` (`loop_value_q`). So the knot `q` and the counting `q` are one `q`.

## Representation

The knot invariants here are finite computations checked by the kernel (`decide +kernel`), so
they need a computable representation; mathlib's `Polynomial` is not computable. A polynomial is
a list of integer coefficients, lowest power first (`QKnots.Poly`), and `QKnots.toPolynomial`
with the lemmas `toPolynomial_padd`, `toPolynomial_pscale`, `toPolynomial_pmul`,
`toPolynomial_ppow`, `toPolynomial_pzeros`, `toPolynomial_trim` says that the list operations are
the operations of `ℤ[X]`.

A diagram is a PD code (Knot Atlas convention `X[a, b, c, d]`: `a` the incoming under-strand,
then counterclockwise). The bracket is returned as `P(A)` with `⟨D⟩ = A ^ (-3 * n) * P(A)`. The
normalised invariant `f = (-A ^ 3) ^ (-w) * ⟨D⟩` (`w` the writhe) is returned as a pair (lowest
power, coefficients).

## Main results

* `trefoilL_f`, `trefoilR_f`, `figEight_f`, `hopf_f`, `cinquefoilL_f`: the invariants.
* `trefoil_chiral`, `figEight_amphichiral`: the trefoil differs from its mirror image, the
  figure-eight does not.
* `jones_at_one`, `jones_at_neg_one`: the dial points `t = 1` (components) and `t = -1`
  (determinant).
* `state_layers`: the states with `k` B-smoothings are counted by `n.choose k`, for every `n`
  and `k` (the original checked `n, k ≤ 5`).
* `reidemeister_one`, `reidemeister_two`, `writhe_cancels`: the local rules in any commutative
  ring.
-/

@[expose] public section

namespace QKnots

open MultAdd Polynomial

/-! ### Polynomials as coefficient lists -/

/-- A polynomial with integer coefficients as its list of coefficients, lowest power first.
This is the computable stand-in for `ℤ[X]`; see `toPolynomial`. -/
abbrev Poly := List ℤ

/-- Sum of two coefficient lists. -/
def padd : Poly → Poly → Poly
  | [], q => q
  | a :: p, [] => a :: p
  | a :: p, b :: q => (a + b) :: padd p q

/-- Multiply a coefficient list by a constant. -/
def pscale (c : ℤ) (p : Poly) : Poly := p.map (c * ·)

/-- Product of two coefficient lists. -/
def pmul : Poly → Poly → Poly
  | [], _ => []
  | a :: p, q => padd (pscale a q) (0 :: pmul p q)

/-- Power of a coefficient list. -/
def ppow (p : Poly) : ℕ → Poly
  | 0 => [1]
  | n + 1 => pmul p (ppow p n)

/-- Multiply by `X ^ k`: put `k` zeros in front. -/
def pzeros (k : ℕ) (p : Poly) : Poly := List.replicate k 0 ++ p

/-- Drop the zero coefficients at the top. -/
def trim (p : Poly) : Poly := p.rdropWhile (· == 0)

/-- The polynomial in `ℤ[X]` with the given list of coefficients. -/
noncomputable def toPolynomial : Poly → ℤ[X]
  | [] => 0
  | a :: p => C a + X * toPolynomial p

@[simp]
theorem toPolynomial_nil : toPolynomial [] = 0 := rfl

@[simp]
theorem toPolynomial_cons (a : ℤ) (p : Poly) :
    toPolynomial (a :: p) = C a + X * toPolynomial p := rfl

theorem toPolynomial_padd : ∀ p q : Poly, toPolynomial (padd p q) = toPolynomial p + toPolynomial q
  | [], q => by simp [padd]
  | a :: p, [] => by simp [padd]
  | a :: p, b :: q => by simp only [padd, toPolynomial_cons, toPolynomial_padd p q, C_add]; ring

theorem toPolynomial_pscale (c : ℤ) (p : Poly) :
    toPolynomial (pscale c p) = C c * toPolynomial p := by
  induction p with
  | nil => simp [pscale]
  | cons a p ih =>
    rw [pscale, List.map_cons, ← pscale, toPolynomial_cons, toPolynomial_cons, ih, C_mul]; ring

theorem toPolynomial_pmul (p q : Poly) :
    toPolynomial (pmul p q) = toPolynomial p * toPolynomial q := by
  induction p with
  | nil => simp [pmul]
  | cons a p ih =>
    simp only [pmul, toPolynomial_padd, toPolynomial_pscale, toPolynomial_cons, ih, C_0]
    ring

theorem toPolynomial_ppow (p : Poly) (n : ℕ) : toPolynomial (ppow p n) = toPolynomial p ^ n := by
  induction n with
  | zero => simp [ppow]
  | succ n ih => rw [ppow, toPolynomial_pmul, ih, pow_succ']

theorem toPolynomial_pzeros (k : ℕ) (p : Poly) :
    toPolynomial (pzeros k p) = X ^ k * toPolynomial p := by
  induction k with
  | zero => simp [pzeros]
  | succ k ih =>
    rw [pzeros, List.replicate_succ, List.cons_append, ← pzeros, toPolynomial_cons, ih]
    simp only [C_0]; ring

theorem toPolynomial_append_zero (p : Poly) : toPolynomial (p ++ [0]) = toPolynomial p := by
  induction p with
  | nil => simp
  | cons a p ih => rw [List.cons_append, toPolynomial_cons, ih, toPolynomial_cons]

/-- Trimming the zeros at the top does not change the polynomial. -/
theorem toPolynomial_trim (p : Poly) : toPolynomial (trim p) = toPolynomial p := by
  induction p using List.reverseRecOn with
  | nil => rfl
  | append_singleton p a ih =>
    by_cases h : a = 0
    · rw [trim, List.rdropWhile_concat_pos _ _ _ (by simpa using h), ← trim, ih, h,
        toPolynomial_append_zero]
    · rw [trim, List.rdropWhile_concat_neg _ _ _ (by simpa using h)]

/-! ### Diagrams, smoothings, loops -/

/-- A crossing `X[a, b, c, d]` of a PD code: four arc labels. -/
abbrev Crossing := ℕ × ℕ × ℕ × ℕ

/-- The A-smoothing (`false`) joins `a`-`b` and `c`-`d`; the B-smoothing (`true`) joins `a`-`d`
and `b`-`c`. -/
def smooth : Crossing → Bool → List (ℕ × ℕ)
  | (a, b, c, d), false => [(a, b), (c, d)]
  | (a, b, c, d), true => [(a, d), (b, c)]

/-- Merge the classes of `p.1` and `p.2` (`comp` sends each arc label to the name of its
class). -/
def relabel (comp : List ℕ) (p : ℕ × ℕ) : List ℕ :=
  if comp.getD p.1 0 == comp.getD p.2 0 then comp
  else comp.map fun c ↦ if c == comp.getD p.2 0 then comp.getD p.1 0 else c

/-- The number of loops a smoothing makes: arc labels `1, …, m`, joined by the pairs. -/
def countLoops (m : ℕ) (pairs : List (ℕ × ℕ)) : ℕ :=
  (((List.range (m + 1)).zip (pairs.foldl relabel (List.range (m + 1)))).filter
    fun ic ↦ decide (1 ≤ ic.1) && ic.2 == ic.1).length

/-- All `2 ^ n` states: the corners of the `n`-cube. -/
def allStates : ℕ → List (List Bool)
  | 0 => [[]]
  | n + 1 => (allStates n).map (false :: ·) ++ (allStates n).map (true :: ·)

/-- The arcs joined by a state: one smoothing for each crossing. -/
def statePairs : List Crossing → List Bool → List (ℕ × ℕ)
  | x :: xs, b :: bs => smooth x b ++ statePairs xs bs
  | _, _ => []

/-- The number of A-smoothings of a state. -/
def numA (s : List Bool) : ℕ := s.count false

/-- The number of B-smoothings of a state. -/
def numB (s : List Bool) : ℕ := s.count true

/-- The term of a corner times `A ^ (3 * n)`: `A ^ (#A - #B) * d ^ (L - 1)` with
`d = -A⁻¹ ^ 2 * (1 + A ^ 4)`, shifted so every power is nonnegative. Its power of `A` is
`3 * n + #A - #B - 2 * (L - 1) = 2 * n + 2 * #A + 2 - 2 * L`. -/
def stateTerm (n m : ℕ) (pd : List Crossing) (s : List Bool) : Poly :=
  pzeros (2 * n + 2 * numA s + 2 - 2 * countLoops m (statePairs pd s))
    (pscale ((-1) ^ (countLoops m (statePairs pd s) - 1))
      (ppow [1, 0, 0, 0, 1] (countLoops m (statePairs pd s) - 1)))

/-- The Kauffman bracket times `A ^ (3 * n)`. -/
def bracketP (n m : ℕ) (pd : List Crossing) : Poly :=
  trim ((allStates n).foldl (fun acc s ↦ padd acc (stateTerm n m pd s)) [])

/-- Strip the zero coefficients at the bottom, keeping track of the lowest power. -/
def dropLead : ℤ → Poly → ℤ × Poly
  | _, [] => (0, [])
  | lo, c :: cs => if c == 0 then dropLead (lo + 1) cs else (lo, c :: cs)

/-- (lowest power, coefficients), with zeros stripped at both ends. -/
def normLP (lp : ℤ × Poly) : ℤ × Poly := dropLead lp.1 (trim lp.2)

/-- `f = (-A ^ 3) ^ (-w) * ⟨D⟩`: the bracket made invariant under all three Reidemeister
moves. -/
def kauffmanF (n m : ℕ) (pd : List Crossing) (w : ℤ) : ℤ × Poly :=
  normLP (-3 * w - 3 * (n : ℤ), pscale (if w % 2 == 0 then 1 else -1) (bracketP n m pd))

/-- The mirror image: `A ↦ A⁻¹`. -/
def mirrorLP (lp : ℤ × Poly) : ℤ × Poly :=
  (-(lp.1 + (lp.2.length : ℤ) - 1), lp.2.reverse)

/-- The value at `A = 1` (so `t = 1`). -/
def evalOne (lp : ℤ × Poly) : ℤ := lp.2.sum

/-- The value at `A ^ 4 = -1` of the coefficients starting at power `e`, when all powers are
multiples of `4`. -/
def evalTnegAux : ℤ → Poly → ℤ
  | _, [] => 0
  | e, c :: cs => (if e % 8 == 0 then c else if e % 8 == 4 then -c else 0) + evalTnegAux (e + 1) cs

/-- The value at `t = -1` (`A ^ 4 = -1`), for invariants whose powers are all multiples of `4`
(knots). -/
def evalTneg (lp : ℤ × Poly) : ℤ := evalTnegAux lp.1 lp.2

/-! ### The diagrams (Knot Atlas PD codes) -/

/-- The unknot drawn with one twist (a Reidemeister I kink), writhe `+1`. -/
def kinkUnknot : List Crossing := [(1, 1, 2, 2)]

/-- The left-handed trefoil `3₁`, writhe `-3`. -/
def trefoilL : List Crossing := [(1, 4, 2, 5), (3, 6, 4, 1), (5, 2, 6, 3)]

/-- The right-handed trefoil: the mirror (the labels of each crossing rotated by one), writhe
`+3`. -/
def trefoilR : List Crossing := [(4, 2, 5, 1), (6, 4, 1, 3), (2, 6, 3, 5)]

/-- The figure-eight knot `4₁`, writhe `0`. -/
def figEight : List Crossing := [(4, 2, 5, 1), (8, 6, 1, 5), (6, 3, 7, 4), (2, 7, 3, 8)]

/-- The Hopf link `L2a1`, oriented with writhe `+2`. -/
def hopf : List Crossing := [(4, 1, 3, 2), (2, 3, 1, 4)]

/-- The cinquefoil `5₁` (the torus knot `T(2, 5)`), writhe `-5`. -/
def cinquefoilL : List Crossing :=
  [(1, 6, 2, 7), (3, 8, 4, 9), (5, 10, 6, 1), (7, 2, 8, 3), (9, 4, 10, 5)]

/-! ### Brackets and invariants (finite computations, checked by the kernel) -/

/-- `⟨kink⟩ = -A ^ 3` (stored times `A ^ 3`, as `-A ^ 6`). -/
theorem kink_bracket : bracketP 1 2 kinkUnknot = [0, 0, 0, 0, 0, 0, -1] := by decide +kernel

/-- The writhe correction removes the kink: the invariant of the unknot is `1`. -/
theorem unknot_f : kauffmanF 1 2 kinkUnknot 1 = (0, [1]) := by decide +kernel

/-- `⟨left trefoil⟩ = A ^ 7 - A ^ 3 - A ^ (-5)` (stored times `A ^ 9`). -/
theorem trefoilL_bracket :
    bracketP 3 6 trefoilL = [0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 1] := by
  decide +kernel

/-- Left trefoil: `f = A ^ 4 + A ^ 12 - A ^ 16`, so `V(t) = t⁻¹ + t⁻³ - t⁻⁴`. -/
theorem trefoilL_f :
    kauffmanF 3 6 trefoilL (-3) = (4, [1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -1]) := by
  decide +kernel

/-- Right trefoil: `f = -A ^ (-16) + A ^ (-12) + A ^ (-4)`, so `V(t) = t + t³ - t⁴`. -/
theorem trefoilR_f :
    kauffmanF 3 6 trefoilR 3 = (-16, [-1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1]) := by
  decide +kernel

/-- The two trefoils are mirror images of each other ... -/
theorem trefoil_mirror : mirrorLP (kauffmanF 3 6 trefoilL (-3)) = kauffmanF 3 6 trefoilR 3 := by
  decide +kernel

/-- ... and the invariant tells them apart: the trefoil is chiral. -/
theorem trefoil_chiral : kauffmanF 3 6 trefoilL (-3) ≠ kauffmanF 3 6 trefoilR 3 := by
  decide +kernel

/-- Figure-eight: `f = A ^ (-8) - A ^ (-4) + 1 - A ^ 4 + A ^ 8`, so
`V(t) = t² - t + 1 - t⁻¹ + t⁻²`. -/
theorem figEight_f :
    kauffmanF 4 8 figEight 0 = (-8, [1, 0, 0, 0, -1, 0, 0, 0, 1, 0, 0, 0, -1, 0, 0, 0, 1]) := by
  decide +kernel

/-- The figure-eight equals its own mirror image (it is amphichiral), unlike the trefoil. -/
theorem figEight_amphichiral :
    mirrorLP (kauffmanF 4 8 figEight 0) = kauffmanF 4 8 figEight 0 := by
  decide +kernel

/-- Hopf link: `f = -A ^ (-10) - A ^ (-2)`, so `V(t) = -t ^ (1 / 2) - t ^ (5 / 2)`. -/
theorem hopf_f : kauffmanF 2 4 hopf 2 = (-10, [-1, 0, 0, 0, 0, 0, 0, 0, -1]) := by
  decide +kernel

/-- Cinquefoil: `V(t) = t⁻² + t⁻⁴ - t⁻⁵ + t⁻⁶ - t⁻⁷`. -/
theorem cinquefoilL_f :
    kauffmanF 5 10 cinquefoilL (-5) =
      (8, [1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -1, 0, 0, 0, 1, 0, 0, 0, -1]) := by
  decide +kernel

/-! ### The dial points of the knot `q` -/

/-- At `t = 1` (the counting point, weight `0`) the Jones polynomial only counts components:
`1` for every knot, `-2` for the two-component Hopf link. -/
theorem jones_at_one :
    evalOne (kauffmanF 3 6 trefoilL (-3)) = 1 ∧ evalOne (kauffmanF 3 6 trefoilR 3) = 1 ∧
    evalOne (kauffmanF 4 8 figEight 0) = 1 ∧ evalOne (kauffmanF 5 10 cinquefoilL (-5)) = 1 ∧
    evalOne (kauffmanF 2 4 hopf 2) = -2 := by
  decide +kernel

/-- At `t = -1` (the parity point, weight `-2`) the Jones polynomial gives the determinant of the
knot up to sign: `3` for the trefoil, `5` for the figure-eight and for the cinquefoil. -/
theorem jones_at_neg_one :
    evalTneg (kauffmanF 3 6 trefoilL (-3)) = -3 ∧ evalTneg (kauffmanF 3 6 trefoilR 3) = -3 ∧
    evalTneg (kauffmanF 4 8 figEight 0) = 5 ∧ evalTneg (kauffmanF 5 10 cinquefoilL (-5)) = 5 := by
  decide +kernel

/-! ### The state cube is the mult-add cube -/

/-- The eight states of the trefoil as (number of B-smoothings, loops). The B-count runs over
the layers of the 3-cube (`1, 3, 3, 1`). -/
theorem trefoilL_cube :
    (allStates 3).map (fun s ↦ (numB s, countLoops 6 (statePairs trefoilL s))) =
      [(0, 3), (1, 2), (1, 2), (2, 1), (1, 2), (2, 1), (2, 1), (3, 2)] := by
  decide +kernel

/-- The cube has `2 ^ n` corners. -/
theorem length_allStates (n : ℕ) : (allStates n).length = 2 ^ n := by
  induction n with
  | zero => rfl
  | succ n ih => simp [allStates, ih, pow_succ, mul_two]

/-- The states with `k` B-smoothings are counted by the layers of the cube (Pascal's triangle),
for every number of crossings. -/
theorem state_layers : ∀ n k : ℕ, (allStates n).countP (fun s ↦ numB s = k) = n.choose k
  | 0, 0 => by simp [allStates, numB]
  | 0, k + 1 => by simp [allStates, numB]
  | n + 1, 0 => by
    have h := state_layers n 0
    simp only [numB, Nat.choose_zero_right] at h
    simp only [allStates, numB, List.countP_append, List.countP_map, Function.comp_def,
      List.count_cons_self, List.count_cons_of_ne Bool.false_ne_true,
      Nat.succ_ne_zero, decide_false, List.countP_false, Nat.choose_zero_right]
    exact h
  | n + 1, k + 1 => by
    have h₁ := state_layers n k
    have h₂ := state_layers n (k + 1)
    simp only [numB] at h₁ h₂
    simp only [allStates, numB, List.countP_append, List.countP_map, Function.comp_def,
      List.count_cons_self, List.count_cons_of_ne Bool.false_ne_true,
      Nat.add_right_cancel_iff, Nat.choose_succ_succ]
    rw [add_comm]
    exact congrArg₂ (· + ·) h₁ h₂

/-! ### The local rules, in any commutative ring (`A * Ai = 1`, so `Ai = A⁻¹`) -/

section Local

variable {R : Type*} [CommRing R]

/-- The loop value is minus a `q`-integer: `d = -A⁻¹ ^ 2 * [2]` with
`[2] = 1 + A ^ 4 = QDial.qnum (A ^ 4) 2`. -/
theorem loop_value_q (A Ai : R) (h : A * Ai = 1) :
    -A ^ 2 - Ai ^ 2 = -Ai ^ 2 * QDial.qnum (A ^ 4) 2 := by
  simp only [QDial.qnum, Finset.sum_range_succ, Finset.sum_range_zero]
  linear_combination (A ^ 2 * (A * Ai + 1)) * h

/-- Reidemeister I: a kink multiplies the bracket by `A * d + A⁻¹ = -A ^ 3`. -/
theorem reidemeister_one (A Ai : R) (h : A * Ai = 1) :
    A * (-A ^ 2 - Ai ^ 2) + Ai = -A ^ 3 := by
  linear_combination (-Ai) * h

/-- Reidemeister II: in the bracket of two crossings, the cup-cap term has coefficient
`A ^ 2 + A * A⁻¹ * d + A⁻¹ ^ 2 = 0` (and the parallel term has coefficient `A * A⁻¹ = 1`). -/
theorem reidemeister_two (A Ai : R) (h : A * Ai = 1) :
    A ^ 2 + A * Ai * (-A ^ 2 - Ai ^ 2) + Ai ^ 2 = 0 := by
  linear_combination (-(A ^ 2 + Ai ^ 2)) * h

/-- The writhe factor undoes the kink: `(-A ^ 3) * (-A⁻¹ ^ 3) = 1`. -/
theorem writhe_cancels (A Ai : R) (h : A * Ai = 1) : -A ^ 3 * -Ai ^ 3 = 1 := by
  linear_combination ((A * Ai) ^ 2 + A * Ai + 1) * h

end Local

end QKnots
