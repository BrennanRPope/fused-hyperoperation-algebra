/-
  Tesseract coordinates for Fused Hyperoperation Algebra
  ================================================

  STATUS: checked. Compiled against mathlib v4.37.0-rc1 on Lean 4.37.0-rc1
  on Brennan's Mac Studio (Nov 25, 2026) with zero errors and no `sorry`.
  An earlier compile failed only because `⊞`, the mult-add symbol used here
  before, is mathlib's global notation for biproducts; mult-add was then
  written with the eight-spoked asterisk (U+2733), which mathlib does not
  use. In the 1.0 rename it became the circle operation `◯` (typed
  \bigcirc, free in mathlib, scoped to `CircleOp`); this file has not been
  rechecked against mathlib since that rename. To recheck, put the file in a
  Lean project that depends on mathlib and run `lake env lean
  TesseractCoords.lean`.

  Notation: Lean already uses `⊕` for sum types and mathlib uses `⊞` for
  biproducts, so the operation is written `◯` here, the circle operation (called mult-add before 1.0) of
  ring theory. a ◯ b = a + b + a*b.

  The picture. Cut a box with sides 1+a, 1+b, 1+c, 1+d at 1 on every axis.
  It falls into 16 pieces, one per corner of the tesseract. A corner is a
  set of switches S (which axes are "on"); its piece has measure
  ∏_{i ∈ S} a i and its grade (dimension) is |S|.
-/
import Mathlib

open Finset BigOperators

namespace CircleOp

variable {R : Type*} [CommRing R]

/-! ## The circle operation on numbers -/

/-- The circle operation: one operation that adds and multiplies at once. -/
def circ (a b : R) : R := a + b + a * b

@[inherit_doc] scoped infixl:65 " ◯ " => circ

/-- The geometric heart: the circle operation is a box with sides 1 + a and 1 + b,
minus the unit square. -/
theorem one_add_circ (a b : R) : 1 + (a ◯ b) = (1 + a) * (1 + b) := by
  unfold circ; ring

theorem circ_comm (a b : R) : a ◯ b = b ◯ a := by
  unfold circ; ring

theorem circ_assoc (a b c : R) : (a ◯ b) ◯ c = a ◯ (b ◯ c) := by
  unfold circ; ring

theorem circ_zero (a : R) : a ◯ 0 = a := by
  unfold circ; ring

/-- Four numbers: the tesseract box with 16 pieces. -/
theorem one_add_circ₄ (a b c d : R) :
    1 + (a ◯ b ◯ c ◯ d) = (1 + a) * (1 + b) * (1 + c) * (1 + d) := by
  simp only [circ]; ring

/-- The circle operation of a whole list (0 is the identity). -/
def circList : List R → R
  | [] => 0
  | x :: xs => x ◯ circList xs

theorem one_add_circList :
    ∀ l : List R, 1 + circList l = (l.map (fun x => 1 + x)).prod
  | [] => by simp [circList]
  | x :: xs => by
      simp only [circList, one_add_circ, List.map_cons, List.prod_cons,
        one_add_circList xs]

/-! ## Tesseract coordinates -/

/-- A corner of the n-cube: the set of switches that are on. -/
abbrev Corner (n : ℕ) := Finset (Fin n)

/-- Tesseract coordinates: one weight for every corner. -/
abbrev Coord (n : ℕ) (R : Type*) := Corner n → R

/-- The grade of a corner is how many switches are on: its dimension
(0 = point, 1 = line, 2 = area, 3 = volume, 4 = hypervolume). -/
def grade {n : ℕ} (S : Corner n) : ℕ := S.card

/-- Corner coordinates of an ordinary point a: corner S holds ∏_{i ∈ S} a i. -/
def ofPoint {n : ℕ} (a : Fin n → R) : Coord n R := fun S => ∏ i ∈ S, a i

/-- Back to ordinary coordinates: read off the grade-1 corners. -/
def toPoint {n : ℕ} (c : Coord n R) : Fin n → R := fun i => c {i}

/-- The evaluation map: weights become one number at x. -/
def eval {n : ℕ} (x : Fin n → R) (c : Coord n R) : R :=
  ∑ S : Corner n, c S * ∏ i ∈ S, x i

/-- Box coordinates (the zeta transform): add up every corner inside T. -/
def zeta {n : ℕ} (c : Coord n R) : Coord n R := fun T => ∑ S ∈ T.powerset, c S

/-- The circle operation of two tesseract points: corners combine by joining their
switches (S ∪ U). This is the union, or OR, product on subsets. -/
def compose {n : ℕ} (v w : Coord n R) : Coord n R :=
  fun T => ∑ S : Corner n, ∑ U : Corner n, if S ∪ U = T then v S * w U else 0

/-- Hirota's N-soliton weights: corner S gets the product of the pair
interactions A i j for i < j inside S. -/
def solitonWeights {n : ℕ} (A : Fin n → Fin n → R) : Coord n R :=
  fun S => ∏ p ∈ (S ×ˢ S).filter (fun p => p.1 < p.2), A p.1 p.2

/-! ## Facts -/

/-- Nothing is lost going in and coming back out. -/
theorem toPoint_ofPoint {n : ℕ} (a : Fin n → R) : toPoint (ofPoint a) = a := by
  funext i; simp [toPoint, ofPoint]

/-- There are n choose k corners of grade k. -/
theorem card_corners_of_grade (n k : ℕ) :
    ((univ : Finset (Corner n)).filter (fun S => grade S = k)).card = n.choose k := by
  have h : (univ : Finset (Corner n)).filter (fun S => grade S = k)
      = powersetCard k (univ : Finset (Fin n)) := by
    rw [powersetCard_eq_filter, powerset_univ]; rfl
  rw [h, card_powersetCard, card_univ, Fintype.card_fin]

/-- The tesseract's layers are 1, 4, 6, 4, 1. -/
example : (powersetCard 0 (univ : Finset (Fin 4))).card = 1 := by
  rw [card_powersetCard, card_univ, Fintype.card_fin]; rfl
example : (powersetCard 1 (univ : Finset (Fin 4))).card = 4 := by
  rw [card_powersetCard, card_univ, Fintype.card_fin]; rfl
example : (powersetCard 2 (univ : Finset (Fin 4))).card = 6 := by
  rw [card_powersetCard, card_univ, Fintype.card_fin]; rfl
example : (powersetCard 3 (univ : Finset (Fin 4))).card = 4 := by
  rw [card_powersetCard, card_univ, Fintype.card_fin]; rfl
example : (powersetCard 4 (univ : Finset (Fin 4))).card = 1 := by
  rw [card_powersetCard, card_univ, Fintype.card_fin]; rfl

/-- The 16 pieces add up to the whole box: ∑_S ∏_{i ∈ S} a i = ∏ (1 + a i). -/
theorem sum_corners_ofPoint {n : ℕ} (a : Fin n → R) :
    ∑ S : Corner n, ofPoint a S = ∏ i, (1 + a i) := by
  rw [Finset.prod_one_add, powerset_univ]
  simp only [ofPoint]

/-- So the sum of all corners is 1 plus the circle operation of the coordinates. -/
theorem sum_corners_eq_one_add_circ {n : ℕ} (a : Fin n → R) :
    ∑ S : Corner n, ofPoint a S = 1 + circList (List.ofFn a) := by
  rw [sum_corners_ofPoint, one_add_circList, List.map_ofFn, List.prod_ofFn]
  rfl

/-- Evaluating "every weight 1" at x gives the box ∏ (1 + x i). -/
theorem eval_ones {n : ℕ} (x : Fin n → R) :
    eval x (fun _ => 1) = ∏ i, (1 + x i) := by
  rw [Finset.prod_one_add, powerset_univ]
  simp [eval]

/-- Box coordinates of a point are sub-box measures ∏_{i ∈ T} (1 + a i). -/
theorem zeta_ofPoint {n : ℕ} (a : Fin n → R) (T : Corner n) :
    zeta (ofPoint a) T = ∏ i ∈ T, (1 + a i) := by
  rw [Finset.prod_one_add]
  rfl

/-- In box coordinates, the circle operation axis by axis is plain multiplication. -/
theorem zeta_ofPoint_circ {n : ℕ} (a b : Fin n → R) (T : Corner n) :
    zeta (ofPoint (fun i => a i ◯ b i)) T = zeta (ofPoint a) T * zeta (ofPoint b) T := by
  simp only [zeta_ofPoint, one_add_circ, Finset.prod_mul_distrib]

/-- With no interaction (every A i j = 1) every soliton weight is 1, so the
tau function is the pure box ∏ (1 + e^{η i}). -/
theorem solitonWeights_one {n : ℕ} :
    solitonWeights (fun _ _ => (1 : R)) = fun (_ : Corner n) => 1 := by
  funext S; simp [solitonWeights]

/-
  Next step (not stated here yet): `zeta (compose v w) = zeta v * zeta w`
  for every pair of points, which makes `compose` exactly the circle operation in box
  coordinates. The Python helpers test it on many random points.
-/

end CircleOp
