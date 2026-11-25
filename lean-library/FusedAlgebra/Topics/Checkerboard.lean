/-
  Fused Hyperoperation Algebra: Feynman's checkerboard and Kac's telegraph as one
  weighted-circle chain (topics/why-i/checkerboard).
  Core Lean 4 only (no mathlib), no `sorry`, no imports.
  Library path: FusedAlgebra/Topics/Checkerboard.lean (compiled, 4 theorems, standard axioms only).

  What is here (integer versions of the identities on the page; the complex and
  limit statements are checks in checks.py):
  1. on 2×2 matrices, which do not commute, 1 + k(x ∘ₖ y) = (1 + kx)(1 + ky);
  2. three steps expand into the 7 ordered sub-chains, each once, order kept;
  3. the swap X and a diagonal "move" do not commute (so the order matters);
  4. path counts: the number of N-step checkerboard paths with R turns ending at x
     matches the Feynman–Hibbs formula, for every x, R and end direction at N = 5.
-/
namespace Checkerboard

/-! ## 1. 2×2 integer matrices and the weighted circle operation -/

structure M2 where
  a : Int
  b : Int
  c : Int
  d : Int
deriving DecidableEq, Repr

def mul (x y : M2) : M2 :=
  ⟨x.a * y.a + x.b * y.c, x.a * y.b + x.b * y.d, x.c * y.a + x.d * y.c, x.c * y.b + x.d * y.d⟩
def add (x y : M2) : M2 := ⟨x.a + y.a, x.b + y.b, x.c + y.c, x.d + y.d⟩
def smul (k : Int) (x : M2) : M2 := ⟨k * x.a, k * x.b, k * x.c, k * x.d⟩
def one : M2 := ⟨1, 0, 0, 1⟩

/-- x ∘ₖ y = x + y + k·x·y, with the matrix product in that order. -/
def circ (k : Int) (x y : M2) : M2 := add (add x y) (smul k (mul x y))

/-- The chain rule of the weighted circle operation holds for non-commuting steps. -/
theorem circ_chain (k : Int) (x y : M2) :
    add one (smul k (circ k x y)) = mul (add one (smul k x)) (add one (smul k y)) := by
  cases x; cases y
  simp only [circ, add, smul, mul, one, M2.mk.injEq]
  refine ⟨?_, ?_, ?_, ?_⟩ <;> grind

/-- Three steps: every non-empty ordered sub-chain once, weight k^(size − 1). -/
theorem circ_three (k : Int) (x y z : M2) :
    circ k (circ k x y) z =
      add (add (add x y) z)
        (add (smul k (add (add (mul x y) (mul x z)) (mul y z))) (smul (k * k) (mul (mul x y) z))) := by
  cases x; cases y; cases z
  simp only [circ, add, smul, mul, M2.mk.injEq]
  refine ⟨?_, ?_, ?_, ?_⟩ <;> grind

/-! ## 2. The turn and the move do not commute -/

/-- The swap of the two directions (a turn). -/
def X : M2 := ⟨0, 1, 1, 0⟩
/-- A stand-in for the move in momentum space: different factors on the two directions. -/
def D : M2 := ⟨2, 0, 0, 3⟩

theorem turn_move_order : mul D X ≠ mul X D := by decide

/-! ## 3. Path counts (Feynman–Hibbs problem 2-6) -/

/-- Walk from position x in direction d (true = right). Each entry says whether to turn
    before the next move. Returns the end position and the last direction. -/
def walk : Bool → Int → List Bool → Int × Bool
  | d, x, [] => (x, d)
  | d, x, t :: ts =>
    let d' := if t then !d else d
    walk d' (if d' then x + 1 else x - 1) ts

def allLists : Nat → List (List Bool)
  | 0 => [[]]
  | n + 1 => (allLists n).flatMap (fun l => [false :: l, true :: l])

def nturns (l : List Bool) : Nat := (l.filter id).length

def count (N : Nat) (x : Int) (R : Nat) (dEnd : Bool) : Nat :=
  ((allLists N).filter (fun l => walk true 0 l == (x, dEnd) && nturns l == R)).length

def choose : Nat → Nat → Nat
  | _, 0 => 1
  | 0, _ + 1 => 0
  | n + 1, k + 1 => choose n k + choose n (k + 1)

/-- r counts the right moves plus the starting direction, l the left moves. -/
def formula (N : Nat) (x : Int) (R : Nat) (dEnd : Bool) : Nat :=
  if (((N : Int) + x) % 2 != 0) then 0 else
  let r := (((N : Int) + x) / 2 + 1).toNat
  let l := N + 1 - r
  if dEnd then
    if R % 2 == 1 then 0 else
    let s := R / 2
    if s == 0 then (if l == 0 then 1 else 0)
    else if r ≥ 1 && l ≥ 1 then choose (r - 1) s * choose (l - 1) (s - 1) else 0
  else
    if R % 2 == 0 then 0 else
    let s := (R - 1) / 2
    if r ≥ 1 && l ≥ 1 then choose (r - 1) s * choose (l - 1) s else 0

theorem counts_five :
    ∀ x ∈ (List.range 11).map (fun i : Nat => (i : Int) - 5), ∀ R ∈ List.range 6, ∀ dEnd ∈ [true, false],
      count 5 x R dEnd = formula 5 x R dEnd := by
  decide

end Checkerboard
