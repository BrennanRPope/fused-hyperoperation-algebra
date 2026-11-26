/-
  FuzzyFrame: the finite facts behind the fuzzy-cut monogamy proof   (DRAFT, ours)
  First compiled in the shared library on Nov 26. Core Lean 4.37, `decide` only
  (`moments_nonneg` uses `decide +kernel`, which the elaborator's `decide` hit a recursion limit on).

  The proof (write-up, part 3) has three finite ingredients, checked here:
    1. the odd frame v₁..v₄ is orthogonal with |vᵢ|² = 4;
    2. negating copy 4 sends the Hadamard frame to (v₁, v₂, v₃, −v₄);
    3. every moment of the rotated single-site set S' is ≥ 0
       (here for all exponents up to 3, which covers every parity pattern;
       on S' a moment only depends on the exponents' parities and on
       whether a single coordinate is used, so this is the whole content).
  The analytic step, Griffiths' first inequality for single-site measures
  with nonnegative moments, is a standard theorem (Ginibre 1970) and is
  not formalized here.
-/

namespace FuzzyFrame

def dot (a b : List Int) : Int := (List.zipWith (· * ·) a b).foldl (· + ·) 0

def oddFrame : List (List Int) := [[1,1,1,-1], [1,1,-1,1], [1,-1,1,1], [-1,1,1,1]]
def hadamard : List (List Int) := [[1,1,1,1], [1,1,-1,-1], [1,-1,1,-1], [1,-1,-1,1]]

theorem oddFrame_orthogonal :
    (oddFrame.map fun a => oddFrame.map fun b => dot a b)
      = [[4,0,0,0], [0,4,0,0], [0,0,4,0], [0,0,0,4]] := by decide

/-- Negate copy 4 (the last entry of each pinned vector). -/
def negLast (a : List Int) : List Int := a.take 3 ++ (a.drop 3).map (fun t => -t)

theorem flip_copy4 :
    hadamard.map negLast = [[1,1,1,-1], [1,1,-1,1], [1,-1,1,1], [1,-1,-1,-1]] := by decide

theorem last_is_neg_v4 : [1,-1,-1,-1] = ([-1,1,1,1] : List Int).map (fun t => -t) := by decide

def signs : List (List Int) :=
  [[1,1,1,1],[1,1,1,-1],[1,1,-1,1],[1,1,-1,-1],[1,-1,1,1],[1,-1,1,-1],[1,-1,-1,1],[1,-1,-1,-1],
   [-1,1,1,1],[-1,1,1,-1],[-1,1,-1,1],[-1,1,-1,-1],[-1,-1,1,1],[-1,-1,1,-1],[-1,-1,-1,1],[-1,-1,-1,-1]]

/-- Rotated coordinates αᵢ(σ) = vᵢ·σ / 2. -/
def alpha (s : List Int) : List Int := oddFrame.map fun v => dot v s / 2

def monomial (a : List Int) (n : List Nat) : Int :=
  (List.zipWith (fun (x : Int) (k : Nat) => x ^ k) a n).foldl (· * ·) 1

def moment (n : List Nat) : Int := (signs.map fun s => monomial (alpha s) n).foldl (· + ·) 0

def exps : List (List Nat) :=
  (List.range 4).flatMap fun a => (List.range 4).flatMap fun b =>
    (List.range 4).flatMap fun c => (List.range 4).map fun d => [a, b, c, d]

theorem moments_nonneg : (exps.all fun n => decide (0 ≤ moment n)) = true := by decide +kernel

/-- A pinned terminal j sits at αⱼ = 2, other coordinates 0. -/
theorem pinned_alpha : oddFrame.map alpha = [[2,0,0,0], [0,2,0,0], [0,0,2,0], [0,0,0,2]] := by decide

end FuzzyFrame
