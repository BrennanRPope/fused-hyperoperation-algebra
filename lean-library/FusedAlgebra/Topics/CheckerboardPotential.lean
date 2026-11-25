/-
  Fused Hyperoperation Algebra: adding a potential to the checkerboard / telegraph
  chain (topics/why-i/checkerboard-potential).
  Core Lean 4 only (no mathlib), no `sorry`, no imports.
  Library path: FusedAlgebra/Topics/CheckerboardPotential.lean (compiled, 3 theorems, standard axioms only).

  What is here (integer versions; the complex and limit statements are checks.py):
  1. a scalar potential element v·I commutes with any 2×2 turn element t, so one
     tick's factor (1 + k v)(1 + k t) equals 1 + k (v ∘ₖ t)   [A1 on the page];
  2. path by path, summing k^|S| ∏_{j∈S} v_j over every subset S of the ticks
     (every way to insert potential elements) gives ∏ (1 + k v_j)   [B2 on the page];
     this is the "each path still counted once" step: grouping the sub-chains of the
     2N-element chain by their turns leaves exactly this sum on each path.
-/
namespace CheckerboardPotential

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
def scalar (v : Int) : M2 := ⟨v, 0, 0, v⟩
def circ (k : Int) (x y : M2) : M2 := add (add x y) (smul k (mul x y))

/-- A potential element v·I commutes with any turn element. -/
theorem scalar_comm (v : Int) (t : M2) : mul (scalar v) t = mul t (scalar v) := by
  cases t; simp only [mul, scalar, M2.mk.injEq, Int.zero_mul, Int.mul_zero, Int.add_zero, Int.zero_add]
  exact ⟨Int.mul_comm _ _, Int.mul_comm _ _, Int.mul_comm _ _, Int.mul_comm _ _⟩

/-- One tick with a potential: (1 + k v)(1 + k t) = 1 + k (v ∘ₖ t). -/
theorem tick_factor (k v : Int) (t : M2) :
    mul (add one (smul k (scalar v))) (add one (smul k t)) =
    add one (smul k (circ k (scalar v) t)) := by
  cases t
  simp only [mul, add, smul, one, scalar, circ, M2.mk.injEq]
  -- core-only normalisation; with mathlib this whole line is `refine ⟨?_, ?_, ?_, ?_⟩ <;> ring`
  refine ⟨?_, ?_, ?_, ?_⟩ <;>
    simp only [Int.mul_add, Int.add_mul, Int.mul_one, Int.one_mul, Int.mul_zero, Int.zero_mul,
      Int.add_zero, Int.zero_add, Int.mul_assoc, Int.mul_comm, Int.mul_left_comm,
      Int.add_assoc, Int.add_comm, Int.add_left_comm]

/-! ## 2. Each path: the sum over insertion subsets is the product -/

/-- Sum over all subsets of the list of k^|S| ∏_{j∈S} v_j, built by choosing,
    for the head tick, "no potential kick" or "a kick". -/
def subsetSum (k : Int) : List Int → Int
  | [] => 1
  | v :: vs => subsetSum k vs + k * v * subsetSum k vs

def prodOnePlus (k : Int) : List Int → Int
  | [] => 1
  | v :: vs => (1 + k * v) * prodOnePlus k vs

theorem subsetSum_eq_prod (k : Int) (vs : List Int) : subsetSum k vs = prodOnePlus k vs := by
  induction vs with
  | nil => rfl
  | cons v vs ih =>
    simp only [subsetSum, prodOnePlus, ih]
    rw [Int.add_mul, Int.one_mul]

/-- Spot check at four ticks. -/
example : subsetSum 3 [2, -1, 5, 7] = prodOnePlus 3 [2, -1, 5, 7] := by decide

end CheckerboardPotential
