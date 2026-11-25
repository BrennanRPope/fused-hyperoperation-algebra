/-
  Combinatorics through mult-add
  ==============================

  Companion page: ../index.html.  Numbers: ../checks.py.
  Core Lean 4.37 only (no mathlib), no `sorry`, standard axioms only.
  Intended library path: FusedAlgebra/Topics/Combinatorics.lean

  1. The cube counts subsets. n copies of 1 mult-added give 2ⁿ − 1 (the nonempty
     subsets); the dial k paints each chosen element one of k colours, (1 + k)ⁿ.
     k = −1 makes the alternating sum vanish, k = −2 gives (−1)ⁿ.
  2. The dial fold of a list, and the wall: one entry on the wall −1/k sends the
     whole fold to the wall.
  3. n! = 1 + (1 ◯ 2 ◯ ⋯ ◯ (n−1)). The cube with axis weights 1, …, n−1 has layers
     equal to the Stirling cycle numbers (permutations sorted by cycles). At k = −1
     the entry 1 sits on the wall, so even and odd permutations balance.
  4. Derangements by inclusion–exclusion (k = −1 on fixed points).
  5. The exponential formula as a recurrence: E applied to "connected" pieces.
     E(nonempty sets) = Bell numbers, E(Bell) = partitions of partitions,
     E(connected graphs) = all graphs.
  6. Rung flow on x: the coefficients a₂(t), a₃(t), a₄(t) of Eᵗ(x) obey the flow's
     group law exactly; t = 1 is e^x − 1, t = 2 gives Bell, t = 1/2 gives a zero
     coefficient and t = 1/4 a negative one.
  7. Stirling numbers of the second kind, the dial exponential, compositions.

  All classical (Euler, Stirling, de Montmort, Pólya, Riordan, Comtet); what is
  ours is only the reading on the cube and the rung flow.
-/
import FusedAlgebra.Base
import FusedAlgebra.Cube

namespace Combinatorics
open CircleOp
open CircleOp.Cube (layer fact)

open scoped CircleOp

set_option linter.unusedSectionVars false

/-! # 1. The cube counts subsets -/

/-- The corners of the n-cube: 2ⁿ sub-lists. -/
theorem groups_length {α : Type} (l : List α) : (groups l).length = 2 ^ l.length := by
  induction l with
  | nil => rfl
  | cons a l ih =>
    simp only [groups, List.length_append, List.length_map, ih, List.length_cons, Nat.pow_succ]
    omega

section Semiring
variable {R : Type} [Lean.Grind.CommSemiring R]

/-- n copies of a, folded by mult-add, is the mult-add power. -/
theorem foldCirc_replicate (a : R) (n : Nat) : foldCirc (List.replicate n a) = mpow a n := by
  induction n with
  | zero => rfl
  | succ n ih => rw [List.replicate_succ, foldCirc, ih, mpow, circ_comm]

/-- **The dial paints subsets.** The corners of the n-cube, each switched-on axis
    weighing k, add up to (1 + k)ⁿ: every element is off, or on in one of k colours. -/
theorem painted_subsets (k : R) (n : Nat) :
    sumL ((groups (List.replicate n k)).map prodL) = (1 + k) ^ n := by
  rw [← one_add_foldCirc_eq_corners, foldCirc_replicate, mpow_lift]

/-- Mult-add of a with itself n times, at a = 1: the lift of 2ⁿ − 1 nonempty subsets. -/
theorem nonempty_subsets_lift (n : Nat) : 1 + foldCirc (List.replicate n (1 : R)) = 2 ^ n := by
  rw [foldCirc_replicate, mpow_lift]
  have : (1 : R) + 1 = 2 := by grind
  rw [this]

end Semiring

/-- On Int: 2ⁿ − 1 nonempty subsets. -/
theorem nonempty_subsets (n : Nat) : foldCirc (List.replicate n (1 : Int)) = 2 ^ n - 1 := by
  have h := nonempty_subsets_lift (R := Int) n
  omega

/-- k = −1: the signed corners cancel (inclusion–exclusion's engine). -/
theorem alternating_corners (n : Nat) :
    sumL ((groups (List.replicate (n + 1) (-1 : Int))).map prodL) = 0 := by
  rw [painted_subsets]; simp [Int.pow_succ]

/-- k = −2: the parity dial gives (−1)ⁿ. -/
theorem parity_corners (n : Nat) :
    sumL ((groups (List.replicate n (-2 : Int))).map prodL) = (-1) ^ n := by
  rw [painted_subsets]; rfl

/-! # 2. The dial fold and its wall -/

section Ring
variable {R : Type} [Lean.Grind.CommRing R]

/-- x₁ ◯ₖ x₂ ◯ₖ ⋯ ◯ₖ xₙ. -/
def wfold (k : R) : List R → R
  | [] => 0
  | a :: l => circK k a (wfold k l)

/-- The dial lift of a whole list. -/
theorem wfold_lift (k : R) (l : List R) :
    1 + k * wfold k l = prodL (l.map (fun x => 1 + k * x)) := by
  induction l with
  | nil => simp only [wfold, List.map, prodL]; grind
  | cons a l ih => simp only [wfold, List.map, prodL]; rw [← ih, circK_lift]

theorem prodL_eq_zero_of_mem (l : List R) (a : R) (h : a ∈ l) (h0 : a = 0) :
    prodL l = 0 := by
  induction l with
  | nil => cases h
  | cons b l ih =>
    simp only [prodL]
    rcases List.mem_cons.mp h with h1 | h1
    · rw [← h1, h0]; grind
    · rw [ih h1]; grind

/-- **The wall swallows the fold.** If one entry sits on the wall (1 + k·a = 0),
    the lift of the whole fold is 0. -/
theorem wall_entry_absorbs (k : R) (l : List R) (a : R) (h : a ∈ l) (hw : 1 + k * a = 0) :
    1 + k * wfold k l = 0 := by
  rw [wfold_lift]
  apply prodL_eq_zero_of_mem _ (1 + k * a)
  · exact List.mem_map.mpr ⟨a, h, rfl⟩
  · exact hw

end Ring

/-! # 3. n! as a mult-add fold; Stirling cycle numbers are the cube's layers -/

/-- 1 ◯ 2 ◯ ⋯ ◯ m. -/
def ffold : Nat → Nat
  | 0 => 0
  | m + 1 => ffold m ◯ (m + 1)

/-- **n! − 1 is a mult-add fold**: 1 + (1 ◯ 2 ◯ ⋯ ◯ m) = (m + 1)!. -/
theorem ffold_fact (m : Nat) : ffold m + 1 = fact (m + 1) := by
  induction m with
  | zero => rfl
  | succ m ih =>
    rw [ffold, circ_lift', ih,
      show fact (m + 1 + 1) = (m + 1 + 1) * fact (m + 1) from rfl, Nat.mul_comm]

theorem ffold_small : (List.range 8).map ffold = [0, 1, 5, 23, 119, 719, 5039, 40319] := by decide

/-- The dial version: 1 ◯ₖ 2 ◯ₖ ⋯ ◯ₖ m. -/
def kfold (k : Int) : Nat → Int
  | 0 => 0
  | m + 1 => circK k (kfold k m) (m + 1)

/-- Rising product (1 + k)(1 + 2k)⋯(1 + mk). -/
def rprod (k : Int) : Nat → Int
  | 0 => 1
  | m + 1 => rprod k m * (1 + k * (m + 1))

theorem kfold_lift (k : Int) (m : Nat) : 1 + k * kfold k m = rprod k m := by
  induction m with
  | zero => simp [kfold, rprod]
  | succ m ih => rw [kfold, rprod, circK_lift, ih]

/-- **Signs of permutations balance.** At k = −1 the entry 1 is on the wall, so the
    signed count Π(1 − i) of permutations of m + 2 things is 0. -/
theorem sign_balance (m : Nat) : 1 + (-1) * kfold (-1) (m + 1) = 0 := by
  rw [kfold_lift]
  induction m with
  | zero => decide
  | succ m ih => rw [rprod, ih]; simp

/-- Stirling cycle numbers c(n, j): permutations of n things with j cycles.
    Insert n + 1 as a new cycle, or after one of the n things already placed. -/
def cyc : Nat → Nat → Nat
  | 0, 0 => 1
  | 0, _ + 1 => 0
  | _ + 1, 0 => 0
  | n + 1, j + 1 => n * cyc n (j + 1) + cyc n j

/-- Layers of the weighted cube: coefficients of (1 + a₁x)(1 + a₂x)⋯ (corner sums by grade). -/
def cubeLayers : List Nat → List Nat
  | [] => [1]
  | a :: l =>
    let c := cubeLayers l
    List.zipWith (· + ·) (c ++ [0]) (0 :: c.map (a * ·))

def row (n : Nat) : List Nat := (List.range (n + 1)).map (cyc n)

/-- **The cube with axis weights 1, …, n−1 has the Stirling cycle numbers as layers**
    (grade j = permutations with n − j cycles), checked for n = 1 … 8. -/
theorem cube_is_stirling :
    (List.range 8).map (fun m => cubeLayers ((List.range m).map (· + 1))) =
    (List.range 8).map (fun m => (row (m + 1)).reverse.take (m + 1)) := by decide

theorem stirling_row4 : row 4 = [0, 6, 11, 6, 1] := by decide

/-- Row sums are n!. -/
theorem cyc_row_sums :
    (List.range 9).map (fun n => (row n).foldl (· + ·) 0) = (List.range 9).map fact := by decide

/-! # 4. Derangements: k = −1 on fixed points -/

/-- Derangements by Euler's recurrence. -/
def der : Nat → Nat
  | 0 => 1
  | 1 => 0
  | n + 2 => (n + 1) * (der (n + 1) + der n)

/-- Inclusion–exclusion: Σⱼ (−1)ʲ C(n, j)(n − j)!. -/
def derIE (n : Nat) : Int :=
  ((List.range (n + 1)).map (fun j => (-1 : Int) ^ j * (layer n j : Int) * (fact (n - j) : Int))).foldl (· + ·) 0

theorem derangements_incl_excl :
    (List.range 11).map (fun n => (der n : Int)) = (List.range 11).map derIE := by decide

theorem derangements_values : (List.range 9).map der = [1, 0, 1, 2, 9, 44, 265, 1854, 14833] := by decide

/-! # 5. The exponential formula: E applied to connected pieces -/

/-- Sum of a list of naturals. -/
def nsum (l : List Nat) : Nat := l.foldl (· + ·) 0

/-- Structures made of a nonempty-or-empty set of connected pieces, by the rooted
    recurrence g(n+1) = Σⱼ C(n, j)·c(j+1)·g(n − j): the piece holding the last point
    has j + 1 points. The list holds g(0), …, g(N). This is 1 + E(C). -/
def expList (c : Nat → Nat) : Nat → List Nat
  | 0 => [1]
  | n + 1 =>
    let g := expList c n
    g ++ [nsum ((List.range (n + 1)).map (fun j => layer n j * c (j + 1) * g.getD (n - j) 0))]

/-- E(nonempty sets) = all set partitions: the Bell numbers (rung 2 of the flow on x). -/
theorem bell_is_E_twice :
    expList (fun _ => 1) 9 = [1, 1, 2, 5, 15, 52, 203, 877, 4140, 21147] := by decide

def bellF (n : Nat) : Nat := (expList (fun _ => 1) 9).getD n 0

/-- E(Bell) = partitions of partitions (rung 3 of the flow on x). -/
theorem rung_three : expList bellF 8 = [1, 1, 3, 12, 60, 358, 2471, 19302, 167894] := by decide

def connF (n : Nat) : Nat := [0, 1, 1, 4, 38, 728, 26704, 1866256].getD n 0

/-- E(connected graphs) = all graphs: 2^C(n,2). -/
theorem graphs_from_connected :
    expList connF 7 = (List.range 8).map (fun n => 2 ^ (n * (n - 1) / 2)) := by decide

/-- Labelled mult-add: (1 + A)(1 + B) − 1 at n is Aₙ + Bₙ + Σ C(n, j)·Aⱼ·B_{n−j}
    over splits with both parts nonempty. For A = linear orders, B = sets this counts
    arrangements Σ C(n, j)·j!. -/
def lsplit (a b : Nat → Nat) (n : Nat) : Nat :=
  a n + b n + nsum ((List.range (n + 1)).map (fun j =>
    if j = 0 ∨ j = n then 0 else layer n j * a j * b (n - j)))

theorem arrangements :
    (List.range 8).map (lsplit fact (fun _ => 1)) = [2, 2, 5, 16, 65, 326, 1957, 13700] := by decide

/-! # 6. Rung flow on x: the coefficients of Eᵗ(x) -/

/-- Eᵗ(x) = x + a₂(t)x² + a₃(t)x³ + a₄(t)x⁴ + ⋯ -/
def a2 (t : Rat) : Rat := t / 2
def a3 (t : Rat) : Rat := (3 * t * t - t) / 12
def a4 (t : Rat) : Rat := t * (2 * t - 1) * (3 * t - 1) / 48

/-- Group law, order 2: Eˢ ∘ Eᵗ = Eˢ⁺ᵗ. -/
theorem flow_order2 (s t : Rat) : a2 t + a2 s = a2 (s + t) := by unfold a2; grind

/-- Group law, order 3: x³ coefficient of f(g(x)) is b₃ + 2a₂b₂ + a₃. -/
theorem flow_order3 (s t : Rat) : a3 t + 2 * a2 s * a2 t + a3 s = a3 (s + t) := by
  unfold a2 a3; grind

/-- Group law, order 4: x⁴ coefficient of f(g(x)) is b₄ + a₂(2b₃ + b₂²) + 3a₃b₂ + a₄. -/
theorem flow_order4 (s t : Rat) :
    a4 t + a2 s * (2 * a3 t + a2 t * a2 t) + 3 * a3 s * a2 t + a4 s = a4 (s + t) := by
  unfold a2 a3 a4; grind

/-- t = 1 is e^x − 1: 1/2, 1/6, 1/24. -/
theorem flow_at_one : a2 1 = 1 / 2 ∧ a3 1 = 1 / 6 ∧ a4 1 = 1 / 24 := by
  unfold a2 a3 a4; refine ⟨?_, ?_, ?_⟩ <;> grind

/-- t = 2 counts set partitions: 2!·a₂ = 2, 3!·a₃ = 5, 4!·a₄ = 15. -/
theorem flow_at_two : 2 * a2 2 = 2 ∧ 6 * a3 2 = 5 ∧ 24 * a4 2 = 15 := by
  unfold a2 a3 a4; refine ⟨?_, ?_, ?_⟩ <;> grind

/-- t = 3 counts partitions of partitions: 2, 12, 60 (with the factorials). -/
theorem flow_at_three : 2 * a2 3 = 3 ∧ 6 * a3 3 = 12 ∧ 24 * a4 3 = 60 := by
  unfold a2 a3 a4; refine ⟨?_, ?_, ?_⟩ <;> grind

/-- t = −1 is log(1 + x): signed cyclic orders −1, 2, −6. -/
theorem flow_at_minus_one : 2 * a2 (-1) = -1 ∧ 6 * a3 (-1) = 2 ∧ 24 * a4 (-1) = -6 := by
  unfold a2 a3 a4; refine ⟨?_, ?_, ?_⟩ <;> grind

/-- The half rung has a zero x⁴ coefficient ... -/
theorem half_rung_x4 : a4 (1 / 2) = 0 := by unfold a4; grind

/-- ... and a quarter rung has a negative x³ coefficient: fractional rungs do not count. -/
theorem quarter_rung_negative : a3 (1 / 4) < 0 := by unfold a3; decide +kernel

/-! # 7. Stirling numbers of the second kind; compositions -/

/-- S(n, j): set partitions of n things into j blocks. -/
def S2 : Nat → Nat → Nat
  | 0, 0 => 1
  | 0, _ + 1 => 0
  | _ + 1, 0 => 0
  | n + 1, j + 1 => (j + 1) * S2 n (j + 1) + S2 n j

/-- Row sums are the Bell numbers, the t = 2 rung. -/
theorem S2_rows_bell :
    (List.range 10).map (fun n => nsum ((List.range (n + 1)).map (S2 n))) =
    expList (fun _ => 1) 9 := by decide

/-- The dial exponential at k: Σⱼ S(n, j)·k^(j−1) (blocks weighed by k). At k = −1. -/
def dialBlocks (k : Int) (n : Nat) : Int :=
  ((List.range n).map (fun j => (S2 n (j + 1) : Int) * k ^ j)).foldl (· + ·) 0

theorem dial_blocks_minus_one :
    (List.range 8).map (fun n => dialBlocks (-1) (n + 1)) = [1, 0, -1, -1, 2, 9, 9, -50] := by decide

/-- Compositions on the cube: m parts, choose which j are nonzero (a corner of the
    m-cube), then a composition of n into j positive parts. Vandermonde. -/
theorem compositions_cube :
    ∀ m : Fin 5, ∀ n : Fin 9,
      nsum ((List.range (m.val + 2)).map (fun j =>
        if j = 0 then 0 else layer (m.val + 1) j * layer n.val (j - 1)))
        = layer (n.val + m.val + 1) m.val := by decide

end Combinatorics
