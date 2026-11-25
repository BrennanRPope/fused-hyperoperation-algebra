/-
  Graph theory through mult-add
  =============================

  Companion page: ../index.html.  Numbers: ../checks.py.
  Core Lean 4.37 only (no mathlib), no `sorry`, standard axioms only.
  Intended library path: FusedAlgebra/Topics/GraphTheory.lean
  Imports FusedAlgebra.Topics.Statistics for its landscape lens law `lens k q a b`
  and the Ali–Mikhail–Haq copula `amh`.

  1. Potts / random-cluster bonds (edge weight v, Q colours). Two bonds in parallel
     join with weight v₁ ◯ v₂ (mult-add); two in series with v₁v₂/(Q + v₁ + v₂).
     In the transmissivity t = v/(Q + v): series multiplies, t = t₁t₂, and parallel
     is the landscape law F_{Q−2, 1−Q}, on the "or" line k + q = −1. Q = 1 is the
     failure union (reliability), Q = 2 velocity addition (Ising, tanh), and the
     law is AMH union with θ = 1 − Q. Walls: a perfect bond t = 1, and the proper
     colouring bond v = −1, t = −1/(Q − 1).
  2. Two-terminal reliability: series ×, parallel ◯₋₁; the failure side swaps them;
     in odds-against coordinates series is plain mult-add. The Wheatstone bridge
     by pivoting.
  3. The diamond hierarchical lattice: the dial pair iterated. Percolation fixed
     point is the golden ratio (√5 − 1)/2; Ising's solves t³ + t² + t = 1; Q = 3's is 1/2.
  4. Products: closed adjacency (A + I) of the strong product is the Kronecker
     product, so entries, degrees and the count 2m + n follow mult-add; edge counts
     of P₃ ⊠ P₃ and C₄ ⊠ C₄ by computation.
  5. Independence polynomials: disjoint union mult-adds I − 1, join adds it;
     checked on C₄ = K̄₂ + K̄₂ by brute force over vertex subsets.

  Classical sources: Moore and Shannon 1956; Fortuin and Kasteleyn 1972; Sokal 2005;
  Berker and Ostlund 1979; Ali, Mikhail and Haq 1978; Gutman and Harary 1983.
  Ours: the placement on the landscape's "or" line and the AMH reading θ = 1 − Q.
  LIBRARY NOTE: in the shared library at FusedAlgebra/Topics/GraphTheory.lean; four exact twins removed (see the LIBRARY comments).
-/
import FusedAlgebra.Base
import FusedAlgebra.Topics.Statistics

namespace GraphTheory
open CircleOp
open StatFused (lens amh)

open scoped CircleOp

set_option linter.unusedSectionVars false

/-! # 1. Potts bonds: parallel mult-adds, series multiplies transmissivity -/

/-- Two parallel bonds join their ends with weight v₁ + v₂ + v₁v₂: one, the other, or both. -/
theorem parallel_join (v₁ v₂ : Rat) : v₁ + v₂ + v₁ * v₂ = v₁ ◯ v₂ := rfl

/-- Series pair through a middle vertex: "ends not joined" weighs Q + v₁ + v₂ (middle alone,
    or tied to one end), "ends joined" weighs v₁v₂. One bond of weight v with the factor
    Q + v₁ + v₂ gives the same two numbers. -/
theorem series_local (Q v₁ v₂ : Rat) (h : Q + v₁ + v₂ ≠ 0) :
    (Q + v₁ + v₂) * 1 = Q + v₁ + v₂ ∧ (Q + v₁ + v₂) * (v₁ * v₂ / (Q + v₁ + v₂)) = v₁ * v₂ := by
  constructor
  · grind
  · grind

theorem add_div' (Q x S : Rat) (hS : S ≠ 0) : Q + x / S = (Q * S + x) / S := by grind

theorem ddiv (x y z : Rat) (hz : z ≠ 0) (hy : y ≠ 0) : x / z / (y / z) = x / y := by grind

theorem cancel_left (Q W : Rat) (hQ : Q ≠ 0) (h : Q + W ≠ 0) :
    Q * W / (Q * (Q + W)) = W / (Q + W) := by grind

/-- Transmissivity of a bond. -/
def tr (Q v : Rat) : Rat := v / (Q + v)

/-- **Series multiplies transmissivity.** -/
theorem series_tr (Q v₁ v₂ : Rat) (h1 : Q + v₁ ≠ 0) (h2 : Q + v₂ ≠ 0) (h : Q + v₁ + v₂ ≠ 0) :
    tr Q (v₁ * v₂ / (Q + v₁ + v₂)) = tr Q v₁ * tr Q v₂ := by
  unfold tr
  rw [add_div' _ _ _ h]
  have e : Q * (Q + v₁ + v₂) + v₁ * v₂ = (Q + v₁) * (Q + v₂) := by grind
  rw [e]; grind

/-- **Parallel is a landscape law.** With tᵢ = vᵢ/(Q + vᵢ), the joined bond v₁ ◯ v₂ has
    t = (t₁ + t₂ + (Q − 2)t₁t₂)/(1 − (1 − Q)t₁t₂) = F_{Q−2, 1−Q}(t₁, t₂). -/
theorem parallel_tr (Q v₁ v₂ : Rat) (h1 : Q + v₁ ≠ 0) (h2 : Q + v₂ ≠ 0)
    (h : Q + (v₁ ◯ v₂) ≠ 0) (hQ : Q ≠ 0) :
    tr Q (v₁ ◯ v₂) = lens (Q - 2) (1 - Q) (tr Q v₁) (tr Q v₂) := by
  unfold tr lens
  have hp : (Q + v₁) * (Q + v₂) ≠ 0 := by grind
  have hn : v₁ / (Q + v₁) + v₂ / (Q + v₂) + (Q - 2) * (v₁ / (Q + v₁)) * (v₂ / (Q + v₂)) =
      Q * (v₁ ◯ v₂) / ((Q + v₁) * (Q + v₂)) := by
    unfold circ; grind
  have hd : 1 - (1 - Q) * (v₁ / (Q + v₁)) * (v₂ / (Q + v₂)) =
      Q * (Q + (v₁ ◯ v₂)) / ((Q + v₁) * (Q + v₂)) := by
    have e : (1 - Q) * (v₁ / (Q + v₁)) * (v₂ / (Q + v₂)) =
        (1 - Q) * v₁ * v₂ / ((Q + v₁) * (Q + v₂)) := by grind
    rw [e]
    unfold circ; grind
  rw [hn, hd, ddiv _ _ _ hp (by grind)]
  exact (cancel_left Q _ hQ h).symm

-- LIBRARY: or_line (every Q lands on the line k + q = −1) was an exact twin of LogicSets.hamacher_line and is removed.

/-- The lift of the parallel law is 1 + v: (1 + (Q − 1)t)/(1 − t) = 1 + v. -/
theorem parallel_lift (Q v : Rat) (h : Q + v ≠ 0) (hQ : Q ≠ 0) :
    (1 + (Q - 1) * tr Q v) / (1 - tr Q v) = 1 + v := by
  unfold tr
  have e : 1 - v / (Q + v) = Q / (Q + v) := by grind
  rw [e]; grind

/-- Q = 1: percolation and reliability, the failure union a ◯₋₁ b. -/
theorem q_one (a b : Rat) : lens (1 - 2) (1 - 1) a b = circK (-1) a b := by
  unfold lens circK; grind

/-- Q = 2: Ising, velocity addition (tanh addition). -/
theorem q_two (a b : Rat) : lens (2 - 2) (1 - 2) a b = (a + b) / (1 + a * b) := by
  unfold lens; grind

/-- Q = 0: the spanning-forest end, the cusp k² = 4q of the landscape. -/
theorem q_zero_cusp : ((0 : Rat) - 2) * (0 - 2) = 4 * (1 - 0) := by grind

/-- **Potts parallel is AMH union**, with dependence θ = 1 − Q: the two bonds "fail"
    together like AMH-dependent events. -/
theorem potts_is_amh (Q a b : Rat) (h : 1 - (1 - Q) * a * b ≠ 0) :
    lens (Q - 2) (1 - Q) a b = 1 - amh (1 - Q) (1 - a) (1 - b) := by
  rw [StatFused.amh_union_is_lens (1 - Q) a b h]
  have e : (-(1 + (1 - Q))) = Q - 2 := by grind
  rw [e]

/-- Wall 1: a perfect bond t = 1 swallows anything put in parallel with it. -/
theorem perfect_bond_absorbs (Q a : Rat) (h : 1 - (1 - Q) * a * 1 ≠ 0) :
    lens (Q - 2) (1 - Q) a 1 = 1 := by
  unfold lens at *; grind

/-- The proper-colouring bond v = −1 has t = −1/(Q − 1). -/
theorem colouring_bond (Q : Rat) : tr Q (-1) = -1 / (Q - 1) := by
  unfold tr; grind

/-- Wall 2: that bond swallows anything in parallel (a forbidden pair stays forbidden). -/
theorem colouring_wall_absorbs (Q a : Rat) (h : Q - 1 ≠ 0)
    (hd : 1 - (1 - Q) * a * (-1 / (Q - 1)) ≠ 0) :
    lens (Q - 2) (1 - Q) a (-1 / (Q - 1)) = -1 / (Q - 1) := by
  unfold lens at *
  have e : 1 - (1 - Q) * a * (-1 / (Q - 1)) = 1 - a := by grind
  rw [e] at hd ⊢
  grind

/-- In weights, the colouring bond is mult-add's own wall −1. -/
theorem colouring_is_wall (v : Rat) : (-1) ◯ v = -1 := neg_one_absorbs v

/-! # 2. Two-terminal reliability: the dial pair (×, ◯₋₁) -/

-- LIBRARY: parallel_rel (parallel works if either works) was an exact twin of InfoFused.union_rule (also ConsensusFused.either) and is removed.

-- LIBRARY: series_fails (a series pair fails by ◯₋₁) was an exact twin of ThermoFused.ratio_series (also InfoFused.kraft_slack) and is removed.

/-- In odds against, u = 1/p − 1, series is plain mult-add. -/
theorem series_odds (p₁ p₂ : Rat) (h1 : p₁ ≠ 0) (h2 : p₂ ≠ 0) :
    1 / (p₁ * p₂) - 1 = (1 / p₁ - 1) ◯ (1 / p₂ - 1) := by
  unfold circ; grind

/-- **The Wheatstone bridge** has no series or parallel pair; pivot on the middle bond:
    contracted, it is two parallel pairs in series; deleted, two series pairs in parallel. -/
theorem wheatstone (p : Rat) :
    p * (circK (-1) p p * circK (-1) p p) + (1 - p) * circK (-1) (p * p) (p * p) =
      2 * p ^ 2 + 2 * p ^ 3 - 5 * p ^ 4 + 2 * p ^ 5 := by
  unfold circK; grind

/-! # 3. The diamond lattice: the dial pair iterated -/

/-- One diamond step: two series pairs in parallel. -/
def diamond (p : Rat) : Rat := circK (-1) (p * p) (p * p)

theorem diamond_poly (p : Rat) : diamond p = 2 * p ^ 2 - p ^ 4 := by
  unfold diamond circK; grind

/-- **The golden fixed point.** If p² = 1 − p (p = (√5 − 1)/2), the diamond step fixes p. -/
theorem golden_fixed (p : Rat) (h : p * p = 1 - p) : diamond p = p := by
  unfold diamond circK
  have h4 : p * p * (p * p) = 2 - 3 * p := by rw [h]; grind
  grind

/-- Slope there: 4p − 4p³ = 4 − 4p (≈ 1.528 > 1, so the fixed point repels). -/
theorem golden_slope (p : Rat) (h : p * p = 1 - p) : 4 * p - 4 * (p * p * p) = 4 - 4 * p := by
  have h3 : p * p * p = p - (1 - p) := by
    have : p * p * p = p * (p * p) := by grind
    rw [this, h]; grind
  rw [h3]; grind

/-- The Potts diamond step at Q colours. -/
def diamondQ (Q t : Rat) : Rat := lens (Q - 2) (1 - Q) (t * t) (t * t)

theorem diamondQ_one (t : Rat) : diamondQ 1 t = diamond t := by
  unfold diamondQ diamond lens circK; grind

/-- Ising (Q = 2): if t³ + t² + t = 1, the diamond step fixes t. -/
theorem ising_fixed (t : Rat) (h : t * t * t + t * t + t = 1) : diamondQ 2 t = t := by
  unfold diamondQ lens
  have h4 : t * t * (t * t) = 2 * t - 1 := by
    have e : t * t * (t * t) = t * (t * t * t) := by grind
    rw [e]
    have e2 : t * t * t = 1 - t * t - t := by grind
    rw [e2]; grind
  have hd : 1 - (1 - 2) * (t * t) * (t * t) = 2 * t := by grind
  rw [hd]
  have ht : t ≠ 0 := by intro h0; rw [h0] at h; grind
  grind

/-- Three colours: the diamond fixed point is exactly 1/2. -/
theorem three_state_fixed : diamondQ 3 (1 / 2) = 1 / 2 := by
  unfold diamondQ lens; grind

/-! # 4. Products: the strong product is the Kronecker mult-add -/

/-- Closed adjacency entries multiply: (a + i)(b + j) − i·j = a·j + i·b + a·b. With
    i = j = 1 (same vertex in both) this is a ◯ b; the three terms are the Cartesian
    pieces A⊗I, I⊗B and the tensor corner A⊗B. -/
theorem strong_entry (a b i j : Int) : (a + i) * (b + j) - i * j = a * j + i * b + a * b := by
  grind

theorem strong_entry_diag (a b : Int) : (a + 1) * (b + 1) - 1 = a ◯ b := by
  unfold circ; grind

/-- Edge counts: 2m + n multiplies. -/
theorem strong_count (n₁ m₁ n₂ m₂ : Nat) :
    2 * (n₁ * m₂ + n₂ * m₁ + 2 * m₁ * m₂) + n₁ * n₂ = (2 * m₁ + n₁) * (2 * m₂ + n₂) := by
  grind

/-- A small graph: vertex count and adjacency test. -/
structure SG where
  n : Nat
  adj : Nat → Nat → Bool

def edgeCount (g : SG) : Nat :=
  ((List.range g.n).map (fun i => ((List.range g.n).filter (fun j => i < j && g.adj i j)).length)).foldl (· + ·) 0

def path (n : Nat) : SG := ⟨n, fun i j => i + 1 == j || j + 1 == i⟩
def cycle (n : Nat) : SG := ⟨n, fun i j => i != j && ((i + 1) % n == j || (j + 1) % n == i)⟩

/-- Strong product: different pairs that are equal-or-adjacent in both coordinates. -/
def strong (g h : SG) : SG :=
  ⟨g.n * h.n, fun x y =>
    let (a, b) := (x / h.n, x % h.n); let (c, d) := (y / h.n, y % h.n)
    x != y && (a == c || g.adj a c) && (b == d || h.adj b d)⟩

theorem p3_strong_p3 : edgeCount (strong (path 3) (path 3)) = 20 := by decide

theorem c4_strong_c4 : edgeCount (strong (cycle 4) (cycle 4)) = 4 * 4 + 4 * 4 + 2 * 4 * 4 := by decide

theorem c5_strong_p2 : edgeCount (strong (cycle 5) (path 2)) = 5 * 1 + 2 * 5 + 2 * 5 * 1 := by decide

/-! # 5. Independence polynomials: union mult-adds, join adds -/

-- LIBRARY: union_circ (disjoint union multiplies independence polynomials, so I − 1 mult-adds) was an exact twin of SolitonsRungs.glue_is_circ and is removed.

/-- A join lets a set use only one side, so I − 1 adds (dial k = 0). -/
theorem join_adds (IG IH : Int) : (IG + IH - 1) - 1 = circK 0 (IG - 1) (IH - 1) := by
  unfold circK; grind

/-- Coefficients of the independence polynomial by brute force over vertex subsets. -/
def subsetsOf : List Nat → List (List Nat)
  | [] => [[]]
  | a :: l => subsetsOf l ++ (subsetsOf l).map (a :: ·)

def isIndep (g : SG) (S : List Nat) : Bool :=
  S.all (fun i => S.all (fun j => !(g.adj i j)))

def indepPoly (g : SG) : List Nat :=
  let ind := (subsetsOf (List.range g.n)).filter (isIndep g)
  (List.range (g.n + 1)).map (fun r => (ind.filter (fun S => S.length == r)).length)

/-- C₄ is the join of two non-edges: I = (1 + x)² + (1 + x)² − 1 = 1 + 4x + 2x². -/
theorem c4_indep : indepPoly (cycle 4) = [1, 4, 2, 0, 0] := by decide

/-- The disjoint union of two edges, 2K₂ = the complement of C₄: I = (1 + 2x)². -/
def twoK2 : SG := ⟨4, fun i j => (i == 0 && j == 1) || (i == 1 && j == 0) || (i == 2 && j == 3) || (i == 3 && j == 2)⟩

theorem two_k2_indep : indepPoly twoK2 = [1, 4, 4, 0, 0] := by decide

theorem c5_indep : indepPoly (cycle 5) = [1, 5, 5, 0, 0, 0] := by decide

end GraphTheory
