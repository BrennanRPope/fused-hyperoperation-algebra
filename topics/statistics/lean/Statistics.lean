/-
  Statistics in Fused Hyperoperation Algebra
  ==========================================

  Companion page: ../index.html.  Numbers: ../checks.py.

  Core Lean 4.37 only (no mathlib), no `sorry`, standard axioms only.
  Library-ready: imports FusedAlgebra.Base, .Cube and .Topics.E, and uses
  their mult-add a ◯ b = a + b + ab, weighted mult-add a ◯[k] b = a + b + k·ab
  (`circK k a b`), powers, list folds and the cube's corner expansion.

  What this file proves, by part:

  1. Indicators.  On 0/1 numbers the corner a·b is AND, k = −1 is OR, k = −2 is
     XOR, and (x − y)² = XOR.  The complement x ↦ 1 − x turns OR into AND (De
     Morgan).  Inclusion–exclusion for n events is the k = −1 fold, i.e. the
     signed corners of the n-cube.
  2. Finite probability spaces.  For events that need not be independent,
     P(A ∪ B) = P(A) ◯[k] P(B) with k = −lift, lift = P(A∩B)/(P(A)P(B)):
     disjoint is k = 0 (Kolmogorov's additivity), independent is k = −1.
     Inclusion–exclusion for n events with any dependence.
  3. The generating-function dial.  G(t) = E[t^X] multiplies under independent
     sums, so with t = 1 + k the "k-moment" a, 1 + k·a = G(1 + k), composes by
     ◯[k].  A coin's k-moment is p at every k.  t = 0 gives P(none) (the union
     rule), t = −1 gives P(even) − P(odd) (piling-up, unequal coins).  The
     binomial is the weighted Pascal layers C(n,j)·pʲ·qⁿ⁻ʲ.
  4. Moments and cumulants.  Moments of X + Y follow the Pascal rows
     1·2·1, 1·3·3·1, 1·4·6·4·1; cumulants κ₂, κ₃, κ₄ add.  For n copies the
     shape numbers fall exactly as 1/n (the exact part of the central limit
     theorem).  Coin and binomial mean, variance, κ₃, κ₄.
  5. Bayes.  Posterior odds = prior odds × likelihood ratio, so odds − 1
     mult-adds; the complement is the mult-add inverse in odds; pooling two
     opinions is velocity addition in ρ = 2p − 1.
  6. The dependence line.  Union under Ali–Mikhail–Haq dependence θ is the
     landscape's lens law at k = −(1 + θ), q = θ: the lens laws with one wall
     at certainty.  A lens law lets certainty absorb (a ⊙ 1 = 1) exactly when
     k + q = −1, so this line is every lens law that can be an "or".  θ = 0 is independence, θ = −1 velocity addition, θ = 1 the
     cusp.

  Not here (needs real numbers or limits): the normal density, the limit in
  the CLT, the Poisson law, logs.  Those are checked in ../checks.py.
-/
import FusedAlgebra.Base
import FusedAlgebra.Cube
import FusedAlgebra.Topics.E

open CircleOp

namespace StatFused

set_option linter.unusedSectionVars false

/-! # 1. Indicators: Boolean logic is weighted mult-add on {0, 1} -/

/-- An indicator: a 0/1 number, x·x = x. -/
def IsInd (x : Rat) : Prop := x * x = x

theorem ind_cases {x : Rat} (h : IsInd x) : x = 0 ∨ x = 1 := by
  unfold IsInd at h
  have h2 : x * (x - 1) = 0 := by grind
  rcases Rat.mul_eq_zero.mp h2 with h3 | h3
  · exact Or.inl h3
  · exact Or.inr (by grind)

theorem and_ind {x y : Rat} (hx : IsInd x) (hy : IsInd y) : IsInd (x * y) := by
  unfold IsInd at *; grind

theorem or_ind {x y : Rat} (hx : IsInd x) (hy : IsInd y) : IsInd (circK (-1) x y) := by
  rcases ind_cases hx with rfl | rfl <;> rcases ind_cases hy with rfl | rfl <;>
    (unfold IsInd circK; grind)

theorem xor_ind {x y : Rat} (hx : IsInd x) (hy : IsInd y) : IsInd (circK (-2) x y) := by
  rcases ind_cases hx with rfl | rfl <;> rcases ind_cases hy with rfl | rfl <;>
    (unfold IsInd circK; grind)

theorem not_ind {x : Rat} (hx : IsInd x) : IsInd (1 - x) := by
  unfold IsInd at *; grind

/-- Truth tables. -/
theorem or_table : circK (-1) (0:Rat) 0 = 0 ∧ circK (-1) (0:Rat) 1 = 1 ∧
    circK (-1) (1:Rat) 0 = 1 ∧ circK (-1) (1:Rat) 1 = 1 := by
  unfold circK; refine ⟨?_, ?_, ?_, ?_⟩ <;> grind

theorem xor_table : circK (-2) (0:Rat) 0 = 0 ∧ circK (-2) (0:Rat) 1 = 1 ∧
    circK (-2) (1:Rat) 0 = 1 ∧ circK (-2) (1:Rat) 1 = 0 := by
  unfold circK; refine ⟨?_, ?_, ?_, ?_⟩ <;> grind

theorem count_table : circK 0 (1:Rat) 1 = 2 := by unfold circK; grind

/-- The squared difference of two indicators is their XOR (k = −2). -/
theorem sq_diff_is_xor {x y : Rat} (hx : IsInd x) (hy : IsInd y) :
    (x - y) * (x - y) = circK (-2) x y := by
  unfold IsInd at *; unfold circK; grind

/-- De Morgan: "not (A or B)" = "not A and not B"; complement x ↦ 1 − x turns ◯[−1] into ×. -/
theorem de_morgan (x y : Rat) : 1 - circK (-1) x y = (1 - x) * (1 - y) := by
  unfold circK; grind

/-- The complement conjugates OR into AND: 1 − ((1−x) ◯[−1] (1−y)) = x·y. -/
theorem de_morgan_and (x y : Rat) : 1 - circK (-1) (1 - x) (1 - y) = x * y := by
  unfold circK; grind

/-- XOR is self-dual under complement of both inputs. -/
theorem xor_complement (x y : Rat) : circK (-2) (1 - x) (1 - y) = circK (-2) x y := by
  unfold circK; grind

/-! ### n events: the k = −1 fold and its signed corners -/

/-- Union of a list: x₁ ◯[−1] x₂ ◯[−1] … -/
def unionL : List Rat → Rat
  | [] => 0
  | a :: l => circK (-1) a (unionL l)

/-- Nothing happens: 1 − union = product of the misses. -/
theorem union_misses (l : List Rat) : 1 - unionL l = prodL (l.map (fun x => 1 - x)) := by
  induction l with
  | nil => simp only [unionL, List.map, prodL]; grind
  | cons a l ih =>
    simp only [unionL, List.map, prodL]; rw [← ih]; unfold circK; grind

/-- The k = −1 fold is mult-add read on the negatives. -/
theorem union_neg_foldCirc (l : List Rat) : unionL l = -(foldCirc (l.map (fun x => -x))) := by
  induction l with
  | nil => simp only [unionL, List.map, foldCirc]; grind
  | cons a l ih =>
    simp only [unionL, List.map, foldCirc]; rw [ih]; unfold circK circ; grind

/-- **Inclusion–exclusion for n events**: 1 − P(union) is the sum over every corner g
    of the n-cube of the signed product Π_{i∈g}(−xᵢ). -/
theorem inclusion_exclusion (l : List Rat) :
    1 - unionL l = sumL ((groups (l.map (fun x => -x))).map prodL) := by
  rw [union_neg_foldCirc, ← one_add_foldCirc_eq_corners]; grind

theorem groups_map {α β : Type} (f : α → β) (l : List α) :
    groups (l.map f) = (groups l).map (List.map f) := by
  induction l with
  | nil => rfl
  | cons a l ih =>
    simp only [List.map, groups, ih, List.map_append, List.map_map]; rfl

theorem incl_excl_two (p q : Rat) : unionL [p, q] = p + q - p * q := by
  simp only [unionL]; unfold circK; grind


/-! # 2. Finite probability spaces: dependence is the weight -/

section Spaces
variable {α : Type}

/-- A finite probability space: outcomes with weights. Expectation of f. -/
def Ex {α : Type} : List (α × Rat) → (α → Rat) → Rat
  | [], _ => 0
  | (a, w) :: ws, f => w * f a + Ex ws f


theorem Ex_add (ws : List (α × Rat)) (f g : α → Rat) :
    Ex ws (fun a => f a + g a) = Ex ws f + Ex ws g := by
  induction ws with
  | nil => simp only [Ex]; grind
  | cons p ws ih => obtain ⟨a, w⟩ := p; simp only [Ex]; rw [ih]; grind

theorem Ex_sub (ws : List (α × Rat)) (f g : α → Rat) :
    Ex ws (fun a => f a - g a) = Ex ws f - Ex ws g := by
  induction ws with
  | nil => simp only [Ex]; grind
  | cons p ws ih => obtain ⟨a, w⟩ := p; simp only [Ex]; rw [ih]; grind

theorem Ex_smul (ws : List (α × Rat)) (c : Rat) (f : α → Rat) :
    Ex ws (fun a => c * f a) = c * Ex ws f := by
  induction ws with
  | nil => simp only [Ex]; grind
  | cons p ws ih => obtain ⟨a, w⟩ := p; simp only [Ex]; rw [ih]; grind

theorem Ex_const (ws : List (α × Rat)) (c : Rat) :
    Ex ws (fun _ => c) = c * Ex ws (fun _ => 1) := by
  induction ws with
  | nil => simp only [Ex]; grind
  | cons p ws ih => obtain ⟨a, w⟩ := p; simp only [Ex]; rw [ih]; grind

theorem Ex_congr (ws : List (α × Rat)) (f g : α → Rat) (h : ∀ a, f a = g a) :
    Ex ws f = Ex ws g := by
  have : f = g := funext h
  rw [this]

/-- Linearity over a finite list of functions. -/
theorem Ex_sumL {β : Type} (ws : List (α × Rat)) (gs : List β) (F : β → α → Rat) :
    Ex ws (fun a => sumL (gs.map (fun g => F g a))) = sumL (gs.map (fun g => Ex ws (F g))) := by
  induction gs with
  | nil =>
    simp only [List.map, sumL]
    rw [Ex_const]; grind
  | cons g gs ih =>
    simp only [List.map, sumL]
    rw [Ex_add, ih]


/-! ## Events that need not be independent -/

/-- Union of two events, dependent or not: P(A ∪ B) = P(A) + P(B) − P(A ∩ B). -/
theorem union_general (ws : List (α × Rat)) (A B : α → Rat) :
    Ex ws (fun a => circK (-1) (A a) (B a)) = Ex ws A + Ex ws B - Ex ws (fun a => A a * B a) := by
  have h : ∀ a, circK (-1) (A a) (B a) = (A a + B a) - A a * B a := by
    intro a; unfold circK; grind
  rw [Ex_congr ws _ _ h, Ex_sub, Ex_add]

/-- **Dependence is the weight.** If the overlap is `lift` times what independence
    predicts, P(A ∩ B) = lift·P(A)·P(B), then P(A ∪ B) = P(A) ◯[k] P(B) with k = −lift. -/
theorem union_lift (ws : List (α × Rat)) (A B : α → Rat) (lift : Rat)
    (h : Ex ws (fun a => A a * B a) = lift * Ex ws A * Ex ws B) :
    Ex ws (fun a => circK (-1) (A a) (B a)) = circK (-lift) (Ex ws A) (Ex ws B) := by
  rw [union_general, h]; unfold circK; grind

/-- Disjoint events (lift 0): k = 0, plain addition, Kolmogorov's axiom. -/
theorem union_disjoint (ws : List (α × Rat)) (A B : α → Rat)
    (h : Ex ws (fun a => A a * B a) = 0) :
    Ex ws (fun a => circK (-1) (A a) (B a)) = circK 0 (Ex ws A) (Ex ws B) := by
  rw [union_general, h]; unfold circK; grind

/-- Exclusive or with dependence: k = −2·lift. -/
theorem xor_lift (ws : List (α × Rat)) (A B : α → Rat) (lift : Rat)
    (h : Ex ws (fun a => A a * B a) = lift * Ex ws A * Ex ws B) :
    Ex ws (fun a => circK (-2) (A a) (B a)) = circK (-2 * lift) (Ex ws A) (Ex ws B) := by
  have h' : ∀ a, circK (-2) (A a) (B a) = (A a + B a) - 2 * (A a * B a) := by
    intro a; unfold circK; grind
  rw [Ex_congr ws _ _ h', Ex_sub, Ex_add, Ex_smul, h]; unfold circK; grind

/-- Covariance of two indicators is the corner's excess over the independent corner. -/
theorem cov_is_corner_excess (ws : List (α × Rat)) (A B : α → Rat) (lift : Rat)
    (h : Ex ws (fun a => A a * B a) = lift * Ex ws A * Ex ws B) :
    Ex ws (fun a => A a * B a) - Ex ws A * Ex ws B = (lift - 1) * Ex ws A * Ex ws B := by
  rw [h]; grind

/-! ### Inclusion–exclusion for n events, with any dependence -/

def sgn : List β → Rat
  | [] => 1
  | _ :: l => -(sgn l)

theorem prodL_neg_sgn (g : List Rat) : prodL (g.map (fun x => -x)) = sgn g * prodL g := by
  induction g with
  | nil => simp only [List.map, prodL, sgn]; grind
  | cons a g ih => simp only [List.map, prodL, sgn]; rw [ih]; grind

theorem union_corners (l : List Rat) :
    1 - unionL l = sumL ((groups l).map (fun g => sgn g * prodL g)) := by
  have key : 1 - unionL l = 1 + foldCirc (l.map (fun x => -x)) := by
    induction l with
    | nil => simp only [unionL, List.map, foldCirc]; grind
    | cons a l ih =>
      simp only [unionL, List.map, foldCirc]; unfold circK circ; grind
  rw [key, one_add_foldCirc_eq_corners]
  have gm : ∀ (m : List Rat), groups (m.map (fun x => -x)) = (groups m).map (List.map (fun x => -x)) := by
    intro m
    induction m with
    | nil => rfl
    | cons a m ih => simp only [List.map, groups, ih, List.map_append, List.map_map]; rfl
  rw [gm, List.map_map]
  congr 1
  apply List.map_congr_left
  intro g _
  exact prodL_neg_sgn g

/-- The indicator of "every event in g happens". -/
def allOf (g : List (α → Rat)) (a : α) : Rat := prodL (g.map (fun A => A a))

theorem sgn_map {γ δ : Type} (f : γ → δ) (g : List γ) : sgn (g.map f) = sgn g := by
  induction g with
  | nil => rfl
  | cons a g ih => simp only [List.map, sgn]; rw [ih]

/-- **Inclusion–exclusion, n events, any dependence.** With total weight 1,
    P(no event happens) = Σ over every corner g of the cube of (−1)^|g|·P(all of g). -/
theorem inclusion_exclusion_general (ws : List (α × Rat)) (As : List (α → Rat)) :
    Ex ws (fun a => 1 - unionL (As.map (fun A => A a))) =
      sumL ((groups As).map (fun g => sgn g * Ex ws (allOf g))) := by
  have pt : ∀ a, 1 - unionL (As.map (fun A => A a)) =
      sumL ((groups As).map (fun g => sgn g * allOf g a)) := by
    intro a
    rw [union_corners, groups_map, List.map_map]
    congr 1
    apply List.map_congr_left
    intro g _
    simp only [Function.comp, allOf]; rw [sgn_map]
  rw [Ex_congr ws _ _ pt, Ex_sumL]
  congr 1
  apply List.map_congr_left
  intro g _
  exact Ex_smul ws (sgn g) (allOf g)

end Spaces

/-! # 3. The generating-function dial t = 1 + k -/

/-! ## Counts and the generating-function dial t = 1 + k -/

abbrev NDist := List (Nat × Rat)

/-- G(t) = E[t^X]. -/
def gen : NDist → Rat → Rat
  | [], _ => 0
  | (x, w) :: d, t => w * pw t x + gen d t

def nshift (x : Nat) (w : Rat) (e : NDist) : NDist := e.map (fun p => (x + p.1, w * p.2))

def nconv : NDist → NDist → NDist
  | [], _ => []
  | (x, w) :: d, e => nshift x w e ++ nconv d e

theorem gen_append (d e : NDist) (t : Rat) : gen (d ++ e) t = gen d t + gen e t := by
  induction d with
  | nil => simp only [List.nil_append, gen]; grind
  | cons p d ih => obtain ⟨x, w⟩ := p; simp only [List.cons_append, gen]; rw [ih]; grind

theorem gen_nshift (x : Nat) (w t : Rat) (e : NDist) : gen (nshift x w e) t = w * pw t x * gen e t := by
  induction e with
  | nil => simp only [nshift, List.map, gen]; grind
  | cons p e ih =>
    obtain ⟨y, v⟩ := p; simp only [nshift, List.map, gen] at ih ⊢; rw [ih, pw_add]; grind

/-- Generating functions multiply under independent sums: the volume reading (1+a)(1+b). -/
theorem gen_conv (d e : NDist) (t : Rat) : gen (nconv d e) t = gen d t * gen e t := by
  induction d with
  | nil => simp only [nconv, gen]; grind
  | cons p d ih =>
    obtain ⟨x, w⟩ := p; simp only [nconv, gen]; rw [gen_append, gen_nshift, ih]; grind

/-- **The dial.** If a and b are the k-moments of X and Y, read through 1 + k·a = G(1 + k),
    then the k-moment of X + Y is a ◯[k] b. -/
theorem kmoment_law (k a b : Rat) (d e : NDist) (ha : 1 + k * a = gen d (1 + k))
    (hb : 1 + k * b = gen e (1 + k)) : 1 + k * circK k a b = gen (nconv d e) (1 + k) := by
  rw [gen_conv, ← ha, ← hb, circK_lift]

def nbern (p : Rat) : NDist := [(0, 1 - p), (1, p)]

/-- A coin's k-moment is p at every k: one number, every weight. -/
theorem bern_kmoment (k p : Rat) : gen (nbern p) (1 + k) = 1 + k * p := by
  simp only [nbern, gen, pw]; grind

/-- Coins with different chances: X = number of heads. -/
def heads : List Rat → NDist
  | [] => [(0, 1)]
  | p :: ps => nconv (nbern p) (heads ps)

def foldW (k : Rat) : List Rat → Rat
  | [] => 0
  | p :: ps => circK k p (foldW k ps)

theorem heads_kmoment (k : Rat) (ps : List Rat) : gen (heads ps) (1 + k) = 1 + k * foldW k ps := by
  induction ps with
  | nil => simp only [heads, foldW, gen, pw]; grind
  | cons p ps ih => simp only [heads, foldW]; rw [gen_conv, bern_kmoment, ih, circK_lift]

/-- t = 0 (k = −1): G(0) = P(no heads) = Π(1 − pᵢ), the union rule's survival. -/
theorem heads_none (ps : List Rat) : gen (heads ps) 0 = 1 - foldW (-1) ps := by
  have h := heads_kmoment (-1) ps
  have e : (1 : Rat) + -1 = 0 := by grind
  rw [e] at h; rw [h]; grind

/-- t = −1 (k = −2): G(−1) = P(even) − P(odd) = 1 − 2·(p₁ ◯[−2] p₂ ◯[−2] …): piling-up for unequal coins. -/
theorem heads_parity (ps : List Rat) : gen (heads ps) (-1) = 1 - 2 * foldW (-2) ps := by
  have h := heads_kmoment (-2) ps
  have e : (1 : Rat) + -2 = -1 := by grind
  rw [e] at h; rw [h]; grind

theorem pw_neg_one_parity (x : Nat) : pw (-1 : Rat) x = (if x % 2 = 0 then 1 else -1) := by
  induction x with
  | zero => rfl
  | succ n ih =>
    rw [pw, ih]
    by_cases h : n % 2 = 0
    · have h' : ¬ (n + 1) % 2 = 0 := by omega
      simp only [h, h', ↓reduceIte]; grind
    · have h' : (n + 1) % 2 = 0 := by omega
      simp only [h, h', ↓reduceIte]; grind

def pEven : NDist → Rat
  | [] => 0
  | (x, w) :: d => (if x % 2 = 0 then w else 0) + pEven d
def pOdd : NDist → Rat
  | [] => 0
  | (x, w) :: d => (if x % 2 = 0 then 0 else w) + pOdd d

theorem gen_neg_one (d : NDist) : gen d (-1) = pEven d - pOdd d := by
  induction d with
  | nil => simp only [gen, pEven, pOdd]; grind
  | cons p d ih =>
    obtain ⟨x, w⟩ := p
    simp only [gen, pEven, pOdd]; rw [ih, pw_neg_one_parity]
    by_cases h : x % 2 = 0
    · simp only [h, ↓reduceIte]; grind
    · simp only [h, ↓reduceIte]; grind

def pZero : NDist → Rat
  | [] => 0
  | (x, w) :: d => (if x = 0 then w else 0) + pZero d

theorem gen_zero (d : NDist) : gen d 0 = pZero d := by
  induction d with
  | nil => rfl
  | cons p d ih =>
    obtain ⟨x, w⟩ := p
    simp only [gen, pZero]; rw [ih]
    cases x with
    | zero => simp [pw]
    | succ n => simp [pw]

/-- n equal coins: the k-moment climbs, 1 + k·a = (1 + kp)ⁿ = 1 + k·p^(◯[k] n) (rung 2). -/
def binom (p : Rat) : Nat → NDist
  | 0 => [(0, 1)]
  | n + 1 => nconv (nbern p) (binom p n)

theorem binom_kmoment (k p : Rat) (n : Nat) : gen (binom p n) (1 + k) = 1 + k * wpow k p n := by
  induction n with
  | zero => simp only [binom, wpow, gen, pw]; grind
  | succ n ih =>
    simp only [binom, wpow]; rw [gen_conv, bern_kmoment, ih, circK_lift]; grind

/-- The k-moment counts lit corners: (1 + k)^x = Σⱼ C(x, j)·kʲ, so
    G(1 + k) = Σ over outcomes of weight × Σⱼ C(x, j)·kʲ (factorial moments). -/
def litCorners : NDist → Rat → Rat
  | [], _ => 0
  | (x, w) :: d, k => w * E.rsum (fun j => (Cube.layer x j : Rat) * k ^ j) (x + 1) + litCorners d k

theorem gen_is_lit_corners (d : NDist) (k : Rat) : gen d (1 + k) = litCorners d k := by
  induction d with
  | nil => rfl
  | cons p d ih =>
    obtain ⟨x, w⟩ := p
    simp only [gen, litCorners]; rw [ih, pw_eq_pow, E.binomial]

/-! ### The binomial distribution is the weighted Pascal layers -/

/-- P(n coins show j heads), built by adding one coin at a time. -/
def bp (p q : Rat) : Nat → Nat → Rat
  | 0, 0 => 1
  | 0, _ + 1 => 0
  | n + 1, 0 => q * bp p q n 0
  | n + 1, j + 1 => q * bp p q n (j + 1) + p * bp p q n j

theorem bp_zero_col (p q : Rat) (n : Nat) : bp p q n 0 = pw q n := by
  induction n with
  | zero => rfl
  | succ n ih => simp only [bp, pw]; rw [ih]; grind

theorem bp_above (p q : Rat) : ∀ n j, n < j → bp p q n j = 0
  | 0, 0, h => absurd h (by omega)
  | 0, _ + 1, _ => rfl
  | _ + 1, 0, h => absurd h (by omega)
  | n + 1, j + 1, h => by
    simp only [bp]; rw [bp_above p q n (j + 1) (by omega), bp_above p q n j (by omega)]; grind

/-- **Binomial = Pascal layers**: P(j heads of n) = C(n, j)·pʲ·qⁿ⁻ʲ. -/
theorem bp_layers (p q : Rat) : ∀ n j, j ≤ n →
    bp p q n j = (Cube.layer n j : Rat) * pw p j * pw q (n - j)
  | n, 0, _ => by rw [bp_zero_col, Cube.layer_zero_right]; simp only [pw, Nat.sub_zero]; grind
  | 0, _ + 1, h => absurd h (by omega)
  | n + 1, j + 1, h => by
    simp only [bp]
    rw [bp_layers p q n j (by omega), Cube.layer_succ_succ, Rat.natCast_add]
    by_cases hj : j + 1 ≤ n
    · rw [bp_layers p q n (j + 1) hj]
      have e1 : n + 1 - (j + 1) = n - (j + 1) + 1 := by omega
      have e2 : n - j = n - (j + 1) + 1 := by omega
      rw [e1, e2]; simp only [pw]; grind
    · have hn : j = n := by omega
      subst hn
      rw [bp_above p q j (j + 1) (by omega), Cube.layer_eq_zero j (j + 1) (by omega)]
      simp only [Nat.sub_self, pw]; grind


/-! # 4. Moments and cumulants: Pascal rows, then the log coordinate -/

/-! Distributions: lists of (value, weight). Independent sum = convolution. -/

abbrev Dist := List (Rat × Rat)

def mass : Dist → Rat
  | [] => 0
  | (_, w) :: d => w + mass d
def m1 : Dist → Rat
  | [] => 0
  | (x, w) :: d => w * x + m1 d
def m2 : Dist → Rat
  | [] => 0
  | (x, w) :: d => w * (x * x) + m2 d
def m3 : Dist → Rat
  | [] => 0
  | (x, w) :: d => w * (x * x * x) + m3 d
def m4 : Dist → Rat
  | [] => 0
  | (x, w) :: d => w * (x * x * x * x) + m4 d

/-- Shift every value by x and scale every weight by w. -/
def shift (x w : Rat) (e : Dist) : Dist := e.map (fun p => (x + p.1, w * p.2))

/-- The distribution of X + Y for independent X, Y. -/
def conv : Dist → Dist → Dist
  | [], _ => []
  | (x, w) :: d, e => shift x w e ++ conv d e

theorem mass_append (d e : Dist) : mass (d ++ e) = mass d + mass e := by
  induction d with
  | nil => simp only [List.nil_append, mass]; grind
  | cons p d ih => obtain ⟨x, w⟩ := p; simp only [List.cons_append, mass]; rw [ih]; grind
theorem m1_append (d e : Dist) : m1 (d ++ e) = m1 d + m1 e := by
  induction d with
  | nil => simp only [List.nil_append, m1]; grind
  | cons p d ih => obtain ⟨x, w⟩ := p; simp only [List.cons_append, m1]; rw [ih]; grind
theorem m2_append (d e : Dist) : m2 (d ++ e) = m2 d + m2 e := by
  induction d with
  | nil => simp only [List.nil_append, m2]; grind
  | cons p d ih => obtain ⟨x, w⟩ := p; simp only [List.cons_append, m2]; rw [ih]; grind
theorem m3_append (d e : Dist) : m3 (d ++ e) = m3 d + m3 e := by
  induction d with
  | nil => simp only [List.nil_append, m3]; grind
  | cons p d ih => obtain ⟨x, w⟩ := p; simp only [List.cons_append, m3]; rw [ih]; grind
theorem m4_append (d e : Dist) : m4 (d ++ e) = m4 d + m4 e := by
  induction d with
  | nil => simp only [List.nil_append, m4]; grind
  | cons p d ih => obtain ⟨x, w⟩ := p; simp only [List.cons_append, m4]; rw [ih]; grind

theorem mass_shift (x w : Rat) (e : Dist) : mass (shift x w e) = w * mass e := by
  induction e with
  | nil => simp only [shift, List.map, mass]; grind
  | cons p e ih =>
    obtain ⟨y, v⟩ := p; simp only [shift, List.map, mass] at ih ⊢; rw [ih]; grind
theorem m1_shift (x w : Rat) (e : Dist) : m1 (shift x w e) = w * (x * mass e + m1 e) := by
  induction e with
  | nil => simp only [shift, List.map, m1, mass]; grind
  | cons p e ih =>
    obtain ⟨y, v⟩ := p; simp only [shift, List.map, m1, mass] at ih ⊢; rw [ih]; grind
theorem m2_shift (x w : Rat) (e : Dist) :
    m2 (shift x w e) = w * (x * x * mass e + 2 * x * m1 e + m2 e) := by
  induction e with
  | nil => simp only [shift, List.map, m2, m1, mass]; grind
  | cons p e ih =>
    obtain ⟨y, v⟩ := p; simp only [shift, List.map, m2, m1, mass] at ih ⊢; rw [ih]; grind
theorem m3_shift (x w : Rat) (e : Dist) :
    m3 (shift x w e) = w * (x*x*x * mass e + 3 * (x*x) * m1 e + 3 * x * m2 e + m3 e) := by
  induction e with
  | nil => simp only [shift, List.map, m3, m2, m1, mass]; grind
  | cons p e ih =>
    obtain ⟨y, v⟩ := p; simp only [shift, List.map, m3, m2, m1, mass] at ih ⊢; rw [ih]; grind
theorem m4_shift (x w : Rat) (e : Dist) :
    m4 (shift x w e) = w * (x*x*x*x * mass e + 4 * (x*x*x) * m1 e + 6 * (x*x) * m2 e
      + 4 * x * m3 e + m4 e) := by
  induction e with
  | nil => simp only [shift, List.map, m4, m3, m2, m1, mass]; grind
  | cons p e ih =>
    obtain ⟨y, v⟩ := p; simp only [shift, List.map, m4, m3, m2, m1, mass] at ih ⊢; rw [ih]; grind

/-- Total weight multiplies: the grid of outcomes. -/
theorem mass_conv (d e : Dist) : mass (conv d e) = mass d * mass e := by
  induction d with
  | nil => simp only [conv, mass]; grind
  | cons p d ih => obtain ⟨x, w⟩ := p; simp only [conv, mass]; rw [mass_append, mass_shift, ih]; grind

/-- Moments of a sum: the Pascal rows 1·1, 1·2·1, 1·3·3·1, 1·4·6·4·1. -/
theorem m1_conv (d e : Dist) : m1 (conv d e) = m1 d * mass e + mass d * m1 e := by
  induction d with
  | nil => simp only [conv, m1, mass]; grind
  | cons p d ih =>
    obtain ⟨x, w⟩ := p; simp only [conv, m1, mass]; rw [m1_append, m1_shift, ih]; grind
theorem m2_conv (d e : Dist) :
    m2 (conv d e) = m2 d * mass e + 2 * m1 d * m1 e + mass d * m2 e := by
  induction d with
  | nil => simp only [conv, m2, m1, mass]; grind
  | cons p d ih =>
    obtain ⟨x, w⟩ := p; simp only [conv, m2, m1, mass]; rw [m2_append, m2_shift, ih]; grind
theorem m3_conv (d e : Dist) :
    m3 (conv d e) = m3 d * mass e + 3 * m2 d * m1 e + 3 * m1 d * m2 e + mass d * m3 e := by
  induction d with
  | nil => simp only [conv, m3, m2, m1, mass]; grind
  | cons p d ih =>
    obtain ⟨x, w⟩ := p; simp only [conv, m3, m2, m1, mass]; rw [m3_append, m3_shift, ih]; grind
theorem m4_conv (d e : Dist) :
    m4 (conv d e) = m4 d * mass e + 4 * m3 d * m1 e + 6 * m2 d * m2 e
      + 4 * m1 d * m3 e + mass d * m4 e := by
  induction d with
  | nil => simp only [conv, m4, m3, m2, m1, mass]; grind
  | cons p d ih =>
    obtain ⟨x, w⟩ := p; simp only [conv, m4, m3, m2, m1, mass]; rw [m4_append, m4_shift, ih]; grind

/-! Cumulants: the additive coordinates. -/
def k2 (d : Dist) : Rat := m2 d - m1 d * m1 d
def k3 (d : Dist) : Rat := m3 d - 3 * m2 d * m1 d + 2 * (m1 d * m1 d * m1 d)
def k4 (d : Dist) : Rat :=
  m4 d - 4 * m3 d * m1 d - 3 * (m2 d * m2 d) + 12 * m2 d * (m1 d * m1 d)
    - 6 * (m1 d * m1 d * m1 d * m1 d)

theorem mean_add (d e : Dist) (hd : mass d = 1) (he : mass e = 1) :
    m1 (conv d e) = m1 d + m1 e := by rw [m1_conv, hd, he]; grind
theorem var_add (d e : Dist) (hd : mass d = 1) (he : mass e = 1) :
    k2 (conv d e) = k2 d + k2 e := by
  unfold k2; rw [m2_conv, m1_conv, hd, he]; grind
theorem k3_add (d e : Dist) (hd : mass d = 1) (he : mass e = 1) :
    k3 (conv d e) = k3 d + k3 e := by
  unfold k3; rw [m3_conv, m2_conv, m1_conv, hd, he]; grind
theorem k4_add (d e : Dist) (hd : mass d = 1) (he : mass e = 1) :
    k4 (conv d e) = k4 d + k4 e := by
  unfold k4; rw [m4_conv, m3_conv, m2_conv, m1_conv, hd, he]; grind


/-! ## n independent copies and the CLT's shrinking shape terms -/

def point0 : Dist := [(0, 1)]

def convPow (d : Dist) : Nat → Dist
  | 0 => point0
  | n + 1 => conv (convPow d n) d

theorem convPow_mass (d : Dist) (hd : mass d = 1) (n : Nat) : mass (convPow d n) = 1 := by
  induction n with
  | zero => simp only [convPow, point0, mass]; grind
  | succ n ih => simp only [convPow]; rw [mass_conv, ih, hd]; grind

theorem convPow_mean (d : Dist) (hd : mass d = 1) (n : Nat) : m1 (convPow d n) = n * m1 d := by
  induction n with
  | zero => simp only [convPow, point0, m1]; grind
  | succ n ih =>
    simp only [convPow]; rw [mean_add _ _ (convPow_mass d hd n) hd, ih, Rat.natCast_add]; grind

theorem convPow_k2 (d : Dist) (hd : mass d = 1) (n : Nat) : k2 (convPow d n) = n * k2 d := by
  induction n with
  | zero => simp only [convPow, point0, k2, m2, m1]; grind
  | succ n ih =>
    simp only [convPow]; rw [var_add _ _ (convPow_mass d hd n) hd, ih, Rat.natCast_add]; grind

theorem convPow_k3 (d : Dist) (hd : mass d = 1) (n : Nat) : k3 (convPow d n) = n * k3 d := by
  induction n with
  | zero => simp only [convPow, point0, k3, m3, m2, m1]; grind
  | succ n ih =>
    simp only [convPow]; rw [k3_add _ _ (convPow_mass d hd n) hd, ih, Rat.natCast_add]; grind

theorem convPow_k4 (d : Dist) (hd : mass d = 1) (n : Nat) : k4 (convPow d n) = n * k4 d := by
  induction n with
  | zero => simp only [convPow, point0, k4, m4, m3, m2, m1]; grind
  | succ n ih =>
    simp only [convPow]; rw [k4_add _ _ (convPow_mass d hd n) hd, ih, Rat.natCast_add]; grind

/-- Rescaling X ↦ c·X. -/
def scale (c : Rat) (d : Dist) : Dist := d.map (fun p => (c * p.1, p.2))

theorem scale_m (c : Rat) (d : Dist) : mass (scale c d) = mass d ∧ m1 (scale c d) = c * m1 d ∧
    m2 (scale c d) = c*c * m2 d ∧ m3 (scale c d) = c*c*c * m3 d ∧
    m4 (scale c d) = c*c*c*c * m4 d := by
  induction d with
  | nil => simp only [scale, List.map, mass, m1, m2, m3, m4]; grind
  | cons p d ih =>
    obtain ⟨x, w⟩ := p
    simp only [scale, List.map, mass, m1, m2, m3, m4] at ih ⊢
    obtain ⟨h0, h1, h2, h3, h4⟩ := ih
    refine ⟨?_, ?_, ?_, ?_, ?_⟩ <;> grind

/-- Cumulant r scales as cʳ: the reason the shape terms die after standardizing. -/
theorem scale_k (c : Rat) (d : Dist) : k2 (scale c d) = c^2 * k2 d ∧
    k3 (scale c d) = c^3 * k3 d ∧ k4 (scale c d) = c^4 * k4 d := by
  obtain ⟨_, h1, h2, h3, h4⟩ := scale_m c d
  unfold k2 k3 k4; rw [h1, h2, h3, h4]
  refine ⟨?_, ?_, ?_⟩ <;> grind

/-- Skewness² and excess kurtosis, the two shape numbers (scale-free). -/
def skewSq (d : Dist) : Rat := k3 d * k3 d / (k2 d * k2 d * k2 d)
def kurt (d : Dist) : Rat := k4 d / (k2 d * k2 d)

theorem skewSq_scale (c : Rat) (hc : c ≠ 0) (d : Dist) : skewSq (scale c d) = skewSq d := by
  obtain ⟨h2, h3, _⟩ := scale_k c d
  unfold skewSq; rw [h2, h3]
  by_cases hk : k2 d = 0
  · rw [hk]; grind
  · have : c ^ 2 ≠ 0 := by grind
    grind

theorem kurt_scale (c : Rat) (hc : c ≠ 0) (d : Dist) : kurt (scale c d) = kurt d := by
  obtain ⟨h2, _, h4⟩ := scale_k c d
  unfold kurt; rw [h2, h4]
  by_cases hk : k2 d = 0
  · rw [hk]; grind
  · grind

/-- **CLT, the exact part.** For a sum of n independent copies the shape numbers
    fall as 1/n: skewness²(Sₙ) = skewness²(X)/n and kurtosis(Sₙ) = kurtosis(X)/n.
    Only the mean and variance survive, which is the normal law. -/
theorem clt_skew (d : Dist) (hd : mass d = 1) (n : Nat) (hn : 0 < n) (hk : k2 d ≠ 0) :
    skewSq (convPow d n) * n = skewSq d := by
  unfold skewSq; rw [convPow_k2 d hd n, convPow_k3 d hd n]
  have : (n : Rat) ≠ 0 := by
    have := Rat.natCast_pos.mpr hn; grind
  grind

theorem clt_kurt (d : Dist) (hd : mass d = 1) (n : Nat) (hn : 0 < n) (hk : k2 d ≠ 0) :
    kurt (convPow d n) * n = kurt d := by
  unfold kurt; rw [convPow_k2 d hd n, convPow_k4 d hd n]
  have : (n : Rat) ≠ 0 := by
    have := Rat.natCast_pos.mpr hn; grind
  grind

/-! ## The coin -/

def bern (p : Rat) : Dist := [(0, 1 - p), (1, p)]

theorem bern_stats (p : Rat) : mass (bern p) = 1 ∧ m1 (bern p) = p ∧
    k2 (bern p) = p * (1 - p) ∧ k3 (bern p) = p * (1 - p) * (1 - 2 * p) ∧
    k4 (bern p) = p * (1 - p) * (1 - 6 * p + 6 * p * p) := by
  unfold k2 k3 k4; simp only [bern, mass, m1, m2, m3, m4]
  refine ⟨?_, ?_, ?_, ?_, ?_⟩ <;> grind

/-- A coin's variance is half its self-parity: Var = ½·(p ◯[−2] p), the chance two copies differ, halved. -/
theorem bern_var_parity (p : Rat) : k2 (bern p) = circK (-2) p p / 2 := by
  rw [(bern_stats p).2.2.1]; unfold circK; grind

/-- Binomial mean np and variance np(1 − p), from additivity alone. -/
theorem binom_mean_var (p : Rat) (n : Nat) :
    m1 (convPow (bern p) n) = n * p ∧ k2 (convPow (bern p) n) = n * (p * (1 - p)) := by
  have hm := (bern_stats p).1
  rw [convPow_mean _ hm, convPow_k2 _ hm, (bern_stats p).2.1, (bern_stats p).2.2.1]
  exact ⟨rfl, rfl⟩

/-- Fair coins: skewness 0 and kurtosis −2/n, falling to the normal's 0. -/
theorem fair_coin_kurt (n : Nat) (hn : 0 < n) : kurt (convPow (bern (1/2)) n) * n = -2 := by
  have h := clt_kurt (bern (1/2)) (bern_stats _).1 n hn (by rw [(bern_stats _).2.2.1]; grind)
  rw [h]; unfold kurt; rw [(bern_stats _).2.2.1, (bern_stats _).2.2.2.2]; grind


/-! # 5. Bayes: odds multiply, so odds − 1 mult-adds -/

def odds (p : Rat) : Rat := p / (1 - p)

/-- Bayes' rule in odds form: posterior odds = prior odds × likelihood ratio. -/
theorem bayes_odds (p lH lN : Rat) (hp0 : p ≠ 0) (hp1 : p ≠ 1) (hH : lH ≠ 0) (hN : lN ≠ 0)
    (hden : p * lH + (1 - p) * lN ≠ 0) :
    odds (p * lH / (p * lH + (1 - p) * lN)) = odds p * (lH / lN) := by
  unfold odds
  have h1 : (1 : Rat) - p ≠ 0 := by grind
  have e : 1 - p * lH / (p * lH + (1 - p) * lN) = (1 - p) * lN / (p * lH + (1 - p) * lN) := by
    grind
  rw [e]; grind

/-- Two independent pieces of evidence: o·L₁·L₂ − 1 = (o − 1) ◯ (L₁ − 1) ◯ (L₂ − 1). -/
theorem evidence_circ (o l1 l2 : Rat) :
    o * l1 * l2 - 1 = circ (circ (o - 1) (l1 - 1)) (l2 - 1) := by
  unfold circ; grind

/-- The same evidence n times: o·Lⁿ − 1 = (o − 1) ◯ (L − 1)^(◯n), the climb (rung 2). -/
theorem repeated_evidence (o l : Rat) (n : Nat) :
    o * l ^ n - 1 = circ (o - 1) (mpow (l - 1) n) := by
  have h := mpow_lift (l - 1) n
  have e : (1 : Rat) + (l - 1) = l := by grind
  rw [e] at h; unfold circ; rw [← h]; grind

/-- The complement is the way back: odds(1 − p) − 1 is the mult-add inverse of odds(p) − 1. -/
theorem complement_is_inverse (p : Rat) (hp0 : p ≠ 0) (hp1 : p ≠ 1) :
    odds (1 - p) - 1 = circInv (odds p - 1) := by
  unfold odds circInv
  have h1 : (1 : Rat) - p ≠ 0 := by grind
  have e : (1 : Rat) - (1 - p) = p := by grind
  rw [e]; grind

theorem complement_cancels (p : Rat) (hp0 : p ≠ 0) (hp1 : p ≠ 1) :
    circ (odds p - 1) (odds (1 - p) - 1) = 0 := by
  rw [complement_is_inverse p hp0 hp1]
  apply circ_circInv
  unfold odds; have h1 : (1 : Rat) - p ≠ 0 := by grind
  grind

/-- Pooling two independent opinions from a flat prior. -/
def pool (p q : Rat) : Rat := p * q / (p * q + (1 - p) * (1 - q))

/-- In the coordinate ρ = 2p − 1 pooling is velocity addition, the lens law at
    landscape point (0, −1), walls at ±1 (certainty). -/
theorem pool_velocity (p q : Rat) (h : p * q + (1 - p) * (1 - q) ≠ 0) :
    2 * pool p q - 1 = ((2*p - 1) + (2*q - 1)) / (1 + (2*p - 1) * (2*q - 1)) := by
  unfold pool
  have e : 1 + (2*p - 1) * (2*q - 1) = 2 * (p * q + (1 - p) * (1 - q)) := by grind
  rw [e]; grind

/-- Pooling multiplies odds: the same law in the odds coordinate. -/
theorem pool_odds (p q : Rat) (hp : p ≠ 1) (hq : q ≠ 1) (h : p * q + (1 - p) * (1 - q) ≠ 0) :
    odds (pool p q) = odds p * odds q := by
  unfold odds pool
  have hp' : (1 : Rat) - p ≠ 0 := by grind
  have hq' : (1 : Rat) - q ≠ 0 := by grind
  have e : 1 - p * q / (p * q + (1 - p) * (1 - q)) =
      (1 - p) * (1 - q) / (p * q + (1 - p) * (1 - q)) := by grind
  rw [e]; grind


/-! # 6. The dependence line: lens laws with a wall at certainty -/

/-- The Ali–Mikhail–Haq copula: P(both) for marginals u, v and dependence θ ∈ [−1, 1]. -/
def amh (θ u v : Rat) : Rat := u * v / (1 - θ * (1 - u) * (1 - v))

/-- The landscape's lens law at weight k and bend q. -/
def lens (k q a b : Rat) : Rat := (a + b + k * a * b) / (1 - q * a * b)

/-- **Union under AMH dependence is a lens law** at k = −(1 + θ), q = θ:
    "neither happens" follows the AMH copula, so "at least one happens" is
    (p + q − (1 + θ)pq)/(1 − θpq). -/
theorem amh_union_is_lens (θ p q : Rat) (h : 1 - θ * p * q ≠ 0) :
    1 - amh θ (1 - p) (1 - q) = lens (-(1 + θ)) θ p q := by
  unfold amh lens
  have e : 1 - θ * (1 - (1 - p)) * (1 - (1 - q)) = 1 - θ * p * q := by grind
  rw [e]; grind

/-- Its walls: 1 + k·x + q·x² = (1 − x)(1 − θx). One wall sits at certainty, x = 1;
    the other at 1/θ. -/
theorem amh_walls (θ x : Rat) : 1 + (-(1 + θ)) * x + θ * (x * x) = (1 - x) * (1 - θ * x) := by
  grind

/-- The lens M(x) = (1 − x)/(1 − θx) turns the law into multiplication, so it is
    associative: two independent-looking "survivals" multiply in this coordinate. -/
theorem amh_lens_multiplies (θ a b : Rat) (h : 1 - θ * a * b ≠ 0) :
    (1 - lens (-(1 + θ)) θ a b) * ((1 - θ * a) * (1 - θ * b)) =
      (1 - a) * (1 - b) * (1 - θ * lens (-(1 + θ)) θ a b) := by
  unfold lens; grind

/-- θ = 0: independence, the union rule a ◯[−1] b. -/
theorem amh_independent (a b : Rat) : lens (-(1 + 0)) 0 a b = circK (-1) a b := by
  unfold lens circK; grind

/-- θ = −1: walls at ±1, velocity addition, the same law as pooling opinions. -/
theorem amh_velocity (a b : Rat) : lens (-(1 + -1)) (-1) a b = (a + b) / (1 + a * b) := by
  unfold lens; grind

/-- θ = 1: a double wall at certainty, the border parabola k² = 4q (the cusp). -/
theorem amh_cusp : (-(1 + (1 : Rat))) * (-(1 + 1)) = 4 * 1 := by grind

/-- Every θ lands on the line k + q = −1 of the landscape plane. -/
theorem amh_line (θ : Rat) : (-(1 + θ)) + θ = -1 := by grind

/-- **Why the line.** An "or" law must let certainty swallow everything: a ⊙ 1 = 1.
    A lens law does that exactly when k + q = −1, i.e. when it has a wall at 1. -/
theorem certainty_absorbs_iff (k q a : Rat) (ha : a ≠ 0) (hd : 1 - q * a * 1 ≠ 0) :
    lens k q a 1 = 1 ↔ k + q = -1 := by
  unfold lens
  constructor
  · intro h
    have h2 : a + 1 + k * a * 1 = 1 - q * a * 1 := by
      have e : (a + 1 + k * a * 1) / (1 - q * a * 1) * (1 - q * a * 1) = a + 1 + k * a * 1 := by
        grind
      rw [h] at e; grind
    have h3 : a * (1 + k + q) = 0 := by grind
    rcases Rat.mul_eq_zero.mp h3 with h4 | h4
    · exact absurd h4 ha
    · grind
  · intro h
    have e : a + 1 + k * a * 1 = 1 - q * a * 1 := by grind
    rw [e]; grind

/-- The AMH copula at θ = 0 is independence: P(both) = uv. -/
theorem amh_zero (u v : Rat) : amh 0 u v = u * v := by unfold amh; grind


end StatFused
