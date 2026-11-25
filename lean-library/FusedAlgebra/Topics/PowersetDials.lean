/-
  The power set with one dial per element
  =======================================

  Part of Fused Hyperoperation Algebra. Companion page and checks:
  math-theory/fractional-powerset-dials/ (index.html, checks.py).
  Core Lean 4.37 only (no mathlib), no `sorry`. Imports Base only.

  Every element x of a finite set gets a dial t_x. A subset A is weighed by
  the product of the dials of its members. What is proved, by part:

  1. The count. The weighed subsets add up to ∏(1 + t_x), which is
     1 + (t₁ ◯ t₂ ◯ … ◯ tₙ). All dials 1 give 2ⁿ, all 0 give 1, dials that are
     0 or 1 give 2^(number of ones): the old split P(A) × P(S∖A).
  2. Merging. Combining two dialled power sets of the same set by union
     (every pair (A, B), weighed by the product, landing on A ∪ B) gives the
     dialled power set whose dials are mult-added, element by element
     (`conv_W`). The same one-element computation (`cell_merge`) says that
     merging two elements into one block, present when either is present,
     gives the block the dial t ◯ s. Dial 0 is the unit of merging.
  3. The half. If r·r = 2 then h = r − 1 has h ◯ h = 1, so the dialled power
     set with every dial h, merged with itself, is exactly P(S) with every
     subset weighed 1 (`half_squared_family`). In general a dial t with
     t^(◯m) = 1 is an m-th root. Halving commutes with merging.
     Each half has a partner: t ◯ (1 − t)/(1 + t) = 1; partnering is an
     involution swapping 0 and 1 and fixing exactly the dials with t ◯ t = 1.
  4. Complement. Taking complements turns dial t into 1/t (`W_compl`), so it
     swaps "merge by union" (mult-add on the dials) with "merge by
     intersection" (mult-add on the reciprocal dials). The only dial that
     complement leaves alone is 1: the full power set, not a half. The union
     half and the intersection half are reciprocal (√2 − 1 and √2 + 1); in
     probabilities they are 1 − r/2 and r/2, which add to 1, and no single
     probability is a half both ways (`no_half_both_ways`).
  5. No dial gives a half step by iteration: whatever the dials, the empty set
     has exactly one subset of weight 1, and that alone rules out
     f (f n) = 2ⁿ (`no_iterate_half`).
  6. Known sets. Primes: the factor identity behind ∏(1 + p^(−s)) = ζ(s)/ζ(2s)
     (`euler_squarefree`); the naturals as merged squarefree layers
     (`fermi_dirac`); the rationals' factor (1 + x)/(1 − x) (`rational_factor`);
     ℕ an exact half of ℤ∖{0} (`integers_square`); tiny dials linearise, so the
     plain half and the merging half agree to second order (`tiny_half`);
     participation complement fixes 1/2 only (`participation_fixed`).
     The growth exponents (Selberg–Delange, Meinardus) are cited, not proved.
-/
import FusedAlgebra.Base

namespace PowersetDials

open CircleOp
open Lean.Grind

set_option linter.unusedSectionVars false

open scoped CircleOp

/-! ## 1. Subsets as bit vectors, weights, and the count -/

section Semiring
variable {R : Type} [Lean.Grind.CommSemiring R]

/-- Every bit vector of length n: the subsets of an n-element set (the corners
    of the cube [0,1]ⁿ). -/
def vecs : Nat → List (List Bool)
  | 0 => [[]]
  | n + 1 => (vecs n).map (fun a => false :: a) ++ (vecs n).map (fun a => true :: a)

theorem length_vecs (n : Nat) : (vecs n).length = 2 ^ n := by
  induction n with
  | zero => rfl
  | succ n ih => simp only [vecs, List.length_append, List.length_map, ih, Nat.pow_succ]; omega

/-- One element's weight: its dial if it is in, 1 if it is out. -/
def cell (t : R) : Bool → R
  | true => t
  | false => 1

/-- The weight of a subset: the product of the dials of its members. -/
def W : List R → List Bool → R
  | t :: ts, b :: bs => cell t b * W ts bs
  | _, _ => 1

theorem sumL_map_append {α : Type} (f : α → R) (l m : List α) :
    sumL ((l ++ m).map f) = sumL (l.map f) + sumL (m.map f) := by
  rw [List.map_append, sumL_append]

theorem sumL_map_mul {α : Type} (k : R) (f : α → R) (l : List α) :
    sumL (l.map (fun a => k * f a)) = k * sumL (l.map f) := by
  induction l with
  | nil => simp only [List.map, sumL]; grind
  | cons a l ih => simp only [List.map, sumL] at ih ⊢; rw [ih]; grind

theorem sumL_vecs_succ (n : Nat) (F : List Bool → R) :
    sumL ((vecs (n + 1)).map F)
      = sumL ((vecs n).map (fun a => F (false :: a))) + sumL ((vecs n).map (fun a => F (true :: a))) := by
  simp only [vecs, sumL_map_append, List.map_map]; rfl

theorem sumL_map_zero {α : Type} (l : List α) :
    sumL (l.map (fun _ => (0 : R))) = 0 := by
  induction l with
  | nil => rfl
  | cons a l ih => simp only [List.map, sumL] at ih ⊢; rw [ih]; grind

theorem sumL_map_add {α : Type} (f g : α → R) (l : List α) :
    sumL (l.map (fun a => f a + g a)) = sumL (l.map f) + sumL (l.map g) := by
  induction l with
  | nil => simp only [List.map, sumL]; grind
  | cons a l ih => simp only [List.map, sumL] at ih ⊢; rw [ih]; grind

/-- **The count.** The dialled subsets of an n-element set add up to ∏(1 + tᵢ). -/
theorem total_W (ts : List R) :
    sumL ((vecs ts.length).map (W ts)) = prodL (ts.map (fun t => 1 + t)) := by
  induction ts with
  | nil => show sumL [W [] []] = prodL [] ; simp only [W, sumL, prodL]; grind
  | cons t ts ih =>
    simp only [List.length_cons, vecs, sumL_map_append, List.map_map]
    have e1 : (W (t :: ts) ∘ fun a => false :: a) = fun a => 1 * W ts a := by
      funext a; simp only [Function.comp, W, cell]
    have e2 : (W (t :: ts) ∘ fun a => true :: a) = fun a => t * W ts a := by
      funext a; simp only [Function.comp, W, cell]
    rw [e1, e2, sumL_map_mul, sumL_map_mul, ih]
    simp only [List.map, prodL]; grind

/-- The count is the mult-add chain plus one: ∏(1 + tᵢ) = 1 + (t₁ ◯ … ◯ tₙ). -/
theorem total_W_chain (ts : List R) :
    sumL ((vecs ts.length).map (W ts)) = 1 + foldCirc ts := by
  rw [total_W, foldCirc_lift]

/-- All dials equal: (1 + t)ⁿ. -/
theorem total_uniform (t : R) (n : Nat) :
    prodL ((List.replicate n t).map (fun x => 1 + x)) = (1 + t) ^ n := by
  induction n with
  | zero => simp only [List.replicate, List.map, prodL]; rw [Semiring.pow_zero]
  | succ n ih =>
    simp only [List.replicate, List.map, prodL] at ih ⊢
    rw [ih, Semiring.pow_succ]; grind

end Semiring

/-- All dials 1: every subset weighs 1, and there are 2ⁿ of them. -/
theorem W_ones (bs : List Bool) (n : Nat) : W (List.replicate n (1 : Nat)) bs = 1 := by
  induction n generalizing bs with
  | zero => cases bs <;> rfl
  | succ n ih => cases bs with
    | nil => rfl
    | cons b bs => simp only [List.replicate, W]; rw [ih]; cases b <;> rfl

theorem total_full (n : Nat) :
    sumL ((vecs n).map (W (List.replicate n (1 : Nat)))) = 2 ^ n := by
  have h := total_W (List.replicate n (1 : Nat))
  rw [List.length_replicate] at h
  rw [h, total_uniform]

/-- All dials 0: only the empty set survives, total 1. -/
theorem total_empty (n : Nat) :
    sumL ((vecs n).map (W (List.replicate n (0 : Nat)))) = 1 := by
  have h := total_W (List.replicate n (0 : Nat))
  rw [List.length_replicate] at h
  rw [h, total_uniform]; simp

/-- Dials that are only 0 or 1 give 2^(number of ones): the classical split
    P(S) = P(A) × P(S∖A), the corners of the box of halves. -/
theorem total_split (bs : List Bool) :
    prodL ((bs.map (fun b => if b then (1 : Nat) else 0)).map (fun t => 1 + t))
      = 2 ^ (bs.filter (fun b => b)).length := by
  induction bs with
  | nil => rfl
  | cons b bs ih =>
    cases b <;> simp only [List.map, prodL, List.filter, List.length_cons] at ih ⊢ <;>
      rw [ih] <;> simp [Nat.pow_succ] <;> omega

/-! ## 2. Merging -/

section Merge
variable {R : Type} [Lean.Grind.CommSemiring R]

/-- Element-wise "or" of two bit vectors: the union of two subsets. -/
def orV : List Bool → List Bool → List Bool
  | a :: as, b :: bs => (a || b) :: orV as bs
  | _, _ => []

/-- Union convolution: every pair (A, B) of subsets, weighed W_t(A)·W_s(B),
    lands on A ∪ B. -/
def conv (n : Nat) (ts ss : List R) (c : List Bool) : R :=
  sumL ((vecs n).map fun a => sumL ((vecs n).map fun b =>
    if orV a b = c then W ts a * W ss b else 0))

/-- **One element merged.** The four ways (in/out, in/out) that land on
    "out" weigh 1; the three that land on "in" weigh t + s + ts = t ◯ s.
    Read as coarse-graining: a block of two elements, present when either is,
    carries the dial t ◯ s. -/
theorem cell_merge (t s : R) (z : Bool) :
    ((if (false || false) = z then cell t false * cell s false else 0) +
     (if (false || true) = z then cell t false * cell s true else 0)) +
    ((if (true || false) = z then cell t true * cell s false else 0) +
     (if (true || true) = z then cell t true * cell s true else 0))
      = cell (t ◯ s) z := by
  cases z <;> simp only [Bool.or_false, Bool.or_true, cell, circ] <;>
    simp <;> grind

theorem summand_split (t s : R) (ts ss : List R) (x y z : Bool) (a b c : List Bool) :
    (if orV (x :: a) (y :: b) = z :: c then W (t :: ts) (x :: a) * W (s :: ss) (y :: b) else 0)
      = (if (x || y) = z then cell t x * cell s y else 0) *
        (if orV a b = c then W ts a * W ss b else 0) := by
  simp only [orV, W, List.cons.injEq]
  by_cases h1 : (x || y) = z <;> by_cases h2 : orV a b = c <;> simp [h1, h2] <;> grind

theorem inner_sum (t s : R) (ts ss : List R) (x z : Bool) (c : List Bool) (n : Nat) (a : List Bool) :
    sumL ((vecs (n + 1)).map fun b =>
        if orV (x :: a) b = z :: c then W (t :: ts) (x :: a) * W (s :: ss) b else 0)
      = (if (x || false) = z then cell t x * cell s false else 0) *
          sumL ((vecs n).map fun b => if orV a b = c then W ts a * W ss b else 0) +
        (if (x || true) = z then cell t x * cell s true else 0) *
          sumL ((vecs n).map fun b => if orV a b = c then W ts a * W ss b else 0) := by
  rw [sumL_vecs_succ]
  have e : ∀ y : Bool, (fun b => if orV (x :: a) (y :: b) = z :: c then
      W (t :: ts) (x :: a) * W (s :: ss) (y :: b) else 0)
      = fun b => (if (x || y) = z then cell t x * cell s y else 0) *
          (if orV a b = c then W ts a * W ss b else 0) := by
    intro y; funext b; exact summand_split t s ts ss x y z a b c
  rw [e false, e true, sumL_map_mul, sumL_map_mul]

/-- **Merging dialled power sets multiplies nothing new in: the dials mult-add.**
    Union-convolving the family with dials ts and the family with dials ss
    gives the family with dials tᵢ ◯ sᵢ. -/
theorem conv_W (n : Nat) (ts ss : List R) (c : List Bool)
    (h1 : ts.length = n) (h2 : ss.length = n) (h3 : c.length = n) :
    conv n ts ss c = W (List.zipWith circ ts ss) c := by
  induction n generalizing ts ss c with
  | zero =>
    match ts, ss, c, h1, h2, h3 with
    | [], [], [], _, _, _ => simp [conv, vecs, orV, W, sumL]; grind
  | succ n ih =>
    match ts, ss, c, h1, h2, h3 with
    | t :: ts, s :: ss, z :: c, h1, h2, h3 =>
      simp only [List.length_cons, Nat.add_right_cancel_iff] at h1 h2 h3
      have key := ih ts ss c h1 h2 h3
      unfold conv at key ⊢
      have e : ∀ x : Bool, (fun a => sumL ((vecs (n + 1)).map fun b =>
          if orV (x :: a) b = z :: c then W (t :: ts) (x :: a) * W (s :: ss) b else 0))
          = fun a => (if (x || false) = z then cell t x * cell s false else 0) *
              sumL ((vecs n).map fun b => if orV a b = c then W ts a * W ss b else 0) +
            (if (x || true) = z then cell t x * cell s true else 0) *
              sumL ((vecs n).map fun b => if orV a b = c then W ts a * W ss b else 0) := by
        intro x; funext a; exact inner_sum t s ts ss x z c n a
      rw [sumL_vecs_succ, e false, e true]
      rw [sumL_map_add, sumL_map_add, sumL_map_mul, sumL_map_mul, sumL_map_mul, sumL_map_mul, key]
      simp only [List.zipWith, W]
      rw [← cell_merge t s z]; grind

/-- Dial 0 is the unit of merging: the family {∅}. -/
theorem zipWith_zero (ss : List R) :
    List.zipWith circ (List.replicate ss.length 0) ss = ss := by
  induction ss with
  | nil => rfl
  | cons s ss ih => simp only [List.length_cons, List.replicate, List.zipWith]; rw [ih, zero_circ]

/-- Merging a constant family with itself mult-adds the dial with itself. -/
theorem zipWith_replicate (t : R) (n : Nat) :
    List.zipWith circ (List.replicate n t) (List.replicate n t) = List.replicate n (t ◯ t) := by
  induction n with
  | zero => rfl
  | succ n ih => simp only [List.replicate, List.zipWith]; rw [ih]

/-- Merging commutes with halving: if r halves t and r' halves t', then r ◯ r'
    halves t ◯ t'. (On the line u = log(1 + t), merging is addition and halving
    is u ↦ u/2.) -/
theorem half_of_merge (r r' : R) :
    (r ◯ r') ◯ (r ◯ r') = (r ◯ r) ◯ (r' ◯ r') := by
  unfold circ; grind

/-- m copies merged: the m-fold mult-add power of the dial. An m-th root of the
    full power set is a dial with t^(◯m) = 1, that is (1 + t)^m = 2. -/
theorem mth_root_of_full (t : R) (m : Nat) (h : (1 + t) ^ m = 1 + 1) :
    1 + mpow t m = 1 + 1 := by
  rw [mpow_lift, h]

end Merge

/-! ## 3. The half -/

section Half
variable {R : Type} [Lean.Grind.CommRing R]

/-- **The half dial.** If r·r = 2, then h = r − 1 (that is √2 − 1) has h ◯ h = 1. -/
theorem half_dial (r : R) (hr : r * r = 2) : (r - 1) ◯ (r - 1) = 1 := by
  unfold circ; grind

/-- W with every dial 1 weighs every subset 1. -/
theorem W_replicate_one (n : Nat) (c : List Bool) : W (List.replicate n (1 : R)) c = 1 := by
  induction n generalizing c with
  | zero => cases c <;> rfl
  | succ n ih => cases c with
    | nil => rfl
    | cons b c => simp only [List.replicate, W]; rw [ih]; cases b <;> simp only [cell] <;> grind

/-- **The half power set, squared, is the power set.** With every dial √2 − 1,
    the family merged with itself weighs every subset exactly 1: it is P(S). -/
theorem half_squared_family (r : R) (hr : r * r = 2) (n : Nat) (c : List Bool)
    (hc : c.length = n) :
    conv n (List.replicate n (r - 1)) (List.replicate n (r - 1)) c = 1 := by
  rw [conv_W n _ _ c (List.length_replicate) (List.length_replicate) hc,
      zipWith_replicate, half_dial r hr, W_replicate_one]

end Half

section Field
variable {K : Type} [Lean.Grind.Field K]

/-- The partner of a dial: the dial that completes it to the full power set. -/
def partner (t : K) : K := (1 - t) / (1 + t)

/-- Partners mult-add to 1: (1 + t)(1 + partner t) = 2. -/
theorem partner_spec (t : K) (h : 1 + t ≠ 0) : t ◯ partner t = 1 := by
  unfold partner circ; grind

/-- Partnering is an involution. -/
theorem partner_partner (t : K) (h : 1 + t ≠ 0) (h2 : (1 + 1 : K) ≠ 0) :
    partner (partner t) = t := by
  unfold partner
  have e : 1 + (1 - t) / (1 + t) = (1 + 1) / (1 + t) := by grind
  grind

/-- The corners swap: dial 0 (leave the element out) pairs with dial 1 (all of it). -/
theorem partner_zero : partner (0 : K) = 1 := by unfold partner; grind
theorem partner_one : partner (1 : K) = 0 := by unfold partner; grind

/-- The fixed points of partnering are exactly the half dials. -/
theorem partner_fixed_iff (t : K) (h : 1 + t ≠ 0) : partner t = t ↔ t ◯ t = 1 := by
  unfold partner circ
  constructor
  · intro e; grind
  · intro e; grind

/-- In probabilities (p = t/(1 + t)), the union half is p = 1 − r/2 = 1 − 1/√2:
    two independent copies, joined, give a fair coin, 1/2. -/
theorem half_prob_union (r : K) (hr : r * r = 2) (h2 : (1 + 1 : K) ≠ 0) :
    circK (-1) (1 - r / 2) (1 - r / 2) = 1 / 2 := by
  unfold circK; grind

/-! ## 4. Complement -/

/-- Complement of one element: in becomes out. Weight t on "in" read from the
    other side is t times weight 1/t on "out". -/
theorem cell_compl (t : K) (ht : t ≠ 0) (b : Bool) : cell t (!b) = t * cell t⁻¹ b := by
  cases b <;> simp only [Bool.not_true, Bool.not_false, cell] <;> grind

/-- **Complement inverts every dial.** The family with dials t, read through
    A ↦ S ∖ A, is (∏ t) times the family with dials 1/t. -/
theorem W_compl (ts : List K) (c : List Bool) (h : ∀ t ∈ ts, t ≠ 0)
    (hc : c.length = ts.length) :
    W ts (c.map (fun b => !b)) = prodL ts * W (ts.map (fun t => t⁻¹)) c := by
  induction ts generalizing c with
  | nil => cases c <;> simp [W, prodL] <;> grind
  | cons t ts ih =>
    cases c with
    | nil => simp at hc
    | cons b c =>
      simp only [List.map, W, prodL]
      have ht : t ≠ 0 := h t (by simp)
      have hts : ∀ x ∈ ts, x ≠ 0 := fun x hx => h x (by simp [hx])
      simp only [List.length_cons, Nat.add_right_cancel_iff] at hc
      rw [ih c hts hc, cell_compl t ht b]; grind

/-- Complement swaps the two mergers: reciprocal dials mult-add under
    intersection exactly as dials mult-add under union.
    (1 + 1/t)(1 + 1/s) is the lift of t·s/(1 + t + s), the intersection's dial. -/
theorem compl_swaps_merge (t s : K) (ht : t ≠ 0) (hs : s ≠ 0) (h : 1 + t + s ≠ 0) :
    t⁻¹ ◯ s⁻¹ = (t * s / (1 + t + s))⁻¹ := by
  unfold circ
  have hts : t * s ≠ 0 := by grind
  grind

/-- The intersection half is the reciprocal of the union half: √2 + 1 = 1/(√2 − 1). -/
theorem halves_reciprocal (r : K) (hr : r * r = 2) : (r - 1) * (r + 1) = 1 := by
  grind

/-- In probabilities the two halves are 1 − r/2 and r/2: mirror images about 1/2. -/
theorem halves_mirror (r : K) : (1 - r / 2) + r / 2 = 1 := by grind

/-- The intersection half: two independent copies, intersected, give 1/2. -/
theorem half_prob_inter (r : K) (hr : r * r = 2) (h2 : (1 + 1 : K) ≠ 0) :
    (r / 2) * (r / 2) = 1 / 2 := by grind

end Field

/-- **No single probability is a half both ways.** Joining two copies gives
    1/2 and intersecting two copies gives 1/2 cannot both hold. -/
theorem no_half_both_ways (p : Rat) :
    ¬ (circK (-1) p p = 1 / 2 ∧ p * p = 1 / 2) := by
  intro ⟨h1, h2⟩
  unfold circK at h1
  have hp : p = 1 / 2 := by grind
  rw [hp] at h2; grind

/-- **The only dial complement leaves alone is 1** (the full power set). -/
theorem compl_fixed_dial (t : Rat) (hpos : 0 < t) (h : t⁻¹ = t) : t = 1 := by
  have ht : t ≠ 0 := by grind
  have e : t * t = 1 := by
    have : t * t⁻¹ = 1 := Rat.mul_inv_cancel t ht
    rw [h] at this; exact this
  have h1 : (t - 1) * (t + 1) = 0 := by grind
  rcases Rat.mul_eq_zero.mp h1 with h1 | h1 <;> grind

/-! ## 5. No half step by iteration -/

/-- **No dial construction is a half step of the power set by iteration.**
    The empty set has exactly one subset (itself), weight 1, so any size
    rule f coming from dials has f 0 = 1. Then f (f n) = 2ⁿ fails at once:
    f 1 = f (f 0) = 1, so f (f 1) = 1, not 2. No monotonicity is assumed. -/
theorem no_iterate_half (f : Nat → Nat) (h0 : f 0 = 1) : ¬ (∀ n, f (f n) = 2 ^ n) := by
  intro h
  have a := h 0; have b := h 1
  rw [h0] at a; simp at a
  rw [a] at b; simp at b; omega

/-- The same over the rationals (real-valued masses): f 0 = 1 forces f 1 = 1. -/
theorem no_iterate_half_rat (f : Rat → Rat) (h0 : f 0 = 1) :
    ¬ (f (f 0) = 1 ∧ f (f 1) = 2) := by
  intro ⟨a, b⟩
  rw [h0] at a; rw [a] at b; grind

/-! ## 6. Known sets: primes, naturals, integers, rationals -/

section KnownSets
variable {R : Type} [Lean.Grind.CommRing R]

/-- **Primes → squarefree numbers.** One prime's dial x = p^(−s) contributes the
    factor 1 + x; times 1 − x it is 1 − x², the same factor at 2s. Over all primes
    this is ∏(1 + p^(−s)) = ζ(s)/ζ(2s). -/
theorem euler_squarefree (x : R) : (1 + x) * (1 - x) = 1 - x * x := by grind

/-- **Naturals as merged squarefree layers (Fermi–Dirac).** The exponent of a
    prime written in binary: (1 + x)(1 + x²)(1 + x⁴)⋯(1 + x^(2^(k−1))) times (1 − x)
    is 1 − x^(2^k). So ζ(s) is the merge of the squarefree power sets at dials
    p^(−s), p^(−2s), p^(−4s), … -/
theorem fermi_dirac (x : R) (k : Nat) :
    (1 - x) * prodL ((List.range k).map (fun j => 1 + x ^ (2 ^ j))) = 1 - x ^ (2 ^ k) := by
  induction k with
  | zero => simp only [List.range_zero, List.map, prodL, Nat.pow_zero]; rw [Semiring.pow_succ, Semiring.pow_zero]; grind
  | succ k ih =>
    rw [List.range_succ, List.map_append]
    have hp : ∀ (l m : List R), prodL (l ++ m) = prodL l * prodL m := by
      intro l m; induction l with
      | nil => simp only [List.nil_append, prodL]; grind
      | cons a l ih2 => simp only [List.cons_append, prodL]; rw [ih2]; grind
    rw [hp]; simp only [List.map, prodL]
    have e : x ^ (2 ^ (k + 1)) = x ^ (2 ^ k) * x ^ (2 ^ k) := by
      rw [Nat.pow_succ, Nat.mul_two, pow_add']
    rw [e]; grind

/-- **Rationals: the exponent of a prime runs over all of ℤ.** Weighing exponent e
    by x^|e|, the first 2m + 1 exponents sum to 1 + 2(x + ⋯ + x^m), and times (1 − x)
    that is 1 + x − 2x^(m+1). The full factor is (1 + x)/(1 − x): the squarefree
    factor merged with the naturals' factor, so ∏ = ζ(s)²/ζ(2s). -/
theorem rational_factor (x : R) (m : Nat) :
    (1 - x) * (1 + 2 * sumL ((List.range m).map (fun k => x ^ (k + 1))))
      = 1 + x - 2 * x ^ (m + 1) := by
  induction m with
  | zero => simp only [List.range_zero, List.map, sumL]; rw [Semiring.pow_succ, Semiring.pow_zero]; grind
  | succ m ih =>
    rw [List.range_succ, List.map_append, sumL_append]
    simp only [List.map, sumL]
    have e : x ^ (m + 1 + 1) = x ^ (m + 1) * x := by rw [Semiring.pow_succ]
    rw [e]; grind

/-- **Integers: ℕ is an exact half of ℤ∖{0}.** Pair n with −n, both carrying the
    same dial: the total is the square of the naturals' total, so the naturals
    (a corner half) and the centre half have exactly the same total. -/
theorem integers_square (ts : List R) :
    prodL (ts.map (fun t => (1 + t) * (1 + t))) = prodL (ts.map (fun t => 1 + t)) * prodL (ts.map (fun t => 1 + t)) := by
  induction ts with
  | nil => simp only [List.map, prodL]; grind
  | cons t ts ih => simp only [List.map, prodL]; rw [ih]; grind

end KnownSets

section Tiny
variable {K : Type} [Lean.Grind.Field K]

/-- **Tiny dials linearise.** Halving a small dial x by plain participation (x/2)
    misses the union half only at second order: (x/2) ◯ (x/2) = x + x²/4. Where the
    dials that matter are tiny (large primes, Poisson points on the line), the
    naive half and the merging half agree. -/
theorem tiny_half (x : K) (h2 : (2 : K) ≠ 0) : (x / 2) ◯ (x / 2) = x + x * x / 4 := by
  unfold circ; grind

/-- Participation complement c ↦ 1 − c fixes exactly c = 1/2. -/
theorem participation_fixed (c : K) (h2 : (2 : K) ≠ 0) : c = 1 - c ↔ c = 1 / 2 := by
  constructor
  · intro e; grind
  · intro e; grind

/-- Exponents of merged pieces add, and complement of participation sends α to
    1 − α; the two halves α = 1/2 and 1 − α merge to the whole. -/
theorem halves_merge_to_whole (α : K) (h : α = 1 / 2) (h2 : (2 : K) ≠ 0) : α + (1 - α) = 1 ∧ α = 1 - α := by
  constructor
  · grind
  · rw [h]; grind

end Tiny

/-- The 32 subsets of the first five primes give 32 different squarefree numbers. -/
theorem squarefree_from_subsets :
    ((groups [2, 3, 5, 7, 11]).map (fun g => g.foldl (· * ·) 1)).Nodup := by decide

end PowersetDials
