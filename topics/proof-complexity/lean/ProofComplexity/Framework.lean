import Mathlib.Algebra.BigOperators.Group.Finset.Powerset
import Mathlib.Algebra.BigOperators.Ring.Finset
import Mathlib.Tactic.Ring
import Mathlib.Tactic.LinearCombination

/-!
# The framework's setup: the weighted circle operation and its chains

`circK k a b = a + b + k·a·b`, chained over a finite set. Everything a proof
"through the operation" is allowed to use is built here, and its cost is counted.
-/

namespace HCT
open Finset

variable {R : Type*} [CommRing R]

/-- The weighted circle operation. -/
def circK (k a b : R) : R := a + b + k * a * b

instance (k : R) : Std.Commutative (circK k) := ⟨fun a b => by unfold circK; ring⟩
instance (k : R) : Std.Associative (circK k) := ⟨fun a b c => by unfold circK; ring⟩

/-- The lift: `x ↦ 1 + k·x` turns the circle operation into multiplication. -/
theorem lift (k a b : R) : 1 + k * circK k a b = (1 + k * a) * (1 + k * b) := by
  unfold circK; ring

variable {ι : Type*} [DecidableEq ι]

/-- A chain `f i₁ ∘ₖ f i₂ ∘ₖ …` over a finite set (order does not matter). -/
def chain (k : R) (s : Finset ι) (f : ι → R) : R := s.fold (circK k) 0 f

@[simp] theorem chain_empty (k : R) (f : ι → R) : chain k ∅ f = 0 := fold_empty

theorem chain_insert (k : R) {s : Finset ι} {a : ι} (f : ι → R) (ha : a ∉ s) :
    chain k (insert a s) f = circK k (f a) (chain k s f) := fold_insert ha

/-- Lifted, a chain is a product. -/
theorem chain_lift (k : R) (s : Finset ι) (f : ι → R) :
    1 + k * chain k s f = ∏ i ∈ s, (1 + k * f i) := by
  induction s using Finset.induction_on with
  | empty => simp
  | insert a s ha ih => rw [chain_insert k f ha, lift, ih, prod_insert ha]

/-- The cube: lifted, a chain is a sum over all corners `t ⊆ s`, a corner of rank `m`
weighted `kᵐ`. Proved by peeling one step off the chain at a time. -/
theorem chain_corners (k : R) (s : Finset ι) (f : ι → R) :
    1 + k * chain k s f = ∑ t ∈ s.powerset, k ^ #t * ∏ i ∈ t, f i := by
  induction s using Finset.induction_on with
  | empty => simp
  | insert a s ha ih =>
    rw [chain_insert k f ha, lift, ih, sum_powerset_insert ha, add_mul, one_mul, mul_sum]
    congr 1
    refine sum_congr rfl fun t ht => ?_
    have hat : a ∉ t := fun h => ha (mem_powerset.1 ht h)
    rw [card_insert_of_notMem hat, prod_insert hat]; ring

/-- The cube without the lift: the non-empty corners, a rank-m corner weighted k^(m-1)
(stated multiplied through by k, so it holds at k = 0 too). -/
theorem chain_corners_nonempty (k : R) (s : Finset ι) (f : ι → R) :
    k * chain k s f = ∑ t ∈ s.powerset with t.Nonempty, k ^ #t * ∏ i ∈ t, f i := by
  have h := chain_corners k s f
  rw [← sum_filter_add_sum_filter_not s.powerset Finset.Nonempty] at h
  have h0 : ∑ t ∈ s.powerset with ¬ t.Nonempty, k ^ #t * ∏ i ∈ t, f i = 1 := by
    rw [sum_eq_single_of_mem ∅ (by simp)] <;> simp_all [not_nonempty_iff_eq_empty]
  rw [h0] at h; linear_combination h

end HCT
