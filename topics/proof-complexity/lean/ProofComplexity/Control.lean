import Mathlib.Algebra.BigOperators.Group.Finset.Powerset
import Mathlib.Combinatorics.Enumerative.InclusionExclusion
import Mathlib.Data.Nat.Choose.Sum
import Mathlib.Algebra.Ring.GeomSum
import Mathlib.Tactic.LinearCombination

/-!
# Controls: the operation route's proof strategy, with the operation taken out

Each control follows the same plan as the matching `op_` proof (peel one element off at
a time, or take the product of complements) but never mentions the circle operation, its
chains or its lift. Comparing `op_` with `ctl_` isolates what the operation itself adds.
-/

namespace HCT.Control
open Finset

theorem ctl_prod_one_add {ι R : Type*} [CommRing R] {f : ι → R} (s : Finset ι) :
    ∏ i ∈ s, (1 + f i) = ∑ t ∈ s.powerset, ∏ i ∈ t, f i := by
  classical
  induction s using Finset.induction_on with
  | empty => simp
  | insert a s ha ih =>
    rw [prod_insert ha, ih, sum_powerset_insert ha, add_mul, one_mul, mul_sum]
    congr 1
    refine sum_congr rfl fun t ht => ?_
    have hat : a ∉ t := fun h => ha (mem_powerset.1 ht h)
    rw [prod_insert hat]

theorem ctl_inclusion_exclusion_card_biUnion {ι α : Type*} [DecidableEq α]
    (s : Finset ι) (S : ι → Finset α) :
    #(s.biUnion S) = ∑ t : s.powerset.filter (·.Nonempty),
      (-1 : ℤ) ^ (#t.val + 1) * #(t.val.inf' (mem_filter.1 t.property).2 S) := by
  classical
  let x : ι → α → ℤ := fun i a => if a ∈ S i then 1 else 0
  have h2 : ∀ a, 1 - ∏ i ∈ s, (1 + -x i a) =
      ∑ t ∈ s.powerset with t.Nonempty, (-1) ^ (#t + 1) * ∏ i ∈ t, x i a := by
    intro a
    rw [ctl_prod_one_add, ← sum_filter_add_sum_filter_not s.powerset Finset.Nonempty]
    have h0 : ∑ t ∈ s.powerset with ¬ t.Nonempty, ∏ i ∈ t, -x i a = 1 := by
      rw [sum_eq_single_of_mem ∅ (by simp)] <;> simp_all [not_nonempty_iff_eq_empty]
    rw [h0]
    simp only [prod_neg]
    rw [sub_add_cancel_right, ← sum_neg_distrib]
    refine sum_congr rfl fun t _ => by ring
  have h1 : ∀ a ∈ s.biUnion S, 1 - ∏ i ∈ s, (1 + -x i a) = 1 := by
    intro a ha
    obtain ⟨i, hi, hai⟩ := mem_biUnion.1 ha
    rw [prod_eq_zero hi (by simp [x, hai])]; ring
  calc (#(s.biUnion S) : ℤ) = ∑ a ∈ s.biUnion S, (1 - ∏ i ∈ s, (1 + -x i a)) := by
        rw [sum_congr rfl h1]; simp
    _ = ∑ t ∈ s.powerset with t.Nonempty,
          (-1) ^ (#t + 1) * ∑ a ∈ s.biUnion S, ∏ i ∈ t, x i a := by
        simp_rw [h2, mul_sum]; exact sum_comm
    _ = _ := by
        rw [← sum_coe_sort]
        refine sum_congr rfl fun ⟨t, ht⟩ _ => ?_
        obtain ⟨hts, ⟨j, hj⟩⟩ := mem_filter.1 ht
        simp only [x, prod_boole, sum_boole]
        congr 3; ext a
        simp only [mem_filter, mem_biUnion, mem_inf']
        exact ⟨fun h => h.2, fun h => ⟨⟨j, mem_powerset.1 hts hj, h j hj⟩, h⟩⟩

theorem ctl_sum_powerset_neg_one_pow_card_of_nonempty {α : Type*} {x : Finset α}
    (h0 : x.Nonempty) : (∑ m ∈ x.powerset, (-1 : ℤ) ^ #m) = 0 := by
  obtain ⟨i, hi⟩ := h0
  have h := ctl_prod_one_add (f := fun _ : α => (-1 : ℤ)) x
  rw [prod_eq_zero hi (by norm_num)] at h
  simp only [prod_const] at h
  exact h.symm

theorem ctl_one_add_pow {R : Type*} [CommRing R] (x : R) (n : ℕ) :
    (1 + x) ^ n = ∑ m ∈ range (n + 1), x ^ m * n.choose m := by
  have h := ctl_prod_one_add (f := fun _ : ℕ => x) (range n)
  simp only [prod_const, card_range] at h
  rw [h, sum_powerset_apply_card (fun m => x ^ m), card_range]
  simp [nsmul_eq_mul, mul_comm]

theorem ctl_geom_sum_mul {R : Type*} [CommRing R] (x : R) (n : ℕ) :
    (∑ i ∈ range n, x ^ i) * (x - 1) = x ^ n - 1 := by
  induction n with
  | zero => simp
  | succ n ih => rw [sum_range_succ, add_mul, ih]; ring

end HCT.Control
