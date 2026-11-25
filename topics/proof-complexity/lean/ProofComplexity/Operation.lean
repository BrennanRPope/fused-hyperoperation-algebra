import ProofComplexity.Framework
import Mathlib.Combinatorics.Enumerative.InclusionExclusion
import Mathlib.Data.Nat.Choose.Sum
import Mathlib.Algebra.Ring.GeomSum
import Mathlib.Analysis.Complex.Trigonometric

/-!
# Known results proved through the circle operation

Each theorem has exactly the statement of a standard result (Mathlib's statement where
Mathlib has one) and is proved by going through `HCT.chain`, its lift or its corners.
-/

namespace HCT.Results
open Finset

/-- R1. The subset expansion of a product (Mathlib: `Finset.prod_one_add`). -/
theorem op_prod_one_add {ι R : Type*} [CommRing R] {f : ι → R} (s : Finset ι) :
    ∏ i ∈ s, (1 + f i) = ∑ t ∈ s.powerset, ∏ i ∈ t, f i := by
  classical
  have h1 := chain_lift (1 : R) s f
  have h2 := chain_corners (1 : R) s f
  simp only [one_mul, one_pow] at h1 h2
  rw [← h1, h2]

/-- R2. Inclusion–exclusion for the size of a union
(Mathlib: `Finset.inclusion_exclusion_card_biUnion`, same statement). On the union, the
k = −1 chain of the indicators is 1; its corners are the alternating intersections. -/
theorem op_inclusion_exclusion_card_biUnion {ι α : Type*} [DecidableEq α]
    (s : Finset ι) (S : ι → Finset α) :
    #(s.biUnion S) = ∑ t : s.powerset.filter (·.Nonempty),
      (-1 : ℤ) ^ (#t.val + 1) * #(t.val.inf' (mem_filter.1 t.property).2 S) := by
  classical
  let x : ι → α → ℤ := fun i a => if a ∈ S i then 1 else 0
  have h1 : ∀ a ∈ s.biUnion S, chain (-1) s (x · a) = 1 := by
    intro a ha
    obtain ⟨i, hi, hai⟩ := mem_biUnion.1 ha
    have h := chain_lift (-1 : ℤ) s (x · a)
    rw [prod_eq_zero hi (by simp [x, hai])] at h
    linarith
  have h2 : ∀ a, chain (-1) s (x · a) =
      ∑ t ∈ s.powerset with t.Nonempty, (-1) ^ (#t + 1) * ∏ i ∈ t, x i a := by
    intro a
    have h := chain_corners_nonempty (-1 : ℤ) s (x · a)
    simp only [pow_succ, mul_neg_one, neg_mul, sum_neg_distrib, ← h]; ring
  calc (#(s.biUnion S) : ℤ) = ∑ a ∈ s.biUnion S, chain (-1) s (x · a) := by
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

/-- R3. The heart of Möbius inversion on the cube: the signed count of the corners of a
non-empty cube is 0 (Mathlib: `Finset.sum_powerset_neg_one_pow_card_of_nonempty`).
Through the operation: −1 is the absorbing element of `∘₁`, since `1 + 1·(−1) = 0`. -/
theorem op_sum_powerset_neg_one_pow_card_of_nonempty {α : Type*} {x : Finset α}
    (h0 : x.Nonempty) : (∑ m ∈ x.powerset, (-1 : ℤ) ^ #m) = 0 := by
  classical
  have h1 := chain_lift (1 : ℤ) x (fun _ => -1)
  have h2 := chain_corners (1 : ℤ) x (fun _ => -1)
  obtain ⟨i, hi⟩ := h0
  rw [prod_eq_zero hi (by norm_num)] at h1
  simp only [one_pow, one_mul, prod_const] at h1 h2
  rw [← h2, h1]

/-- R4. The binomial theorem for `(1 + x)ⁿ`. Through the operation: n equal steps,
lifted, and the cube's corners counted by rank. -/
theorem op_one_add_pow {R : Type*} [CommRing R] (x : R) (n : ℕ) :
    (1 + x) ^ n = ∑ m ∈ range (n + 1), x ^ m * n.choose m := by
  have h1 := chain_lift (1 : R) (range n) (fun _ => x)
  have h2 := chain_corners (1 : R) (range n) (fun _ => x)
  simp only [one_mul, one_pow, prod_const, card_range] at h1 h2
  rw [← h1, h2, sum_powerset_apply_card (fun m => x ^ m), card_range]
  simp [nsmul_eq_mul, mul_comm]

/-- R5. The geometric series (Mathlib: `geom_sum_mul`, same statement). The chain of n
copies of `x − 1` unrolls into the geometric sum, and its lift is `xⁿ`. -/
theorem op_geom_sum_mul {R : Type*} [CommRing R] (x : R) (n : ℕ) :
    (∑ i ∈ range n, x ^ i) * (x - 1) = x ^ n - 1 := by
  have hc : ∀ n, chain (1 : R) (range n) (fun _ => x - 1) = (∑ i ∈ range n, x ^ i) * (x - 1) := by
    intro n
    induction n with
    | zero => simp
    | succ n ih =>
      rw [sum_range_succ', range_add_one, chain_insert _ _ notMem_range_self, ih]
      simp only [pow_succ', ← mul_sum, pow_zero]
      unfold circK; ring
  have h1 := chain_lift (1 : R) (range n) (fun _ => x - 1)
  simp only [one_mul, prod_const, card_range, add_sub_cancel] at h1
  rw [← hc]; linear_combination h1

end HCT.Results
