/-! Round 11 (DRAFT): small facts for README section 15.
    NOT BUILT this round: the Lean toolchain could not be installed in this session. These are written to be
    core-only like NewTies.lean and should be checked with `lean Round11.lean` before anyone relies on them.
    Built with Lean 4.37 (core only) in the repo on Nov 27: all five check, no proof changed. In the library as
    FusedAlgebra.Topics.NewTies11 (namespace CircleOp.NewTies11). -/

namespace CircleOp.NewTies11

/-- y-model: N(y) = (y+1)/19 and K(y) = (y+1)/10, so 19·N = 10·K (N is 10/19 times K). -/
theorem ymodel_N_K (y n k : Int) (hn : 19 * n = y + 1) (hk : 10 * k = y + 1) : 19 * n = 10 * k := by omega

/-- 19 never divides 200^a. -/
theorem pow200_mod19 (a : Nat) : 200 ^ a % 19 ≠ 0 := by
  induction a with
  | zero => decide
  | succ n ih =>
    rw [Nat.pow_succ]
    omega

/-- The stretch is never exactly 1: λ^a = μ^b is impossible for a ≥ 1 (19-adic count). -/
theorem rho_ne_one (a b : Nat) (ha : 1 ≤ a) : 19 ^ a * 20 ^ b ≠ 200 ^ a := by
  intro h
  have h1 : 19 ^ a * 20 ^ b % 19 = 0 := by
    obtain ⟨n, rfl⟩ : ∃ n, a = n + 1 := ⟨a - 1, by omega⟩
    rw [Nat.pow_succ, Nat.mul_assoc, Nat.mul_comm 19, ← Nat.mul_assoc]
    exact Nat.mul_mod_left _ _
  have h2 := pow200_mod19 a
  rw [h] at h1
  exact h2 h1

/-- q = λ¹⁴/μ¹¹ lies in (0.998, 0.999): the bad pairs (14k, 11k) drift through the N window slowly. -/
theorem family_ratio :
    998 * 200 ^ 14 < 1000 * (19 ^ 14 * 20 ^ 11) ∧ 1000 * (19 ^ 14 * 20 ^ 11) < 999 * 200 ^ 14 := by decide

/-- (14, 11) is a bad pair: λ¹⁴/μ¹¹ lies in the N window [9/10, 10/9]. -/
theorem bad_pair_14_11 :
    9 * 200 ^ 14 ≤ 10 * (19 ^ 14 * 20 ^ 11) ∧ 9 * (19 ^ 14 * 20 ^ 11) ≤ 10 * 200 ^ 14 := by decide

end CircleOp.NewTies11
