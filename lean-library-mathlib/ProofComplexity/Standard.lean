import Mathlib.Data.Nat.Choose.Sum

/-! Standard proofs where Mathlib states the result slightly differently. -/

namespace HCT.Results
open Finset

/-- R4, standard route: Mathlib's binomial theorem `add_pow`, specialised. -/
theorem std_one_add_pow {R : Type*} [CommRing R] (x : R) (n : ℕ) :
    (1 + x) ^ n = ∑ m ∈ range (n + 1), x ^ m * n.choose m := by
  rw [add_comm, add_pow]; simp

end HCT.Results
