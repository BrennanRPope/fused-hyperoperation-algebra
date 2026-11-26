import Mathlib.Analysis.SpecificLimits.Normed
import MathlibMap.Walls

/-!
# Near dial zero there are no walls (Neumann series)

In a complete normed algebra (non-unital is fine), Mathlib's `IsQuasiregular.of_norm_lt_one`
says every element of norm below 1 is quasiregular: its quasi-inverse is the Neumann series.
Through `hasWayBack_iff_isQuasiregular` this becomes a statement about our dial:

* `hasWayBack_of_norm_smul_lt_one`: if `‖k • a‖ < 1` then `a` has a way back at dial `k`.
* `hasWayBack_of_norm_mul_lt_one`: so with `‖k‖ · ‖a‖ < 1` no wall can sit within `1/‖a‖` of
  dial zero.
-/

namespace MathlibMap

variable {𝕜 A : Type*} [NontriviallyNormedField 𝕜] [NonUnitalNormedRing A] [NormedSpace 𝕜 A]
  [IsScalarTower 𝕜 A A] [SMulCommClass 𝕜 A A] [CompleteSpace A]

/-- A way back at dial `k` whenever `‖k • a‖ < 1` (Mathlib's Neumann series). -/
theorem hasWayBack_of_norm_smul_lt_one {k : 𝕜} (hk : k ≠ 0) {a : A} (h : ‖k • a‖ < 1) :
    HasWayBack k a :=
  (hasWayBack_iff_isQuasiregular hk a).mpr (IsQuasiregular.of_norm_lt_one h)

/-- No wall within `1/‖a‖` of dial zero. -/
theorem hasWayBack_of_norm_mul_lt_one {k : 𝕜} {a : A} (h : ‖k‖ * ‖a‖ < 1) : HasWayBack k a := by
  by_cases hk : k = 0
  · subst hk; exact hasWayBack_zero a
  · exact hasWayBack_of_norm_smul_lt_one hk ((norm_smul_le k a).trans_lt h)

end MathlibMap
