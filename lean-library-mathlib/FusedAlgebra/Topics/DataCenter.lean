module

public import FusedAlgebra.Basic

/-!
# Data center energy as mult-add

A stage's overhead odds are `r = (input - output) / output`, so `input = output * (1 + r)`, and
a chain of stages has overhead the mult-add of its stages.

## Main results

* `multAddList_add_one`, `multAddList_swap`, `multAddList_split`: a power chain's overhead is the
  mult-add of its stages, in any order, and splits at any meter (`TUE = PUE * ITUE`).
* `pue_corner`, `additive_books_corner`, `loss_costs_one_add`: the corner `p * c` of
  `PUE = (1 + p) * (1 + c)` is cooling the electrical losses.
* `unit_gap`, `add_unit_iff`, `equal_share`: paralleled units with loss `a + c * x ^ 2` each; add
  a unit exactly when `X ^ 2` passes `(a / c) * k * (k + 1)`.
* `cube_layers`, `spread`, `spare_fan`, `spare_fan_helps`: fans, with power proportional to
  flow cubed.
* `stack`: kept fractions multiply.

Every statement is a polynomial identity or inequality. The core-Lean original stated them over
`ℕ` and `ℤ` ("whole units"); here they hold in any commutative (semi)ring, ordered where an
inequality is involved, so in particular over `ℝ`. The product of lifts `lifts rs` of the
original is `(rs.map (· + 1)).prod`.
-/

@[expose] public section

namespace DataCenterFused

open MultAdd

/-! ### Lead 1. A power chain's overhead is the mult-add of its stages, in any order -/

section Chain

variable {R : Type*}

/-- The list lift with the `1` on the right: `(x₁ ✳ ⋯ ✳ xₙ) + 1 = (x₁ + 1) * ⋯ * (xₙ + 1)`. -/
theorem multAddList_add_one [Semiring R] (rs : List R) :
    multAddList rs + 1 = (rs.map (· + 1)).prod := by
  rw [add_comm, one_add_multAddList]
  exact congrArg List.prod (List.map_congr_left fun x _ ↦ add_comm 1 x)

variable [CommSemiring R]

/-- Swapping two neighbouring stages never changes the total: blame by position is not
physical. -/
theorem multAddList_swap (xs ys : List R) (a b : R) :
    multAddList (xs ++ a :: b :: ys) = multAddList (xs ++ b :: a :: ys) :=
  multAddList_perm ((List.Perm.swap b a ys).append_left xs)

/-- Splitting the chain at the IT meter: outside mult-added to inside is the whole chain
(`TUE = PUE * ITUE`). -/
theorem multAddList_split (xs ys : List R) :
    multAddList (xs ++ ys) = multAddList xs ✳ multAddList ys :=
  multAddList_append xs ys

/-- Position independence: scaling one stage's factor `f` by `k` scales the facility total by
`k`, whatever the product `x` of the stages before it and `y` of the stages after it. -/
theorem scale_anywhere (x y f k : R) : x * (f * k) * y = k * (x * f * y) := by ring

/-! ### Lead 2. PUE overhead is electrical mult-add cooling -/

/-- The corner `p * c` is cooling the losses. -/
theorem pue_corner (p c : R) : (1 + p) * (1 + c) = 1 + p + c + p * c := by ring

/-- The additive breakdown `PUE = 1 + PLF + CLF` books the corner as cooling:
`CLF = c + p * c`. -/
theorem additive_books_corner (it p c : R) :
    it * (1 + p) * (1 + c) = it + it * p + (it * c + it * p * c) := by ring

/-- One watt of electrical loss removed saves `1 + c` watts at the meter. -/
theorem loss_costs_one_add (it L c : R) :
    (it + (L + 1)) * (1 + c) = (it + L) * (1 + c) + (1 + c) := by ring

end Chain

/-! ### Lead 3. Paralleled units with loss `a + c * x ^ 2` each

`k` units sharing load `X` lose `k * a + c * X ^ 2 / k`. Multiply through by `k * (k + 1)` to
avoid division. -/

section Units

variable {R : Type*} [CommRing R]

/-- Exact gap: `k * (k + 1) * (cost k - cost (k + 1)) = c * X ^ 2 - k * (k + 1) * a`. -/
theorem unit_gap (k a c X : R) :
    (k * k * (k + 1) * a + (k + 1) * c * X ^ 2) - (k * (k + 1) * (k + 1) * a + k * c * X ^ 2)
      = c * X ^ 2 - k * (k + 1) * a := by ring

/-- Add a unit exactly when `X ^ 2` passes `(a / c) * k * (k + 1)`: the switch point is the
geometric mean of `k` and `k + 1` in units of the sweet spot `√(a / c)`. -/
theorem add_unit_iff [PartialOrder R] [IsOrderedRing R] (k a c X : R) :
    k * (k + 1) * (k + 1) * a + k * c * X ^ 2 < k * k * (k + 1) * a + (k + 1) * c * X ^ 2
      ↔ k * (k + 1) * a < c * X ^ 2 := by
  rw [← sub_pos, ← sub_pos (b := k * (k + 1) * a), unit_gap]

/-- Equal sharing is best for a square loss: shifting `d` from one unit to the other costs
`2 * d ^ 2`. -/
theorem equal_share (y d : R) : (y + d) ^ 2 + (y - d) ^ 2 = 2 * y ^ 2 + 2 * d ^ 2 := by ring

example : (1 * 2 * 2 * 2 + 1 * 1 * 3 ^ 2 : ℤ) < 1 * 1 * 2 * 2 + 2 * 1 * 3 ^ 2 := by norm_num

example : ¬((1 * 2 * 2 * 2 + 1 * 1 * 1 ^ 2 : ℤ) < 1 * 1 * 2 * 2 + 2 * 1 * 1 ^ 2) := by norm_num

/-- Kept fractions multiply, so savings compose by mult-add of the (negative) odds. -/
theorem stack (K s₁ s₂ : R) : (K - s₁) * (K - s₂) = K * K - (s₁ * K + s₂ * K - s₁ * s₂) := by
  ring

end Units

/-! ### Lead 4. Fans: power is proportional to flow cubed -/

section Fans

variable {R : Type*} [CommSemiring R]

/-- Layers `1, 3, 3, 1` of the fused cube `(1 + d) ^ 3 - 1 = 3 * d + 3 * d ^ 2 + d ^ 3`. -/
theorem cube_layers (d : R) : (1 + d) ^ 3 = 1 + 3 * d + 3 * d ^ 2 + d ^ 3 := by ring

/-- `k` fans each moving `q`: total power is `Q ^ 3 / k ^ 2` for total flow `Q = k * q`. -/
theorem spread (k q : R) : k * q ^ 3 * k ^ 2 = (k * q) ^ 3 := by ring

/-- `N + 1` fans at speed `N` versus `N` fans at speed `N + 1` (same total flow): the power
ratio is `N ^ 2 : (N + 1) ^ 2`. -/
theorem spare_fan (N : R) : ((N + 1) * N ^ 3) * (N + 1) ^ 2 = (N * (N + 1) ^ 3) * N ^ 2 := by
  ring

/-- Running the spare fan always helps. -/
theorem spare_fan_helps [PartialOrder R] [IsStrictOrderedRing R] {N : R} (h : 0 < N) :
    (N + 1) * N ^ 3 < N * (N + 1) ^ 3 := by
  have h1 : N ^ 2 < (N + 1) ^ 2 := pow_lt_pow_left₀ (lt_add_one N) h.le two_ne_zero
  have hp : 0 < N * (N + 1) := by positivity
  calc (N + 1) * N ^ 3 = N * (N + 1) * N ^ 2 := by ring
    _ < N * (N + 1) * (N + 1) ^ 2 := mul_lt_mul_of_pos_left h1 hp
    _ = N * (N + 1) ^ 3 := by ring

end Fans

end DataCenterFused
