module

public import FusedAlgebra.Rungs

/-!
# Rulecraft on the mult-add cube

Readings of mult-add `a ✳ b = a + b + a * b = (1 + a) * (1 + b) - 1` in the words of Rulecraft.

* **Scope matches spillover** (`boundary`). Split the parties a decision touches into an inside
  list `S` and an outside list `T`. The whole effect is `X_S ✳ X_T = X_S + X_T + X_S * X_T`. The
  cross term is exactly the corners of the cube that touch both sides: the spillover a body with
  scope `S` cannot see. (`MultAdd.one_add_multAddList_eq_sum_sublists` is the corner expansion.)
* **Order-free reforms.** Mult-add is commutative and associative, so the mult-add of any
  reshuffled list is the same (`MultAdd.multAddList_perm`). The fused exponent is not commutative
  (`fexp_not_comm_values`), so path dependence needs a higher rung.
* **The event horizon.** `-1` absorbs (`MultAdd.neg_one_multAdd`); everywhere else there is an
  exact way back (`MultAdd.multAdd_multAddInv`). Losing `n` parts of `n + 1` needs a gain of `n`
  to recover (`recovery`).
* **The thermostat.** Undoing a drift `d` by subtracting `d` leaves `-d ^ 2`
  (`thermostat_residue`).
* **Wisdom as an area.** The area of a panel is the members' own corners plus the cross corners
  (`panel_grid`); for a fixed total the balanced square is largest (`balance`, `balance_le`).
* **Voice pricing with a free floor vote** (`fused_price_step`, `fused_price_expand`).

All of this is classical algebra (Harsanyi 1959 for the corner decomposition of the value of a
group; Möbius inversion on the cube); `balance` is mathlib's `sq_sub_sq`. What is recorded here
is the reading. The original stated `balance_le` over `ℤ`; here it holds in any linearly ordered
commutative ring.
-/

@[expose] public section

namespace MultAdd.Rulecraft

variable {R : Type*}

/-! ### 1. Readings of mult-add in any commutative ring -/

section CommRing

variable [CommRing R]

/-- The thermostat: undoing `+d` by `-d` leaves the corner square `-d ^ 2`. -/
theorem thermostat_residue (d : R) : d ✳ -d = -d ^ 2 := by
  unfold multAdd; ring

/-- A panel of two: area = own corners + cross corners. -/
theorem panel_grid (e₁ e₂ i₁ i₂ : R) :
    (e₁ + e₂) * (i₁ + i₂) = (e₁ * i₁ + e₂ * i₂) + (e₁ * i₂ + e₂ * i₁) := by
  ring

/-- Balance: moving `d` away from the even split loses exactly `d ^ 2`. This is `sq_sub_sq`. -/
theorem balance (s d : R) : (s + d) * (s - d) = s ^ 2 - d ^ 2 := (sq_sub_sq s d).symm

/-- Fused voice price: total voice `V` costs `V ^ 2 - 1`, so the next unit of voice costs
`2 * V + 1`. -/
theorem fused_price_step (V : R) : ((V + 1) ^ 2 - 1) - (V ^ 2 - 1) = 2 * V + 1 := by
  ring

/-- Plain and fused prices for `k` extra votes differ by the `2 * k` strips. -/
theorem fused_price_expand (k : R) : (1 + k) ^ 2 - 1 = k ^ 2 + 2 * k := by
  ring

end CommRing

/-- The boundary theorem: the effect on the whole is the inside effect, plus the outside effect,
plus the cross term that a body with scope `S` cannot see. -/
theorem boundary [NonUnitalSemiring R] (S T : List R) :
    multAddList (S ++ T) = multAddList S + multAddList T + multAddList S * multAddList T :=
  multAddList_append S T

/-! ### 2. The way back needs division -/

/-- Losing `n` parts of `n + 1` needs a gain of `n`: `90%` lost (`n = 9`) needs `+900%`. -/
theorem recovery {K : Type*} [Field K] (n : K) (h : n + 1 ≠ 0) : -(n / (n + 1)) ✳ n = 0 := by
  unfold multAdd; field_simp; ring

/-! ### 3. Balance: the even split is largest -/

theorem balance_le [CommRing R] [LinearOrder R] [IsStrictOrderedRing R] (s d : R) :
    (s + d) * (s - d) ≤ s ^ 2 := by
  rw [balance]; exact sub_le_self _ (sq_nonneg d)

/-! ### 4. The fused exponent rung does not commute -/

theorem fexp_not_comm_values : fexp 1 2 = 7 ∧ fexp 2 1 = 8 := by decide

/-! ### 5. The worked numbers -/

/-- Wisdom panel: expert `(1, 9)` and resident `(9, 1)`. Own corners `18`, cross corners `82`,
area `100`. -/
example : (1 + 9) * (9 + 1) = (1 * 9 + 9 * 1) + (1 * 1 + 9 * 9) := by decide

/-- Two experts: area `36 = 10 ^ 2 - 8 ^ 2`. -/
example : (1 + 1) * (9 + 9) = 10 * 10 - 8 * 8 := by decide

/-- Corner counting: of the `2 ^ 5 - 1 = 31` corners of five parties, a body holding `2` of them
sees `3`, and `21` straddle its edge (`3 * 7`). -/
example : (2 ^ 2 - 1) ✳ (2 ^ 3 - 1) = 2 ^ 5 - 1 := by decide

end MultAdd.Rulecraft
