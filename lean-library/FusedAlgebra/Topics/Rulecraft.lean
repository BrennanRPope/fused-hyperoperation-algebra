/-
  Rulecraft on the circle-operation cube (called mult-add cube before 1.0)
  ==============================

  Part of Fused Hyperoperation Algebra. Companion to the page
  "Rulecraft, in fused algebra" and to ../checks.py (numbers) and
  ../figures.py (pictures).

  STATUS: compiled in the shared library with core Lean 4.37.0-rc1, no
  `sorry`. (The original file was written unchecked; one proof step needed
  a fix, and the generic circle-operation laws, the list fold, the corner expansion
  and the way back now come from FusedAlgebra.Base.)

  The circle operation: a ◯ b = a + b + a·b = (1 + a)(1 + b) - 1.

  What it proves, in Rulecraft's words:

  * Scope matches spillover (the boundary theorem). Split the parties a
    decision touches into an inside list S and an outside list T. The
    whole effect is X_S ◯ X_T = X_S + X_T + X_S·X_T. The cross term X_S·X_T
    is exactly the corners of the cube that touch both sides: the
    spillover a body with scope S cannot see. The corner expansion
    (`one_add_foldCirc_eq_corners`) shows the circle operation of a list is the sum,
    over every nonempty group of parties, of the product of their effects.

  * Order-free reforms. The circle operation is commutative and associative, so the
    circle operation of any reshuffled list is the same (`foldCirc_perm`). Percent
    changes to one quantity can never be path-dependent. The fused
    exponent rung is not commutative (`fexp_not_comm`), so path dependence
    needs a higher rung or more than one quantity.

  * The event horizon. -1 (losing everything) absorbs: (-1) ◯ b = -1 for
    every b, so nothing ever composes back to 0. Everywhere else there is
    an exact way back, -a/(1+a). Losing n parts of n+1 needs a gain of n
    to recover (`recovery`): 90% lost needs +900%.

  * The thermostat. Undoing a drift d by subtracting d leaves -d², the
    corner square of the drift (`thermostat_residue`).

  * Wisdom as an area, and mixed panels. A panel's area (Σe)(Σi) is the
    members' own corners plus the cross corners (`panel_grid`). For a fixed
    total, the balanced square is largest: (s+d)(s-d) = s² - d² ≤ s²
    (`balance`, `balance_le`). This is also why a tension held in balance
    keeps the most corner.

  * Coupling weight. a ◯[A] b = a + b + A·a·b is the ordinary circle operation seen
    through x ↦ A·x (`weighted_conj`). A > 0 is Benedict's high synergy,
    A < 0 a zero-sum corner, A = 0 no spillover at all.

  * Voice pricing with a free floor vote. Pricing total voice V as V² - 1
    (the floor vote is the 1) makes the next unit cost 2V + 1, growing
    with total voice (`fused_price_step`).

  Compounding by layers ((1 + x)ⁿ = Σ C(n,k) xᵏ, the debt-trap example) is
  proved in ../../topology/lean/CubeProofs.lean (`binomial`) and is not
  repeated here.

  All of this is classical algebra (Harsanyi 1959 for the corner
  decomposition of a group's value; Möbius inversion on the cube). What is
  new is only the reading on Rulecraft and the formal statements.
-/

import FusedAlgebra.Base
import FusedAlgebra.Rungs

namespace CircleOp.Rulecraft

open Lean.Grind

set_option linter.unusedSectionVars false

/-! ## 1. Readings of the circle operation on any commutative ring -/

section Ring
variable {R : Type} [CommRing R]

/-- The thermostat: undoing +d by -d leaves the corner square -d². -/
theorem thermostat_residue (d : R) : circ d (-d) = -(d * d) := by
  unfold circ; grind

/-- A panel of two: area = own corners + cross corners. -/
theorem panel_grid (e₁ e₂ i₁ i₂ : R) :
    (e₁ + e₂) * (i₁ + i₂) = (e₁ * i₁ + e₂ * i₂) + (e₁ * i₂ + e₂ * i₁) := by
  grind

/-- Balance: moving d away from the even split loses exactly d². -/
theorem balance (s d : R) : (s + d) * (s - d) = s * s - d * d := by
  grind

/-- Fused voice price: total voice V costs V² - 1, so the next unit of
    voice costs 2V + 1. -/
theorem fused_price_step (V : R) :
    ((V + 1) * (V + 1) - 1) - (V * V - 1) = 2 * V + 1 := by
  grind

/-- Plain and fused prices for k extra votes differ by the 2k strips. -/
theorem fused_price_expand (k : R) : (1 + k) * (1 + k) - 1 = k * k + 2 * k := by
  grind

/-! ### The boundary theorem (foldCirc, foldCirc_append, foldCirc_perm and the corner
  expansion one_add_foldCirc_eq_corners are proved once in FusedAlgebra.Base) -/

theorem boundary (S T : List R) :
    foldCirc (S ++ T) = foldCirc S + foldCirc T + foldCirc S * foldCirc T := by
  rw [foldCirc_append]; rfl

end Ring

/-! ## 2. The way back needs division (circ_circInv and
  no_circInv_at_wall are in Base) -/

section Field
variable {K : Type} [Field K]

/-- Losing n parts of n+1 needs a gain of n: 90% lost (n = 9) needs +900%. -/
theorem recovery (n : K) (h : n + 1 ≠ 0) : circ (-(n / (n + 1))) n = 0 := by
  unfold circ; grind

end Field

/-! ## 3. Balance on the integers: the even split is largest -/

theorem sq_nonneg' (d : Int) : 0 ≤ d * d := by
  rcases Int.le_total 0 d with h | h
  · exact Int.mul_nonneg h h
  · have h' : 0 ≤ -d := by omega
    have := Int.mul_nonneg h' h'
    rwa [Int.neg_mul_neg] at this

theorem balance_le (s d : Int) : (s + d) * (s - d) ≤ s * s := by
  have h1 : (s + d) * (s - d) = s * s - d * d := by
    rw [Int.add_mul, Int.mul_sub, Int.mul_sub, Int.mul_comm d s]; omega
  have h2 := sq_nonneg' d
  omega

/-! ## 4. The fused exponent rung does not commute -/

theorem fexp_not_comm_values : fexp 1 2 = 7 ∧ fexp 2 1 = 8 := by decide

/-! ## 5. The worked numbers -/

/-- Wisdom panel: expert (1, 9) and resident (9, 1). Own corners 18,
    cross corners 82, area 100. Two experts: area 36 = 10² - 8². -/
example : (1 + 9) * (9 + 1) = (1 * 9 + 9 * 1) + (1 * 1 + 9 * 9) := by decide
example : (1 + 1) * (9 + 9) = 10 * 10 - 8 * 8 := by decide

/-- Corner counting: of the 2^5 - 1 = 31 corners of five parties, a body
    holding 2 of them sees 3 and 21 straddle its edge (3 · 7). -/
example : (2 ^ 2 - 1) + (2 ^ 3 - 1) + (2 ^ 2 - 1) * (2 ^ 3 - 1) = 2 ^ 5 - 1 := by decide

end CircleOp.Rulecraft
