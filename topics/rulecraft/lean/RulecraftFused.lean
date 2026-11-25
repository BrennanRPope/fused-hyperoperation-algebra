/-
  Rulecraft on the mult-add cube
  ==============================

  Part of Fused Hyperoperation Algebra. Companion to the page
  "Rulecraft, in fused algebra" and to ../checks.py (numbers) and
  ../figures.py (pictures).

  STATUS: UNCHECKED. Written in core Lean 4 only (no mathlib, no `sorry`),
  in the same style as ../../tesseract/ToricCharts.lean, which was compiled
  with Lean 4.37.0-rc1. Installing Lean was refused in the session that
  wrote this file, so it has not been compiled yet. Every identity below
  was checked symbolically with sympy in ../checks.py (section 0), so a
  failure would be a proof-script problem, not a false statement. Run it
  with `lean RulecraftFused.lean`.

  Mult-add: a ⊞ b = a + b + a·b = (1 + a)(1 + b) - 1.

  What it proves, in Rulecraft's words:

  * Scope matches spillover (the boundary theorem). Split the parties a
    decision touches into an inside list S and an outside list T. The
    whole effect is X_S ⊞ X_T = X_S + X_T + X_S·X_T. The cross term X_S·X_T
    is exactly the corners of the cube that touch both sides: the
    spillover a body with scope S cannot see. The corner expansion
    (`one_add_foldCirc_eq_corners`) shows the mult-add of a list is the sum,
    over every nonempty group of parties, of the product of their effects.

  * Order-free reforms. Mult-add is commutative and associative, so the
    mult-add of any reshuffled list is the same (`foldCirc_perm`). Percent
    changes to one quantity can never be path-dependent. The fused
    exponent rung is not commutative (`fexp_not_comm`), so path dependence
    needs a higher rung or more than one quantity.

  * The event horizon. -1 (losing everything) absorbs: (-1) ⊞ b = -1 for
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

  * Coupling weight. a ⊞_A b = a + b + A·a·b is ordinary mult-add seen
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

namespace CircleOp.Rulecraft

open Lean.Grind

set_option linter.unusedSectionVars false

/-! ## 1. Mult-add on any commutative ring -/

section Ring
variable {R : Type} [CommRing R]

/-- Mult-add: a ⊞ b = a + b + a·b. -/
def circ (a b : R) : R := a + b + a * b

theorem circ_lift (a b : R) : 1 + circ a b = (1 + a) * (1 + b) := by
  unfold circ; grind

theorem circ_comm (a b : R) : circ a b = circ b a := by
  unfold circ; grind

theorem circ_assoc (a b c : R) :
    circ (circ a b) c = circ a (circ b c) := by
  unfold circ; grind

theorem zero_circ (a : R) : circ 0 a = a := by
  unfold circ; grind

/-- The event horizon: -1 absorbs everything composed with it. -/
theorem neg_one_absorbs (b : R) : circ (-1) b = -1 := by
  unfold circ; grind

/-- The thermostat: undoing +d by -d leaves the corner square -d². -/
theorem thermostat_residue (d : R) : circ d (-d) = -(d * d) := by
  unfold circ; grind

/-- Weighted mult-add: the corner carries a coupling weight A. -/
def circK (A a b : R) : R := a + b + A * a * b

/-- Weighted mult-add is mult-add seen through x ↦ A·x. -/
theorem weighted_conj (A a b : R) : A * circK A a b = circ (A * a) (A * b) := by
  unfold circK circ; grind

/-- A = 0: no corner, plain addition (nothing spills over). -/
theorem circK_zero (a b : R) : circK 0 a b = a + b := by
  unfold circK; grind

/-- A = 1: plain mult-add (pure compounding). -/
theorem circK_one (a b : R) : circK 1 a b = circ a b := by
  unfold circK circ; grind

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

/-! ### The boundary theorem -/

/-- The mult-add of a list of effects: x₁ ⊞ x₂ ⊞ … ⊞ xₙ. -/
def foldCirc : List R → R
  | [] => 0
  | a :: l => circ a (foldCirc l)

/-- Scope matches spillover: the whole is the inside ⊞ the outside, and
    the cross term X_S·X_T is what crosses the boundary. -/
theorem foldCirc_append (S T : List R) :
    foldCirc (S ++ T) = circ (foldCirc S) (foldCirc T) := by
  induction S with
  | nil => simp only [List.nil_append, foldCirc]; rw [zero_circ]
  | cons a S ih =>
    simp only [List.cons_append, foldCirc]
    rw [ih, circ_assoc]

theorem boundary (S T : List R) :
    foldCirc (S ++ T) = foldCirc S + foldCirc T + foldCirc S * foldCirc T := by
  rw [foldCirc_append]; rfl

/-- Order-free reforms: any reshuffle of the same effects gives the same
    total. -/
theorem foldCirc_perm {l l' : List R} (h : List.Perm l l') : foldCirc l = foldCirc l' := by
  induction h with
  | nil => rfl
  | cons x _ ih => simp only [foldCirc]; rw [ih]
  | swap x y l =>
    simp only [foldCirc]
    unfold circ; grind
  | trans _ _ ih₁ ih₂ => exact ih₁.trans ih₂

/-! ### The corners: every nonempty group of parties -/

/-- All sub-lists (groups of parties), the corners of the cube. -/
def groups {α : Type} : List α → List (List α)
  | [] => [[]]
  | a :: l => groups l ++ (groups l).map (fun g => a :: g)

/-- Product of a group's effects: one corner's weight. -/
def prodL : List R → R
  | [] => 1
  | a :: l => a * prodL l

/-- Sum of a list. -/
def sumL : List R → R
  | [] => 0
  | a :: l => a + sumL l

theorem sumL_append (l m : List R) : sumL (l ++ m) = sumL l + sumL m := by
  induction l with
  | nil => simp only [List.nil_append, sumL]; grind
  | cons a l ih => simp only [List.cons_append, sumL]; rw [ih]; grind

theorem sumL_map_cons (a : R) (gs : List (List R)) :
    sumL ((gs.map (fun g => a :: g)).map prodL) = a * sumL (gs.map prodL) := by
  induction gs with
  | nil => simp only [List.map, sumL]; grind
  | cons g gs ih => simp only [List.map, sumL, prodL] at ih ⊢; rw [ih]; grind

/-- The corner expansion: 1 + (x₁ ⊞ … ⊞ xₙ) is the sum over every group
    of parties (the empty group gives the 1) of the product of their
    effects. Layer k of the cube is the groups of size k. -/
theorem one_add_foldCirc_eq_corners (l : List R) :
    1 + foldCirc l = sumL ((groups l).map prodL) := by
  induction l with
  | nil => simp only [foldCirc, groups, List.map, prodL, sumL]; grind
  | cons a l ih =>
    simp only [foldCirc, groups, List.map_append]
    rw [sumL_append, sumL_map_cons, ← ih]
    unfold circ; grind

end Ring

/-! ## 2. The way back needs division -/

section Field
variable {K : Type} [Field K]

/-- Everywhere but -1 there is an exact way back: -a/(1+a). -/
theorem way_back (a : K) (h : 1 + a ≠ 0) : circ a (-a / (1 + a)) = 0 := by
  unfold circ; grind

/-- At -1 there is none: nothing composes with -1 to give 0. -/
theorem no_way_back (b : K) : circ (-1) b ≠ 0 := by
  unfold circ; grind

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

/-- Rung 3: (1 + a)^(1 + b) - 1. -/
def fexp (a b : Nat) : Nat := (a + 1) ^ (b + 1) - 1

theorem fexp_not_comm : fexp 1 2 = 7 ∧ fexp 2 1 = 8 := by decide

/-! ## 5. The worked numbers -/

/-- Wisdom panel: expert (1, 9) and resident (9, 1). Own corners 18,
    cross corners 82, area 100. Two experts: area 36 = 10² - 8². -/
example : (1 + 9) * (9 + 1) = (1 * 9 + 9 * 1) + (1 * 1 + 9 * 9) := by decide
example : (1 + 1) * (9 + 9) = 10 * 10 - 8 * 8 := by decide

/-- Corner counting: of the 2^5 - 1 = 31 corners of five parties, a body
    holding 2 of them sees 3 and 21 straddle its edge (3 · 7). -/
example : (2 ^ 2 - 1) + (2 ^ 3 - 1) + (2 ^ 2 - 1) * (2 ^ 3 - 1) = 2 ^ 5 - 1 := by decide

end CircleOp.Rulecraft
