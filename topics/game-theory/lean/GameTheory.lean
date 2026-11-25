/-
  Game theory through the circle operation (called mult-add before 1.0)
  ============================

  Part of Fused Hyperoperation Algebra. Page and checks: ../ (game-theory).
  Core Lean 4 only (no mathlib), no `sorry`. Checked with Lean 4.37.0-rc1.

  The dial fold x₁ ◯ₖ x₂ ◯ₖ … ◯ₖ xₙ lifts to the product of the (1 + k·xᵢ). Proved:
  * Kelly: bets at fraction f compose by the dial with k = f; ruin is the dial's wall; Kelly's
    fraction makes the lifted wealth p/π on a win and q/(1 − π) on a loss.
  * Multiplicative weights (two experts): the total weight is 2·∏(1 + ε rₜ), each expert's weight
    is ∏(1 + ε gₜ), so the best expert's lifted gain is at most twice the algorithm's (the lifted
    regret bound, before logs).
  * Elo / Bradley–Terry: odds − 1 combines by the circle operation; the edge (O − 1)/(O + 1) composes like relativistic
    velocities; the odds are the edge's Doppler lift; win chances compose by p₁p₂/(p₁p₂ + q₁q₂).
  * The dial game v(S) = dial fold of the members' values: Sugeno's rule, marginal values,
    Shapley values for 2 and 3 players in closed form (Harsanyi dividends = cube corners, split
    equally), efficiency, Owen's multilinear extension = dial fold of t·x, Banzhaf = half
    participation, and who gets the larger bonus ratio.
-/
import FusedAlgebra.Base

namespace GameTheory

open CircleOp
open scoped CircleOp

set_option linter.unusedSectionVars false

/-! ## 1. The dial fold -/

section Ring
variable {R : Type} [Lean.Grind.CommRing R]

/-- x₁ ◯ₖ x₂ ◯ₖ … ◯ₖ xₙ. -/
def dfold (k : R) : List R → R
  | [] => 0
  | a :: l => circK k a (dfold k l)

theorem dfold_lift (k : R) (l : List R) :
    1 + k * dfold k l = prodL (l.map (fun x => 1 + k * x)) := by
  induction l with
  | nil => simp only [dfold, List.map, prodL]; grind
  | cons a l ih => simp only [dfold, List.map, prodL]; rw [circK_lift, ih]

theorem dfold_append (k : R) (S T : List R) :
    dfold k (S ++ T) = circK k (dfold k S) (dfold k T) := by
  induction S with
  | nil => simp only [List.nil_append, dfold]; unfold circK; grind
  | cons a S ih => simp only [List.cons_append, dfold]; rw [ih, circK_assoc]

theorem dfold_one (l : List R) : dfold 1 l = foldCirc l := by
  induction l with
  | nil => rfl
  | cons a l ih => simp only [dfold, foldCirc]; rw [ih, circK_one]

/-! ## 2. Kelly betting -/

/-- Two bets at fraction f: the returns combine by the dial with k = f. -/
theorem kelly_compose (f r₁ r₂ : R) :
    (1 + f * r₁) * (1 + f * r₂) = 1 + f * circK f r₁ r₂ := (circK_lift f r₁ r₂).symm

/-- Many bets: wealth is 1 + f·(r₁ ◯_f r₂ ◯_f … ◯_f rₙ). -/
theorem kelly_many (f : R) (rs : List R) :
    prodL (rs.map (fun r => 1 + f * r)) = 1 + f * dfold f rs := (dfold_lift f rs).symm

/-- In money terms the stakes' returns f·r combine by the circle operation plainly. -/
theorem kelly_stake (f r₁ r₂ : R) : f * circK f r₁ r₂ = (f * r₁) ◯ (f * r₂) := circK_scale f r₁ r₂

/-- Ruin is the wall: one outcome with f·r = −1 leaves nothing, whatever follows. -/
theorem kelly_ruin (f r : R) (h : f * r = -1) (rs : List R) : 1 + f * dfold f (r :: rs) = 0 := by
  rw [dfold_lift]; simp only [List.map, prodL]; rw [h]; grind

/-! ## 3. Multiplicative weights: one round -/

/-- One round on two experts: the new total is the old total plus ε times the algorithm's gain. -/
theorem hedge_step (ε u v g h : R) :
    u * (1 + ε * g) + v * (1 + ε * h) = (u + v) + ε * (u * g + v * h) := by grind

/-- The expert's weight after many rounds is its own dial fold. -/
theorem hedge_expert (ε : R) (gs : List R) :
    prodL (gs.map (fun g => 1 + ε * g)) = 1 + ε * dfold ε gs := (dfold_lift ε gs).symm

/-! ## 4. The dial game: v(S) = dial fold of the members' values -/

/-- Sugeno's rule: disjoint coalitions combine by the dial. The excess over adding is the corner. -/
theorem dial_game_union (k : R) (S T : List R) :
    dfold k (S ++ T) - dfold k S - dfold k T = k * dfold k S * dfold k T := by
  rw [dfold_append]; unfold circK; grind

/-- Joining coalition S is worth x·(1 + k·v(S)). -/
theorem marginal (k x : R) (S : List R) : dfold k (x :: S) - dfold k S = x * (1 + k * dfold k S) := by
  simp only [dfold]; unfold circK; grind

/-- The three-player worth expanded over the cube: corners weighted by k^(size − 1). -/
theorem corners3 (k x y z : R) :
    dfold k [x, y, z] = (x + y + z) + k * (x * y + x * z + y * z) + k * k * (x * y * z) := by
  simp only [dfold]; unfold circK; grind

/-- Owen's multilinear extension, two players: everyone present with chance t gives the
    dial fold of t·x. -/
theorem multilinear2 (k t x y : R) :
    (1 - t) * t * x + t * (1 - t) * y + t * t * circK k x y = circK k (t * x) (t * y) := by
  unfold circK; grind

/-- Three players. -/
theorem multilinear3 (k t x y z : R) :
    (1 - t) * (1 - t) * t * (x + y + z)
      + (1 - t) * t * t * (dfold k [x, y] + dfold k [x, z] + dfold k [y, z])
      + t * t * t * dfold k [x, y, z] = dfold k [t * x, t * y, t * z] := by
  simp only [dfold]; unfold circK; grind

end Ring

/-! ## 5. Fields: Kelly's fraction, Elo, Shapley and Banzhaf -/

section Field
/-- Numbers for this section: rationals (the 1/2, 1/3, 1/6 weights need characteristic 0). -/
abbrev K := Rat

/-- Kelly's fraction f = p − q/b makes the win lift p(1 + b) (= p/π, π = 1/(1 + b)). -/
theorem kelly_win_lift (p b : K) (hb : b ≠ 0) : 1 + b * (p - (1 - p) / b) = p * (1 + b) := by
  grind

/-- … and the loss lift q(1 + b)/b (= q/(1 − π)). -/
theorem kelly_loss_lift (p b : K) (hb : b ≠ 0) :
    1 - (p - (1 - p) / b) = (1 - p) * (1 + b) / b := by grind

/-- The bookmaker's chance π = 1/(1 + b): the win lift is your chance over theirs. -/
theorem kelly_ratio (p b : K) (hb : b ≠ 0) :
    1 + b * (p - (1 - p) / b) = p / (1 / (1 + b)) := by
  rw [kelly_win_lift p b hb]; grind

/-- First-order condition: d/df [p log(1 + bf) + q log(1 − f)] = pb/(1 + bf) − q/(1 − f) vanishes. -/
theorem kelly_foc (p b : K) (hb : b ≠ 0) (h1 : 1 + b ≠ 0) (hp : p ≠ 0) (hq : 1 - p ≠ 0) :
    p * b / (1 + b * (p - (1 - p) / b)) - (1 - p) / (1 - (p - (1 - p) / b)) = 0 := by
  rw [kelly_win_lift p b hb, kelly_loss_lift p b hb]
  have e1 : p * b / (p * (1 + b)) = b / (1 + b) := by grind
  have e2 : (1 - p) / ((1 - p) * (1 + b) / b) = b / (1 + b) := by grind
  rw [e1, e2]; grind

/-- Odds multiply along a chain, so odds − 1 combines by the circle operation. -/
theorem odds_circ (O₁ O₂ : K) : (O₁ - 1) ◯ (O₂ - 1) = O₁ * O₂ - 1 := by
  unfold circ; grind

/-- The edge 2P − 1 written in odds. -/
def edge (O : K) : K := (O - 1) / (O + 1)

/-- The win chance written in odds. -/
def prob (O : K) : K := O / (1 + O)

theorem edge_prob (O : K) (h : O + 1 ≠ 0) : edge O = 2 * prob O - 1 := by
  unfold edge prob
  have : 1 + O = O + 1 := by grind
  rw [this]; grind

/-- Edges compose like relativistic velocities: landscape point (k, q) = (0, −1). -/
theorem edge_compose (O₁ O₂ : K) (h1 : O₁ + 1 ≠ 0) (h2 : O₂ + 1 ≠ 0) (h12 : O₁ * O₂ + 1 ≠ 0) :
    edge (O₁ * O₂) = (edge O₁ + edge O₂) / (1 + edge O₁ * edge O₂) := by
  unfold edge
  have hd : 1 + (O₁ - 1) / (O₁ + 1) * ((O₂ - 1) / (O₂ + 1))
      = 2 * (O₁ * O₂ + 1) / ((O₁ + 1) * (O₂ + 1)) := by grind
  have hn : (O₁ - 1) / (O₁ + 1) + (O₂ - 1) / (O₂ + 1)
      = 2 * (O₁ * O₂ - 1) / ((O₁ + 1) * (O₂ + 1)) := by grind
  rw [hd, hn]
  have h3 : (O₁ + 1) * (O₂ + 1) ≠ 0 := by grind
  grind

/-- The odds are the edge's lift (1 + s)/(1 − s): relativity's Doppler factor squared. -/
theorem odds_is_lift (O : K) (h : O + 1 ≠ 0) : (1 + edge O) / (1 - edge O) = O := by
  unfold edge
  have e1 : 1 + (O - 1) / (O + 1) = 2 * O / (O + 1) := by grind
  have e2 : 1 - (O - 1) / (O + 1) = 2 / (O + 1) := by grind
  rw [e1, e2]; grind

/-- Win chances compose by p₁p₂/(p₁p₂ + q₁q₂). -/
theorem prob_compose (O₁ O₂ : K) (h1 : 1 + O₁ ≠ 0) (h2 : 1 + O₂ ≠ 0) (h12 : 1 + O₁ * O₂ ≠ 0) :
    prob (O₁ * O₂) = prob O₁ * prob O₂ /
      (prob O₁ * prob O₂ + (1 - prob O₁) * (1 - prob O₂)) := by
  unfold prob
  have e : O₁ / (1 + O₁) * (O₂ / (1 + O₂)) + (1 - O₁ / (1 + O₁)) * (1 - O₂ / (1 + O₂))
      = (1 + O₁ * O₂) / ((1 + O₁) * (1 + O₂)) := by grind
  rw [e]
  have h3 : (1 + O₁) * (1 + O₂) ≠ 0 := by grind
  grind

/-! ### Shapley values of the dial game -/

/-- Two players: each gets their own value plus half the corner. -/
def shapley2 (k x y : K) : K := ((dfold k [x] - 0) + (dfold k [x, y] - dfold k [y])) / 2

theorem shapley2_eq (k x y : K) : shapley2 k x y = x + k * x * y / 2 := by
  unfold shapley2; simp only [dfold]; unfold circK; grind

theorem shapley2_efficient (k x y : K) : shapley2 k x y + shapley2 k y x = dfold k [x, y] := by
  rw [shapley2_eq, shapley2_eq]; simp only [dfold]; unfold circK; grind

/-- k = −1, "some cause strikes": each cause keeps its own chance minus half the overlap
    (the actuarial UDD split). -/
theorem shapley2_union (x y : K) : shapley2 (-1) x y = x - x * y / 2 := by
  rw [shapley2_eq]; grind

/-- Three players, by the permutation average (orders weighted 2, 1, 1, 2 out of 6). -/
def shapley3 (k x y z : K) : K :=
  (2 * dfold k [x] + (dfold k [x, y] - dfold k [y]) + (dfold k [x, z] - dfold k [z])
    + 2 * (dfold k [x, y, z] - dfold k [y, z])) / 6

/-- Closed form: x·(1 + k·e₁/2 + k²·e₂/3) over the other two: each corner shared equally. -/
theorem shapley3_eq (k x y z : K) :
    shapley3 k x y z = x * (1 + k * (y + z) / 2 + k * k * (y * z) / 3) := by
  unfold shapley3; simp only [dfold]; unfold circK; grind

theorem shapley3_efficient (k x y z : K) :
    shapley3 k x y z + shapley3 k y x z + shapley3 k z x y = dfold k [x, y, z] := by
  rw [shapley3_eq, shapley3_eq, shapley3_eq]; simp only [dfold]; unfold circK; grind

/-- Corner by corner: singletons whole, each pair corner k·xy halved, the triple k²xyz in thirds. -/
theorem shapley3_corners (k x y z : K) :
    shapley3 k x y z = x + (k * x * y) / 2 + (k * x * z) / 2 + (k * k * x * y * z) / 3 := by
  rw [shapley3_eq]; grind

/-- Bonus ratios φᵢ/xᵢ differ by (x₂ − x₁)(k/2 + k²x₃/3): it is your partners' size that counts. -/
theorem bonus_gap (k x y z : K) :
    (1 + k * (y + z) / 2 + k * k * (y * z) / 3) - (1 + k * (x + z) / 2 + k * k * (x * z) / 3)
      = (y - x) * (k / 2 + k * k * z / 3) := by grind

/-- Banzhaf value, three players: the dial at half participation. -/
def banzhaf3 (k x y z : K) : K :=
  (dfold k [x] + (dfold k [x, y] - dfold k [y]) + (dfold k [x, z] - dfold k [z])
    + (dfold k [x, y, z] - dfold k [y, z])) / 4

theorem banzhaf3_eq (k x y z : K) : banzhaf3 k x y z = x * (1 + k * y / 2) * (1 + k * z / 2) := by
  unfold banzhaf3; simp only [dfold]; unfold circK; grind

end Field

/-- With synergy (k ≥ 0) the bigger player gets the smaller bonus ratio. -/
theorem bigger_smaller_bonus {k x y z : Rat} (hk : 0 ≤ k) (hz : 0 ≤ z) (hxy : y ≤ x) :
    (1 + k * (y + z) / 2 + k * k * (y * z) / 3) ≤ (1 + k * (x + z) / 2 + k * k * (x * z) / 3) := by
  have h := bonus_gap k x y z
  have hc : 0 ≤ k / 2 + k * k * z / 3 := by
    have := Rat.mul_nonneg (Rat.mul_nonneg hk hk) hz; grind
  have := Rat.mul_nonneg (by grind : (0:Rat) ≤ x - y) hc
  grind

/-! ## 6. Multiplicative weights over many rounds (two experts, Rat) -/

/-- Run the weights through a list of rounds (gain of expert 1, gain of expert 2). -/
def run (ε : Rat) : List (Rat × Rat) → Rat × Rat → Rat × Rat
  | [], w => w
  | (g, h) :: rest, (u, v) => run ε rest (u * (1 + ε * g), v * (1 + ε * h))

/-- The algorithm's growth: the product of (1 + ε·rₜ), rₜ its expected gain in round t. -/
def algFactor (ε : Rat) : List (Rat × Rat) → Rat × Rat → Rat
  | [], _ => 1
  | (g, h) :: rest, (u, v) =>
      (1 + ε * ((u * g + v * h) / (u + v))) * algFactor ε rest (u * (1 + ε * g), v * (1 + ε * h))

/-- Expert 1's own growth: the product of (1 + ε·gₜ). -/
def expertFactor (ε : Rat) : List (Rat × Rat) → Rat
  | [] => 1
  | (g, _) :: rest => (1 + ε * g) * expertFactor ε rest

/-- Every round keeps both lifts positive (|ε·g| < 1). -/
def Good (ε : Rat) (rounds : List (Rat × Rat)) : Prop :=
  ∀ p ∈ rounds, 0 < 1 + ε * p.1 ∧ 0 < 1 + ε * p.2

theorem run_total (ε : Rat) : ∀ (rounds : List (Rat × Rat)) (u v : Rat), 0 < u → 0 < v →
    Good ε rounds → (run ε rounds (u, v)).1 + (run ε rounds (u, v)).2 = (u + v) * algFactor ε rounds (u, v)
  | [], u, v, _, _, _ => by simp only [run, algFactor]; grind
  | (g, h) :: rest, u, v, hu, hv, hg => by
    have hgh := hg (g, h) (by simp)
    have hrest : Good ε rest := fun p hp => hg p (by simp [hp])
    have hu' : 0 < u * (1 + ε * g) := Rat.mul_pos hu hgh.1
    have hv' : 0 < v * (1 + ε * h) := Rat.mul_pos hv hgh.2
    simp only [run, algFactor]
    rw [run_total ε rest _ _ hu' hv' hrest]
    have hne : u + v ≠ 0 := by grind
    have : u * (1 + ε * g) + v * (1 + ε * h) = (u + v) * (1 + ε * ((u * g + v * h) / (u + v))) := by
      grind
    rw [this]; grind

theorem run_expert (ε : Rat) : ∀ (rounds : List (Rat × Rat)) (u v : Rat),
    (run ε rounds (u, v)).1 = u * expertFactor ε rounds
  | [], u, v => by simp only [run, expertFactor]; grind
  | (g, h) :: rest, u, v => by
    simp only [run, expertFactor]; rw [run_expert ε rest]; grind

theorem run_pos (ε : Rat) : ∀ (rounds : List (Rat × Rat)) (u v : Rat), 0 < u → 0 < v →
    Good ε rounds → 0 < (run ε rounds (u, v)).2
  | [], u, v, _, hv, _ => by simp only [run]; exact hv
  | (g, h) :: rest, u, v, hu, hv, hg => by
    have hgh := hg (g, h) (by simp)
    simp only [run]
    exact run_pos ε rest _ _ (Rat.mul_pos hu hgh.1) (Rat.mul_pos hv hgh.2)
      (fun p hp => hg p (by simp [hp]))

/-- The lifted regret bound: starting from equal weights, the best expert's growth is at most
    2 = n times the algorithm's. Taking logs gives Σ log(1 + ε g) ≤ Σ log(1 + ε r) + log 2. -/
theorem lifted_regret (ε : Rat) (rounds : List (Rat × Rat)) (hg : Good ε rounds) :
    expertFactor ε rounds ≤ 2 * algFactor ε rounds (1, 1) := by
  have t := run_total ε rounds 1 1 (by decide) (by decide) hg
  have e := run_expert ε rounds 1 1
  have p := run_pos ε rounds 1 1 (by decide) (by decide) hg
  grind

end GameTheory
