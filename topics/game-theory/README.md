# Game theory through the circle operation (Nov 25, 2026)

## Files
- index.html: the page (self-contained, four live SVG figures: multiplicative weights, Kelly growth,
  Elo curve with rulers, the 3-player cube with Shapley bars).
- checks.py -> checks_output.txt: 31 checks (numpy/scipy simulations, exact Fractions), all pass.
- lean/GameTheory.lean: 36 theorems, compiles clean (see lean/README.md).

## Findings
1. Multiplicative weights (Littlestone-Warmuth 1994, Hedge 1997): each expert's weight is
   1 + ε·(g1 ✳ε g2 ✳ε …), the total is n·∏(1 + ε r_t), so the learner's lifted fortune is at least the
   leader's over n. Logs give the lifted regret bound; the corner (εg)² is the only cost of returning
   to plain gains, giving regret ≤ log n/ε + εT. Known; the dial reading is ours.
2. Kelly (1956): bets at fraction f compose by the dial with k = f; ruin is the dial's wall r = −1/f.
   At f* = p − q/b the lift is p/π on a win and q/(1 − π) on a loss (π the bookmaker's chance), growth =
   KL(p‖π). Known.
3. Elo/Bradley-Terry: odds − 1 combines by the circle operation ∘₁ (called mult-add before 1.0); the edge 2P − 1 composes like relativistic velocity
   (landscape point (0, −1)); rating gap = rapidity, odds = Doppler lift (1 + s)/(1 − s). The identity is
   the classical tanh law; the dictionary is ours as stated.
4. The dial game v(S) = (∏(1 + k x_i) − 1)/k (Sugeno's λ-measure, 1974): Harsanyi dividends are the cube's
   corners k^(|S|−1)∏x_S; Shapley = x_i Σ_m k^m e_m(others)/(m+1) = x_i ∫∏(1 + k t x_j)dt (exact for
   n ≤ 6); Owen's multilinear extension is the dial fold of t·x; Banzhaf = x_i ∏(1 + k x_j/2). Bonus
   ratio depends on partners: with synergy small players get the largest ratio, with overlap large
   players keep the most. k = −1 is the actuarial UDD split. Closed form not found written; follows from
   known pieces (Harsanyi 1963, Owen 1972/1975).
5. Dropped: best-of-n series is a binomial tail, not a fold.

## Verdict
Not a mathematical lead: every proof was already a log-of-a-product proof. A small practical lead as a
game-design tool: the dial game is a one-knob synergy model with a closed-form fair credit rule and a
cube picture of where every point comes from.
