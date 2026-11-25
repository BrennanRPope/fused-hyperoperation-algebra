# Point processes on the dial (Nov 25, 2026)

Finding: the α-determinantal processes (Shirai–Takahashi 2003) are the dial with k = α: k = −1 fermions, 0 Poisson,
+1 bosons. On a region with kernel eigenvalues μᵢ, P(no point) = (1 + k·(μ₁ ∘ₖ μ₂ ∘ₖ …))^(−1/k) = exp(−Σ ℓₖ(μᵢ)),
ℓₖ the dial's logarithm. Existence: the wall −1/k must stay off the spectrum, and for k < 0 it must be a whole number.

## Files
- index.html: the page (self-contained; SVG drawn by script, one slider). Built by `python3 build_page.py page_template.html`,
  which embeds results.json and the check count.
- checks.py -> checks_output.txt: 101 of 101 checks pass (numpy, scipy). Exact count laws by convolving independent modes;
  Monte Carlo at k = −1, −1/2, −1/3 (superposed fermion copies, spectral sampler), 0 (Poisson), 1/2, 1, 2 (Cox processes
  with m squared real Gaussian fields). Writes results.json (simulation summaries, sample pictures).
- lean/DialPointProcesses.lean: 40 theorems, compiles clean, no sorry, standard axioms. See lean/README.md.

## Findings
1. Void = dial fold: 1 + k·fold = det(I + kK); −log P(empty) = Σ ℓₖ(μ), so voids look Poisson with mean Σℓₖ(μ). Matches
   simulation at 7 dial values (whole count law, chi-square). Known formula (Shirai–Takahashi); dial reading ours.
2. Existence: one mode is a real count iff k ≥ 0, or k = −1/m with μ ≤ m (Shirai–Takahashi; Maunoury 2016 Cor. 7).
   Picture: blue stalks at k = −1/m ending on the wall curve μ = −1/k. Same "wall outside the range of values" rule as
   Sugeno's λ-measures (measure page), plus a whole-number lattice, and the process may touch the wall (sure point,
   zero variance; rigid count checked). For k > 0 the full process needs an extra kernel condition (Maunoury), always met
   at k = 2/m (real kernels).
3. Variance: each mode adds μ(1 + kμ), value times lift, zero at 0 and at the wall; Var = tr K + k tr K², a straight Fano
   line in k (fit to 7 simulations: slope 0.579 vs 0.577, intercept 1.001). Known (Shirai–Takahashi; Morris 1982).
   Superposition: m copies of (k, K) = (k/m, mK), from c·(a ∘ₖ b) = (ca) ∘_{k/c} (cb). Known; identity in Lean.
4. Bend: tilting a mode by t^N, in the ruler (mean − μ)/Var, composes by the landscape lens law at p = 1 + 2kμ,
   d = kμ(1 + kμ): always on p² − 4d = 1. Rare mode = the circle operation ∘₁ (called mult-add before 1.0) (1,0), full fermion mode = union of failures (−1,0),
   half-filled = velocity addition with limit 2 (0, −1/4); shape (2F − 1)², F = Fano. Placement ours; the bend is pinned by
   the filling, so it gives no new processes.
5. Fractional rungs: nothing genuine (shape matching only); dropped.

## Verdict
No further lead. The dial reading is exact and gives clean pictures, but every count statement is already known
(Shirai–Takahashi, HKPV, Maunoury, Morris), and the two open questions (bend, fractional rungs) produce no new family or
prediction.
