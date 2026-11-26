# The dial's potential

Question: if the circle-operation dial k (a ∘₁ b = a + b + k·ab) behaves like a Higgs field, with a hump at k = 0 and dips at ±k0, what fixes its potential and its coupling? Approach (Brennan's): make 20 guesses, keep the closest, iterate.

Page: https://claude.ai/artifact/WpQzqeezDjmRmUcbgTA7RT

## Result

- 71 guesses over 4 rounds (20, 20, 15, 16), scored in bits against measured numbers, charged for every free parameter and for the size of the search, with held-out tests.
- The Higgs reading misses. Its best version, a dial that is flat at the Planck scale, predicts a Higgs mass of 128.7 ± 1.0 GeV against 125.2 (3.3σ off); it is the known 1995/2009 flat-potential idea, not a new one. Coleman-Weinberg, Veltman, unification and the round-number mass formulas are ruled out; the formulas that land close do so at chance rates (26 to 39% of random targets).
- One lead worth following: the dial as a dark-energy field that rolls. Rescaling forces its potential to be an exponential in ln k. Fitted to supernovae and BAO it beats a cosmological constant by Δχ² ≈ 4.85 (about 2σ), and it survives cross-validation (fit on one data set, score on the other). This is the same hint DESI reported in 2025. Slope best 0.82 (1σ 0.61 to 0.97); the conformal value √(2/3) = 0.816 sits inside, but 14 of 43 simple √(p/q) values do too, so that match means nothing yet, and the conformal reading would need screening to pass solar-system tests (it gives γ = 5/7).
- The dial cannot also be light's coupling: a light dial changing α would need a coupling below 2.4×10⁻⁵.

## Files

- `search.py`: the guesses, rounds, scoring, chance and held-out checks. Writes `results.json`, `search_output.txt`, `running.npz`.
- `sm_rge.py`: Standard Model running (2-loop gauge, leading 2-loop top Yukawa and Higgs quartic, Buttazzo et al. 2013 matching), checked against their NNLO λ(M_Pl).
- `quintessence.py`: rolling-field cosmology (exponential, hilltop, cosine) fitted to Pantheon+ and DESI DR2 BAO from `../dark-sector/data`.
- `lean/DialPotential.lean`: 24 theorems, no sorry, merged into the shared library as `Topics/DialPotential`.
- `build_page.py`, `page_template.html`, `index.html`: the page.

Rerun: `python3 -I search.py` (several minutes), then `python3 -I build_page.py`.

Polish (Nov 26): `polish.py` (numpy only, 40 s CPU; output in `polish_output.txt`) found that the hump shapes were stopped by a shooting-bracket bug in `quintessence.py` `history()` (the probe at 20× the starting field fails, so whole regions scored 1e9), not by the data. Fixed in `quintessence.py` (the top of the bracket shrinks until the probe runs) and `fit_grid` now polishes the shape parameter with Brent between the best grid point's neighbours; `search.py` now polishes the slope profile (Brent for the best slope, root-finds for the 1σ ends) and root-finds the scale where λ turns negative. Rerun of `search.py` with these: 184 s.

- Exponential slope: 0.80 (0.61 to 0.96, Δχ² 4.8) → 0.818 (0.612 to 0.974, Δχ² 4.85); 14 of 43 simple slopes inside (was 13).
- Cosine (dial on a circle): period 1.036 → 0.890, Δχ² 3.33 → 4.80, held-out +13.5 → +14.6. Still alive.
- Hilltop (the sombrero): width 1.055 → 0.497, Δχ² 0.14 → 4.53, but held-out +0.56 → −14.7: fitted to BAO alone it rolls fast (w0 ≈ −0.40) and misses the supernovae by 28 in χ². Now ruled out (it was alive). Its BAO-only fit sits near the grid's lower edge (0.43 against 0.4), so a wider grid could only make that held-out score worse.
- λ turns negative at 3.67×10¹⁰ → 3.62×10¹⁰ GeV.
- Page rebuilt from `results.json` with `build_page.py`.
