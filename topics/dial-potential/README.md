# The dial's potential

Question: if the mult-add dial k (a ⊞ b = a + b + k·ab) behaves like a Higgs field, with a hump at k = 0 and dips at ±k0, what fixes its potential and its coupling? Approach (Brennan's): make 20 guesses, keep the closest, iterate.

Page: https://claude.ai/artifact/WpQzqeezDjmRmUcbgTA7RT

## Result

- 71 guesses over 4 rounds (20, 20, 15, 16), scored in bits against measured numbers, charged for every free parameter and for the size of the search, with held-out tests.
- The Higgs reading misses. Its best version, a dial that is flat at the Planck scale, predicts a Higgs mass of 128.7 ± 1.0 GeV against 125.2 (3.3σ off); it is the known 1995/2009 flat-potential idea, not a new one. Coleman-Weinberg, Veltman, unification and the round-number mass formulas are ruled out; the formulas that land close do so at chance rates (26 to 39% of random targets).
- One lead worth following: the dial as a dark-energy field that rolls. Rescaling forces its potential to be an exponential in ln k. Fitted to supernovae and BAO it beats a cosmological constant by Δχ² ≈ 4.8 (about 2σ), and it survives cross-validation (fit on one data set, score on the other). This is the same hint DESI reported in 2025. Slope best 0.80 (1σ 0.61 to 0.96); the conformal value √(2/3) = 0.816 sits inside, but 13 of 43 simple √(p/q) values do too, so that match means nothing yet, and the conformal reading would need screening to pass solar-system tests (it gives γ = 5/7).
- The dial cannot also be light's coupling: a light dial changing α would need a coupling below 2.4×10⁻⁵.

## Files

- `search.py`: the guesses, rounds, scoring, chance and held-out checks. Writes `results.json`, `search_output.txt`, `running.npz`.
- `sm_rge.py`: Standard Model running (2-loop gauge, leading 2-loop top Yukawa and Higgs quartic, Buttazzo et al. 2013 matching), checked against their NNLO λ(M_Pl).
- `quintessence.py`: rolling-field cosmology (exponential, hilltop, cosine) fitted to Pantheon+ and DESI DR2 BAO from `../dark-sector/data`.
- `lean/DialPotential.lean`: 24 theorems, no sorry, merged into the shared library as `Topics/DialPotential`.
- `build_page.py`, `page_template.html`, `index.html`: the page.

Rerun: `python3 -I search.py` (several minutes), then `python3 -I build_page.py`.
