# Lee–Yang zeros on the dial (Nov 25, 2026)

> **Names:** some words or symbols on this page (⊕, rung, wall) also mean something else in standard math, physics or computing, or inside this project. [NAMING-COLLISIONS.md](../../NAMING-COLLISIONS.md) says which meaning is in play.

Page: https://claude.ai/artifact/UQQzmkzekxe2hL4m8F8NQL

- dialgas.py: dial gas Z = sum_S z^|S| k^e(S) (one edge = 1 + x (+)_k y); exact counts for rings and tori by transfer matrix.
- potts_strips.py: chromatic polynomials of cyclic square-lattice strips (dial k = Q - 1).
- checks.py -> checks_output.txt: 190 checks, all pass. figures.py -> figures/. page/ is the published page.
- lean/LeeYang.lean: 39 theorems, core Lean 4.37.0-rc1, no sorry, standard axioms. Intended FusedAlgebra/Topics/LeeYang.lean,
  sent to the library thread to merge (compiled here against lean-library/build).
- Pickles: torus_counts.pkl, torus_zeros.pkl, fisher_zeros.pkl, potts_zeros.pkl.

Findings: ring zeros fill s = zeta+1/zeta from -2 (wall) to 2-4/k; odd rings vanish on the wall exactly; only k=1 and k=0 change the shape;
2D transitions at k = 3 +- 2 sqrt2. Potts reading: k + 1 = number of colours; k=-1,0,1 walls (Beraha), k=-2 never a zero.
Next test: do landscape laws cluster at Beraha dial values k = 1 + 2cos(2pi/n)?

## Cluster test (Nov 25)
Question: do landscape laws sit at the colour-count predictions k_n = 1 + 2cos(2π/n) (n ≥ 5: 1.618, 2, 2.247, …, → 3)?
- `cluster_test.py` → `cluster/cluster_output.txt`, census `cluster/landscape_census.csv` (85 fixed placements from landscape, thermodynamics, simple-laws, statistics; 11 parameter families listed, not scored). Figure `cluster_figure.py` → `figures/f4_cluster.png`.
- Result: 3 hits, all at k = 2 (leverage, Euler cusp, rung 2); 0 at the seven non-integer predictions. Random-set null p = 0.099; simple-number null (denominator ≤ 4) matches it 80% of the time. Whole dial: integers −2..2 beat Beraha n = 2..6 (80 vs 76 laws); k = −2 (ring never vanishes) holds 4 laws. Units a ↦ λa move λk to k, so positive dial values on d = 0 are unit choices.
- Lean: `lean/LeeYangCluster.lean`, 10 theorems (Beraha minimal polynomials have no integer root; circK_units; leverage_is_compounding; units_bend). Merged into the library (1626 theorems, 55 modules).
- Verdict: lead closed.

## Fill-then-spill test (Nov 25)
Question: do phase changes fill like a cup then spill (latent-heat plateau), the way the rungs might?
- Rungs (`fillspill_rungs.py`): rung curve (s, s(s−1)/4) has shape number s (constant pace); flow E^s(x) speeds up monotonically, rung-locked wobble ~1e-5 vs ±0.1 for other valid interpolations. No fill/spill.
- Border (`fillspill_border.py`): repeating a ⊕ x settles in a wall on the two-walls side; past the border it lingers by the ghost of the double wall, then spills; linger ~ Δ^(−1/2), Δ = 4d − k². Generic fold null gives the same exponent.
- Heat (`fillspill_heat.py`): Potts latent heat (Baxter) zero for k ≤ 3, 0.696 at k = 9 (MC 0.695); dial gas coexistence plateau, finite-size onset 3.06/3.41/3.69 (L = 4/5/6) toward k_c = 5.83. Latent heat → 0 at threshold while border linger → ∞; edge laws stay at (k, 0), never cross the border.
- Figure `fillspill_figure.py` → `figures/f5_fillspill.png`. Lean `lean/LeeYangSpill.lean`, 10 theorems, merged into the library (1682 theorems, 60 modules).
- Verdict: not a lead; same picture, different mechanisms.
