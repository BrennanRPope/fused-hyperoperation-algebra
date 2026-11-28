# Cusp size: with the cusp fitted out, finite chains give Knauf's limit on the line  (DRAFT)

**Plain result.** Knauf's chain sums only converge to ζ(β−1)/ζ(β) for Re β > 2. On the line
Re β = 3/2 a single chain of length 33 still misses that function by about 0.025. Fitting out the
cusp tail (cusp-tail: Z_n = F + A·n^(1−β) + B·n^(−β), over lengths 17 to 33) gives an F that
matches ζ(β−1)/ζ(β) to **0.08%** (median relative error) everywhere from height 5 to 45. That's
about 30 times closer than the single chain, with the dips at every zeta zero in the right place
(picture A). So the chain's finite sums already contain the continued function on the line. What
stands between them is mostly the cusp.

**The cusp's size A(β)** is not zero at the zeta zeros and doesn't match any simple ζ formula tried
(1/ζ(β), ζ(β−1), 1/(β−1) and products; picture B). Working it out from the chain's fractions is
still open: the run structure gives a renewal equation whose kernel mixes the cusp with the bulk.

## Our fuller toolkit on the cusp

| Tool | What it says here | Status |
|---|---|---|
| Parabolic fingerprints (parabolic-fingerprints) | The Farey map's cusp branch x/(1−x) is the Möbius map in that table: fingerprint (0, 0), the class that sits in a continuous flow. The rung step e^x − 1 is (1/3, 64.15). Prediction: the tail is a clean series in powers of n with no log terms. | Tested; inconclusive at these lengths (the extra log columns are not stably small) |
| Homothety dial k on the chain step, (u, v) → (u ∘ₖ v, v) = (u + v + kuv, v) | k = 0 is the chain; any k ≠ 0 makes a run grow geometrically, so the cusp turns from parabolic to hyperbolic and the n^(1−β) tail becomes geometric. A rewind from k > 0 to k = 0 could remove the cusp by construction. | Not tested |
| Roundness dial p | Temperature is p (knauf-field); the line Re β = 3/2 is its ℓ² edge. | Known from earlier pieces |
| Order dial q | q-deforming the chain would make the cusp's run a q-number [m]_q, and the ½-area came from q → 1. | Not tested |

## Grades

| | worth following? | could reach a proof? |
|---|---|---|
| Cusp size | **Yes.** It shows the finite chains already hold the continued limit on the line once the cusp is removed, and it narrows the open part to one function, A(β), plus a small remainder. | Not as it stands. A proof would need A(β) and the remainder controlled at every height; this is numerical, lengths up to 33, heights up to 45. |

## Pictures

- `figures/figA_limit_recovered.png`: A1, the fitted F on top of |ζ(β−1)/ζ(β)| along the line; A2, the miss for one chain against the cusp fit.
- `figures/figB_cusp_size.png`: |A| along the line.

## Known vs ours

Known: Knauf's chain and its limit for β > 2; the Farey map's parabolic point and its role in the
β = 2 transition (Prellberg, Kleban–Özlük, Fiala–Kleban–Özlük); the Möbius map's flow embedding.
Ours (not found in a quick search): recovering ζ(β−1)/ζ(β) on Re β = 3/2 from finite chains by
fitting out the cusp tail; the toolkit table.

## Files

`scan.py` → `scan.json` (fit along the line), `compare.py` → `compare.json`, `shape.py` →
`shape.json` (A against ζ formulas; log-term test), `figs.py`, `checks.py` (6 checks). Lengths above
20 come from `../length-smoothing/stream.py` (rerun it to 33 first).
