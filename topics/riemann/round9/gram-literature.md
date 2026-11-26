# Who has studied where and why Gram's law fails (literature check, DRAFT)

*Nov 26 2026. Part of round 9 (the string probe). Sources read are listed at the end; anything I could not read myself is marked.*

## The answer in plain words

Yes, this has been studied for a hundred years, and the "where" and the "why" are both well understood.

- **Where:** about 73% of Gram intervals obey Gram's law at the heights computed so far. The share drifts down very slowly, and random-matrix models predict it settles near 66%.
- **Why:** the number of zeros in a Gram interval is exactly 1 plus the jump in S(t), the wobble term in the zero count. Gram's law fails whenever S(t) jumps across an interval. Rosser's rule fails when S(t) swings by 2 inside a Gram block. S(t) is unbounded, so failures can never stop.
- **Still open:** nobody has proved that Gram's law, or Rosser's rule, holds even infinitely often. Failures are proved; successes are not.

## Who proved what

| Year | Who | What |
|---|---|---|
| 1903 | Gram | Introduced Gram points; found the first 15 zeros and saw one zero per interval. He said he did not expect the pattern to hold in general. |
| 1925 | Hutchinson | Named it "Gram's law"; computed 138 zeros and found the first failure: [g125, g126) has no zero and [g126, g127) has two. (Korolev's survey indexes the same event as G127.) |
| 1934 | Titchmarsh | The weak Gram law (an odd number of zeros per interval) holds infinitely often. |
| 1935 | Titchmarsh | Gram's law fails infinitely often. With Comrie he found 43 exceptions among the first 1041 zeros. |
| 1946 | Selberg (stated without proof) | Gram's law and the weak law each fail for a positive proportion, and the weak law holds for a positive proportion. |
| 1969 | Rosser, Yohe, Schoenfeld | Rosser's rule: a Gram block of length l usually holds l zeros. They doubted it would hold forever. |
| 1970 | Lehman | The weak Rosser rule, and so Rosser's rule, fails infinitely often (using Selberg's result that S(t) is unbounded below). He conjectured the first failure lies beyond 10^7. |
| 1979 | Brent | First failure of Rosser's rule, at Gram point 13,999,525, confirming Lehman. 75 million zeros checked. |
| 1982–86 | Brent, van de Lune, te Riele, Winter | Up to 1.5 billion zeros; tables of Gram block types and Rosser exceptions. |
| 1987 | Fujii | Random-matrix (GUE) model for Gram intervals ("Fujii's conjecture"); also stated positive-proportion results without proof. |
| ~2004 | Gourdon | 10^13 zeros; lists Rosser violations through the 320,624,341st. |
| 2011 (preprint 2008) | Trudgian, Acta Arith. 148 | First full proofs: Gram's law and the weak law fail for a positive proportion; the weak law holds for a positive proportion; the weak Rosser rule, and so Rosser's rule, fails for a positive proportion. |
| 2010–13 | Korolev, Izvestiya | Short-interval versions of these results; proves Selberg's formulas; studies adjacent zeros that fail together; survey with a table of what is open. |
| 2019–20 | Hanga and Hughes (York) | Probabilistic models for Gram's law; an SU(N) correction to Fujii's model that explains the slow convergence. |
| 2026 | Jerby (arXiv 2610.21952) | Good and bad Gram points sit in different regions when prime sums are sampled at Gram points. An observation, not a predictor. |

## The numbers

Gram intervals with k zeros at the first 10^9 Gram points (Brent, van de Lune, te Riele, Winter, as tabulated by Hanga and Hughes), against the GUE limit from Fujii's model:

| zeros in the interval | observed | GUE limit |
|---|---|---|
| 0 | 13.7% | 17.0% |
| 1 (Gram's law holds) | 72.8% | 66.1% |
| 2 | 13.4% | 16.6% |
| 3 | 0.18% | — |

Hanga and Hughes quote 72.61% for the full 1.5 billion intervals.

## Why it fails (the mechanism)

At a Gram point, θ(g_n) = nπ, so the zero count is N(g_n) = n + 1 + S(g_n) exactly. The zeros in [g_n, g_{n+1}) number 1 + S(g_{n+1}) − S(g_n). So:

- Gram's law holds exactly when S(t) is the same at both ends of the interval.
- A Gram block holds as many zeros as its length when S(t) stays small across it; the literature (Lehman's analysis, restated in arXiv 2512.23064) ties Rosser failures to S(t) reaching size 2 inside the block.
- Selberg's central limit theorem says S(t) is Gaussian with variance about (1/2π²) log log t. It grows without bound but extremely slowly, which is why the first Rosser failure waits until 1.4 × 10^7 and why the Gram share drifts down so slowly toward the GUE value.

## Next-zero predictors

I found no published rule that beats Gram's law or Rosser's rule at predicting the next zero in a fair held-out test. A 2018 student project (Kampe and Vysogorets) fit neural networks to zero sequences; a later paper notes that next-term prediction can become trivial through overfitting. A 2025 machine-learning paper was revised after a data-leakage concern. Rigorous zero verification does not use a predictor at all: it uses Turing's method (1953), which works because the average of S(t), not S(t) itself, is bounded.

## What is open, and where our work could add something

Open (per Korolev's survey and Hanga's thesis):
- whether Gram's law holds infinitely often, or for a positive proportion;
- whether Rosser's rule, or even the weak Rosser rule, holds infinitely often;
- Fujii's conjecture (the 66% limit);
- whether a zero ever falls exactly at a Gram point.

Where round 9 could add something, plainly: a held-out test of predictors (Gram, Rosser, Riemann–Siegel corrections, and our candidate strings) on the same zeros would be a clean comparison, but it would be a numerical study, not new mathematics. A rule that mattered for a proof would have to come with an argument that it never fails at any height. That is what the open problems above ask for, and no one has one.

## Corrections to what was said on the call

- Hutchinson 1925, Titchmarsh (1935) and Lehman (1970) were right as stated. The first Rosser failure was found by Brent (1979).
- Trudgian proved the failures (and the weak Gram law's successes) for a positive proportion. He did not prove that Gram's law or Rosser's rule holds for a positive proportion; that is open. Selberg (1946) and Fujii (1987) had stated the failure results without proof.
- Hiary and Odlyzko's 2012 paper is about moments of ζ against random matrices, not Gram points. Odlyzko's link to Gram points is his unpublished conjecture that zeros become unrelated to Gram points at large heights (cited by Arias de Reyna, arXiv 1403.7126).

## Sources

- Trudgian, Gram's law fails a positive proportion of the time, arXiv 0811.0883; Acta Arith. 148 (2011) 225–256 (abstract page only; the theorem statements are from Korolev's and Hanga's summaries).
- Korolev, Gram's law and the argument of the Riemann zeta function, arXiv 1106.0516 (survey; source for Selberg, Fujii, Titchmarsh, and the open-problem table).
- Hanga and Hughes, Probabilistic models for Gram's law, arXiv 1911.03190.
- Hanga, Random matrix models for Gram's law, PhD thesis, York (2020).
- OEIS A216700 (Rosser violations; Lehman, Brent, Gourdon, Trudgian references).
- van de Lune, te Riele, Winter, Math. Comp. 46 (1986) (abstract only).
- Brent, Math. Comp. 33 (1979) (abstract only).
- Amberger, arXiv 2512.23064 (Gram block definition after Lehman).
- Jerby, arXiv 2610.21952 (abstract).
- Arias de Reyna, arXiv 1403.7126 (Odlyzko's conjecture).
- Kampe and Vysogorets, 2018 REU preprint (abstract-level).
- Not read: Odlyzko's own site (blocked here), Gourdon's report, the full Brent and van de Lune papers.
