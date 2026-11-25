# Measure theory and probability through circle-operation

Page: index.html (self-contained). Checks: checks.py → checks_output.txt (38/38 pass). Lean: lean/MeasureProbability.lean (36 theorems, compiles clean), lean/README.md.

## Findings
- Sugeno's λ-measure (1974), μ(A∪B) = μA + μB + λμAμB for disjoint sets, is the dial with k = λ. Known.
- Its Möbius masses are the circle-operation cube's corners, λ^(|S|−1)·∏g_S; λ ≥ 0 gives belief functions, λ < 0 plausibilities (known: Banon 1981). Cube reading ours.
- Normalization ∏(1+λgᵢ) = 1+λ has one nonzero root above −1; λ > 0 iff Σg < 1 (Sugeno). Readings: λ < 0 = normalized failure union of independently firing atoms; λ > 0 = compound growth.
- log(1+λμ)/log(1+λ) is an ordinary probability P, and μ = [P]_q with q = 1+λ = k+1 (distortion known; q-number naming ours).
- Inclusion–exclusion becomes multiplicative (known); dual weight −λ/(1+λ) is the circle operation's inverse, q* = 1/q (known fact, reading ours).
- λ > −1 is "the wall −1/λ must lie outside [0,1]". On the landscape, decomposable measures (Weber 1984) normalize on finite spaces iff 1+pt+qt² has no root on [0,1]; edge = the statistics page's "or" line k+q = −1 up to the cusp (−2,1), then the border parabola, tangent there. Placement ours.
- E[X ∘ₖ Y] = EX ∘ₖ EY for independent X, Y; gap = k·Cov otherwise. Odds of an independent "or" combine by the circle operation ∘₁ (called mult-add before 1.0). Multiplicative martingales; Yor's formula as step-wise the circle operation; Itô drag = half the summed self-corners (Doléans-Dade 1970, Yor 1976).
- Dropped: Bernoulli convolutions (no structure beyond the generic product), odds of an "and" (fail; odds against work).

## Verdict
A small lead: the dial has a named home in measure theory (λ-measures). Next step worth taking: Choquet integrals against bent landscape measures (e.g. the tangent measure μ = tan(πP/4)) and whether they reproduce known aggregation rules. No new probability predicted.

Earlier pages cited, not redone: statistics (FwYdFUqfFh4DEAjHMqNEsv), information theory (SBsWmduwSNADpPVbeRDq2v), actuarial (AYgNMio88Y3VE9gtx11ipR), q (9AtexUZrzvsPfzSA7g91dL).
