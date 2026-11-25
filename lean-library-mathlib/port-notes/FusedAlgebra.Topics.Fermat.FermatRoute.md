# FusedAlgebra.Topics.Fermat.FermatRoute

- The file does not prove Fermat's Last Theorem. It proves `NoFreyCurve → FermatLastTheorem`,
  where `NoFreyCurve` (the Wiles / Taylor–Wiles / Mazur / Ribet step) is a hypothesis.
- New, with no original: `noFreyCurve_of_flt`, the easy converse, which shows that `NoFreyCurve`
  is equivalent to `FermatLastTheorem` (so the hypothesis is exactly as strong as the theorem
  for primes `p ≥ 5`; the reduction is elementary bookkeeping, not progress on the deep step).
- Most helper lemmas of the original became mathlib one-liners; `cop_left` and the content of
  `normalize_core'` moved into `FreyPieces` (see that module's notes).
