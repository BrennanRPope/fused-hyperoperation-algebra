# Statistics, Lean part (for the library thread to fold in)

- `Statistics.lean`: 98 theorems, namespace `StatFused`, core Lean 4.37.0-rc1, no mathlib, no sorry,
  standard axioms only. Imports `FusedAlgebra.Base`, `FusedAlgebra.Cube`, `FusedAlgebra.Topics.E`.
  Copy to `FusedAlgebra/Topics/Statistics.lean`, add `FusedAlgebra.Topics.Statistics` to `build.sh`
  after `Topics.E` and to `FusedAlgebra.lean`. Compiled that way in a copy of the library on Nov 25.
- `CheckStatistics.lean`: axiom check for this module alone (writes build/statistics_theorems.tsv).
- `theorems.tsv`: the 98 theorems with their axioms.
- Reused from Base (not re-proved): circK lift/comm/assoc, mpow/wpow lift, foldCirc corner expansion,
  circInv. From E: binomial, rsum. From Cube: layer.
