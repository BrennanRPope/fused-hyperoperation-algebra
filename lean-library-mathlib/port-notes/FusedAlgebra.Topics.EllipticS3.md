# FusedAlgebra.Topics.EllipticS3

- New: `gen`, `trace_gen`, `det_gen`, `quatMatrix`, `trace_quatMatrix` (the original only said
  "trace U = 2w" in a comment).
- The moments are only checked for `n < 5`, as in the original (which listed `n = 1..4`). The
  general statements (`walks (2n) 0 = centralBinom n`, `dyck (2n) 0 = catalan n`) are true and
  could be linked to mathlib's `DyckWord`, but were not in the original and are not proved here.
