# Remarks on the port of `FusedAlgebra.Topics.NSContinuum`

- The original imported `FusedAlgebra.Topics.Topology.CubeProofs` only for `binomial` and the
  hand-made `bsum`. The port proves `annulus_graded` directly from
  `MultAdd.Cube.mpow_eq_sum_choose`, so it imports only `FusedAlgebra.Cube`.
- `annulus_graded` drops the factor `1 ^ (k + 1)` from the summand (it is `1`).
- `stretch` is `Matrix.mulVec`; the planar vorticity `fun k ↦ if k = 2 then w else 0` is
  `Pi.single 2 w`.
