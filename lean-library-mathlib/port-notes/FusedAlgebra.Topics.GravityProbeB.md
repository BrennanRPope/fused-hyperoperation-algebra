# Port notes: FusedAlgebra.Topics.GravityProbeB

- The file-local `✳` notation of the original is now the library's scoped notation.
- Identities are stated in any commutative ring or field; those that divide a numeral
  (`weak_split`, `ppn_geodetic`, `thirds`, `polar_average_*`) in any field of characteristic
  zero. The four theorems about measured Gravity Probe B numbers stay in `ℚ`.
- `gamma_sq` lost the hypothesis `1 - 3 * x ≠ 0`, which the proof never needed.
