# FusedAlgebra.Topics.Landscape

Theorems: 27. GENERAL 27 (ℤ → any CommRing; `channel_floor` → any linearly ordered commutative ring, with the same three hypotheses). No findings, no low-confidence items.

- `num` = `wAdd p` and `den q a b = 1 - q*a*b` are unchanged. The point names `multAdd_point` and `velocity_point` plug in the same (u, v) values.

Route changes: `lens_lift` is now derived from `lens_factor_u` and `lens_factor_v`, which are his own lemmas, and `parabola_den` is the term `lens_factor_u u u`. These are not meaningful changes.

New port theorems with no original: 2 (`num_def`, `den_def`).

Rung index note: no theorem in this module uses `MultAdd.rung`. Where the docstrings say "rung 1" or "rung s" they mean the fractional height s (s = 0 addition, s = 1 mult-add), the originals' own convention in this topic; that differs from Rungs.lean / `MultAdd.rung`, where mult-add is index 2. Same in original and port, so not a port finding.
