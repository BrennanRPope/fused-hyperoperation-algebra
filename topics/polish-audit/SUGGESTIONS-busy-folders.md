# DRAFT · Polish suggestions for folders other threads own (Nov 26)

These folders were being changed by other threads, so nothing in them was edited. Each line is a suggested change.

- **qsharp-functors**: already has a polish step (checks P1–P3: exact T-counts in fractions of pi, one Newton step on the
  Euler angles, Richardson extrapolation of the Trotter ratio, and an exact rational k-series replacing the T4 fit). Nothing
  more suggested.
- **shape-dial-survey/new_laws.py**: already polishes its grid maxima with a tight Nelder–Mead (`polish_2d`). Nothing more.
- **shape-dial-survey/checks.py**, B2 (`homog_residual`, line ~124) and D2 (`isometry_residual`, line ~226): both call
  `optimize.least_squares` at default tolerances (ftol = xtol = 1e-8), then compare a "mismatch ~ 0" against nonzero ones.
  Suggest `least_squares(f, x0, method="trf", xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=4000)`, or a second
  `least_squares(..., method="lm")` started from `sol.x`. Cost: a few more evaluations per start (well under a second).
  Effect: the "yes" cases should drop from about 1e-8 to about 1e-14, widening the gap to the "no" cases; verdicts should
  not change.
- **shape-dial-survey/checks.py**, Lee–Yang check (line ~440): `X` is picked by argmax on a 400-point grid but only sets
  a contour edge, so no polish is needed.
- **exceptional-shapes, hypersphere, bridge-survey**: no numerical searches or fits found.
