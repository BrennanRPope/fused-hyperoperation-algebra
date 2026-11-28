# The diamond hierarchical lattice on the landscape (Nov 25, 2026)

> **Names:** some words or symbols on this page (wall, Tsallis q) also mean something else in standard math, physics or computing, or inside this project. [NAMING-COLLISIONS.md](../../NAMING-COLLISIONS.md) says which meaning is in play.

Follow-up of the graph-theory lead: the diamond map t -> F_{Q-2,1-Q}(t^2, t^2) as Q slides along the "or" line k + d = -1.

Files
- index.html: the page (built by build.py from template.html; self-contained, figures as data: URIs, ~240 KB).
- dynamics.py: the map, basins, fixed points, nu(Q). figures.py -> figures/ (13 Julia-set frames, Q-plane picture, figdata.json).
- checks.py -> checks_output.txt: 70 checks, all pass (sympy identities, brute-force Potts sums on the 4- and 16-bond diamonds,
  level-2 chromatic polynomial, real-axis basin scan, Beraha number-field test, picture facts).
- lean/DiamondLattice.lean: 48 theorems, compiles clean; lean/README.md.

Findings
1. One diamond step = square in t (series), then square in the lift phi(t) = (1-t)/(1+(Q-1)t) = e^{-K} (parallel). Four walls:
   series 0, inf; parallel 1 and 1/(1-Q) (the colouring bond v = -1). Q only moves that one wall; the bend d = 1 - Q is its
   reciprocal. Ours (landscape reading); the map itself is Berker-Ostlund 1979, t = v/(Q+v) is Tsallis-Levy 1981.
2. p^2 - 4d = Q^2: walls always real, merge only at Q = 0, where the cold point turns parabolic (T'(1) = 1). The bend's sign
   never flips for real Q, so it carries no dynamical information. Ours (elementary).
3. Fixed points: t = 1/y with Q = (y-1)^2(y+1), e^K = y^2; slope at every one is 4/(1+t)^2 = (2y/(y+1))^2, so
   nu = ln2 / (2 ln(2/(1+t*))). Closed forms: Q = 1 golden (nu 1.6353), Q = 2 tribonacci (1.3383), Q = 32/27 t* = 3/5,
   Q = 3 t* = 1/2 (16/9, nu 1.2047), Q = 16 t* = 1/3, Q = 45 t* = 1/4. The fixed cubic and the closed forms at 1, 3, 16 are in
   Chang-Roeder-Shrock 2020 (q^2 + 2qv = v^3); the y-parametrisation and the slope formula we did not find written (elementary).
4. Fold at 32/27 (known, CRS 2020): Q(y) peaks at y = -1/3, slope -1; below it an attracting antiferro point exists.
5. Beraha test: cubic has a root in Q(4cos^2(pi/n)) only for n = 3, 6 (n <= 30); real switch points of the colouring point
   3, 1.6389, 1.41, 1.32, ... -> 32/27 (known, CRS) contain one Beraha value, Q = 3, explained by T being even and
   1/(1-Q) = -t* exactly when y = 2. nu(Q) is smooth through all Beraha values. No clustering.
6. Julia sets drawn for Q = 0, 0.6, 1, 32/27, 1.41, 1.64, 2, phi^2, 3, 3.247, 4, 8, 16 (t-plane, window +-3). Q = 2 has the
   extra symmetry t -> 1/t (walls harmonic). Pictures are standard (Derrida-De Seze-Itzykson 1983 in the x-plane).

Verdict: lead closed. The landscape gives a clean explanation (two multiplications, one moving wall) and a tidy slope formula,
but predicts nothing beyond Chang-Roeder-Shrock 2020 and Chio-Roeder 2021; Beraha values play no special role.
Sources checked: arXiv 1911.04012 (Chang, Roeder, Shrock, J. Math. Phys. 61, 073301, 2020); arXiv 1904.02195 (Chio, Roeder,
AIHPD 8(4), 2021).
