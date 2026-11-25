# Differential geometry through mult-add (Nov 25, 2026)

Files
- index.html: the page (conic-and-chart figure, equal-distance ticks, mult-add in the upper half plane,
  landscape shaded by curvature, Klein and Poincaré disks).
- checks.py → checks_output.txt: 34 of 34 checks pass (sympy, mpmath, scipy, numpy).
- lean/DiffGeometry.lean: 32 theorems, compiled clean (see lean/README.md).

Findings
- Turning the conic t² + Kx² = 1 and reading the central chart r = x/t gives (r1 + r2)/(1 - K r1 r2): the landscape
  law with q = K. Circle = tan addition (sphere), hyperbola = velocity addition (hyperbolic plane), K = 0 flat.
  Classical Cayley–Klein / Yaglom.
- The law's ruler is distance: ds = dr/(1 + K r²) on the straight-line (gnomonic / Beltrami–Klein) chart, whose
  2-D metric has curvature K (sympy, Brioschi). Walls are the absolute; ruler = log cross-ratio / √(p² - 4q)
  (Cayley 1859, Klein 1871). Circle-type lines close after π/√K (elliptic plane).
- Every landscape point has a curvature: K_proj = q - p²/4 (straight charts), K_conf = 4q - p² (angle-true charts),
  both unchanged by chart changes; doubling 2a/(1 - qa²) links the two (Poincaré → Klein). Ours as a placement.
- Mult-add lives in the upper half plane: height y = 1 + x, mult-add = dilation (isometry of (dx²+dy²)/y², K = -1),
  ruler log(1 + x) = distance, wall -1 = horizon. Dial k: curvature -k², radius 1/|k|; flat at k = 0.
- Velocity addition: K = -1/c², velocity space hyperbolic of radius c (Varićak 1910). Einstein addition = Klein-ball
  translation, Möbius addition = Poincaré; doubling carries one to the other (Ungar), checked in 2-D.
- Gauss–Bonnet (checked on one hyperbolic triangle) is only mentioned, as the reason for the twisted-boosts rotation;
  that page has no published link in the files, so it is named as "the twisted boosts page" without a URL.

Verdict: no independent lead. The geometry is sure and classical and gives every landscape point (and mult-add) a
picture, but in one dimension curvature has no measurable effect; the 2-D effect (holonomy) is the twisted boosts work.
