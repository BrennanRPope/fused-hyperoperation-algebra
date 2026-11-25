# Dynamical systems through the circle operation (Nov 25, 2026)

Files
- checks.py -> checks_output.txt: 53/53 pass (about 2 minutes; the horn-map section runs
  mpmath at 34 digits). Writes series.json (rung-flow vector field and Abel-ruler coefficients,
  horn-map data, c₁) used by the page.
- lean/Dynamics.lean: 31 theorems, compiles clean (see lean/README.md).
- index.html: the page (inline SVG drawn by script; the Abel ruler is computed live in JS).

Findings
- x ∘₁ x is z² through the lift (classical): wall −1 superattracting ((−1+h)∘₁(−1+h) = −1 + h²),
  identity 0 repelling with multiplier 2, Julia set |1+x| = 1 with angle doubling (Lyapunov log 2).
  Every dial square x ∘ₖ x is the same map in other units; k = −2 is the logistic map r = 2.
- Logistic r = 4 is x ∘₁ x folded by X(w) = w·ι(w)/4, ι the circle-operation inverse (2-to-1; circle folds
  to [0,1] as sin²(θ/2)). Classical (Ulam–von Neumann) via Joukowski; the ι reading is ours.
  Chebyshev Tₙ = zⁿ through Joukowski; tent map = folded doubling.
- Koenigs ruler of x ∘₁ x at 0 is log(1+x) (the project's log ruler). E = e^x − 1 is parabolic:
  no Koenigs ruler; Fatou/Abel ruler α = −2/x + ⅓log|x| − x/36 + x²/540 + …, rung-flow vector field
  v = x²/2 − x³/12 + x⁴/48 − x⁵/180 + 11x⁶/8640 − …; E^(1/2)(1) = 1.2710274 by the ruler.
- Horn map (Écalle–Voronin): α₊ − α₋ = c₁e^(2πiα) + … with c₁ ≈ −10.748807 − 63.246298i
  (|c₁| ≈ 64.153, ruler normalized with zero constant term; the size depends on that
  normalization, its being nonzero does not). So e^x − 1 has no analytic flow through 0 (Baker 1958).
  Computed here; not found in print.
- Landscape laws = time laws of Riccati flows y' = 1 + py + qy² (already on the landscape page).
- Feigenbaum: no genuine connection found; left out.

Verdict: mostly classical dynamics through the lift (a dictionary, not a discovery). One modest
lead: the horn map of the rung flow (compute c₂, c₃, …; cross-check with Dudko–Sauzin's Borel route).
