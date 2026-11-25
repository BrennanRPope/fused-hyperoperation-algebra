# Linear algebra from mult-add (Nov 25, 2026)
Page: https://claude.ai/artifact/6MmBPSMK179Q7ZRFTjz3L3 (slug 2b63cb86-ecf4-43fb-bbfe-a127715d93e0, version 1795589645-7228)
- checks.py: 37/37 pass (checks_output.txt)
- lean/LinearAlgebra.lean: 51 theorems, core Lean style, NOT COMPILED (no toolchain in the writing session).
  Fold in as FusedAlgebra/Topics/LinearAlgebra.lean. Heaviest proof: M3.det_mul (grind, maxHeartbeats 2000000).
- Left out of Lean: n×n; non-triangular eigenvalues; exterior/symmetric powers and alpha-permanent dial;
  3x3 Cauchy factoring; log/exp. 3x3 dial homomorphism proved over fields only (rings: sympy).
Findings: see page ledger. Lead: KdV Hirota cube = det(I+M), M Cauchy-type (Sato; Poppe 1984).

## Dial sweep (Nov 25, 07:08Z) — page version 1795590492-ba79
- checks_dial.py: 31/31 pass (checks_dial_output.txt)
- lean/DialSweep.lean: 28 theorems, imports Topics.LinearAlgebra, NOT COMPILED. Fold after LinearAlgebra.
- Finding: every nonzero real k is one law (units); only k=0 (trace) and k→∞ (det) differ. Walls = eigenvalues
  (det(I+kA) = reversed char poly; Fredholm λ = −k). k=i: Hermitian→unitary Cayley.
- Biggest: dial matrix U=[[0,−q],[1,k]] gives landscape lens law; det(I+aU)=1+ka+qa²; walls −1/eig;
  countries = discriminant; or line = eigenvalue −1; κ units = rotation/shear/boost (light cone = wall). Known math
  (Möbius, Cayley–Klein/Yaglom), reframing.
- Left out of Lean: limits k→0, k→∞; complex k (no ℂ in core; units as real 2×2); n×n Sylvester.
