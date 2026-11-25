# Complex analysis through mult-add (Nov 25, 2026)

Finding: in the complex plane z ✳ w is a turn and stretch about the wall −1 by 1 + w; the wall is the branch point of
log(1 + z), the fixed point of every landscape Möbius flow, and the place a product's pieces land at its zeros.

## Files
- index.html: the page (self-contained; phase picture embedded as a data: JPEG; all other figures drawn as SVG by script).
- checks.py -> checks_output.txt: 50 of 50 checks pass (numpy, sympy, mpmath). results.json: the ripple coefficient a1,
  Abel-ruler coefficients, the five ripple estimates.
- figures.py -> figures/phase.png, phase.jpg, phase_b64.txt (domain colouring of 1 + z, hue = angle seen from the wall).
- lean/ComplexAnalysis.lean: 39 theorems, compiles clean, no sorry, standard axioms. See lean/README.md.

## Findings
1. z ✳ w + 1 = (1 + w)(z + 1): mult-add by w is the similarity about −1 (turn arg(1+w), stretch |1+w|); triangles
   (−1, 0, w) and (−1, z, z ✳ w) are similar. Affine maps fixing −1 are exactly z ↦ z ✳ c. Known (Wessel, Argand).
2. log(1 + z) turns ✳ into + mod 2πi; the branch point is the wall, not 0. e^w − 1 (the rung step) sends the plain grid to the
   polar web about the wall. Classical.
3. Mult-add roots of 0 are ζ − 1, a regular n-gon on |1 + z| = 1 through 0; all n of them mult-add to 0 (n odd) or −2 (n even).
4. Landscape laws as Möbius maps: translations compose by the law (projective matrices), det = wall polynomial, fixed points =
   walls, tr² − 4det = a²(p² − 4q). Over C two walls -> conjugate to mult-add (Cayley for tan, Doppler for tanh); one double
   wall -> addition. Textbook one-parameter Möbius groups.
5. e^z − 1 at 0 is parabolic: Abel ruler α = −2/z + (1/3)log z − z/36 + z²/540 + 1/7776 z³ − 71/435456 z⁴ ...; the attracting and
   repelling rulers differ by a1·e^{2πiw}, a1 ≈ 10.748807 + 63.246298i (|a1| ≈ 64.153), decay rate 2π confirmed. Nonzero
   ripple = no analytic fractional rung on a neighbourhood of 0 (Écalle–Voronin; Baker 1958). The number itself was computed
   here in this normalisation; not found published.
6. Products are mult-add folds; zeros = wall hits (sine at ±n, 1/Γ at −n); |fold| ≤ fold of sizes (Weierstrass bound, Lean);
   (−1)^n/√n shows Σa_n convergent is not enough. Argument principle with f'/(1+f) counts wall hits (sine fold: 6 in |z| < 3.5).
7. Short linking sections: Lee–Yang page (not redone) and the elliptic page's Tate torus (gluing respects ✳, Lean).

## Verdict
No new mathematical lead: this is the group of nonzero complex numbers seen through the lift, all classical. It gives the
project a sure geometric footing (the wall is a genuine puncture). One thread with content: the rung step's ripple a1, which
measures how far fractional rungs are from one analytic flow; worth following only if a canonical half rung is wanted (none exists).

Links used: Lee–Yang UQQzmkzekxe2hL4m8F8NQL, linear algebra 6MmBPSMK179Q7ZRFTjz3L3, fractional rungs section on NS page
35PaQq5t56ZxNhHjKdq8pJ, calculus VBRCmFrQVwQiXjMxZkY5SY, elliptic JhsmEQaHNRhZ9t5Ph5RZ6y.
