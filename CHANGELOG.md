# Changelog

Each `## vX.Y.Z` section below becomes the notes of that GitHub release (and its Zenodo version).

## v0.3.0

Third release, one day after v0.2.0 (Nov 27 to Nov 28, 2026). As before, almost everything new is marked **DRAFT**: checked by its own scripts and by CI, not reviewed by anyone outside the project, and not claimed as a discovery unless proved.

### Headline: what time could be (topics/time-study/, DRAFT)

- One shared time axis survives a passing gravitational wave, provided each place's winding rate may change in time: a four-leaf pattern of faster and slower winding then carries the wave's whole tidal pull, and a rate fixed in time fails.
- Read as a composition chain (each point a function that passes its output to the next), the law gives the time axis's tilt, the arrow of time and the speed limit (tilt/). From the rate pattern, the wave's direction and c alone, the twist comes out as exactly i × the rate pattern, with no field equations (twist/). The four-leaf shape follows with one input added, "balanced for every observer", which is what Einstein's vacuum equation says for a wave (shape/).
- Everything here matches general relativity; there is no new prediction yet, and the parts marked ours have not yet been checked against the published literature. 35 checks across the four folders, all run in CI.

### Also new since v0.2.0 (all DRAFT)

- Unimodular gravity as the one place in gravity where the why-i derivation of i goes through unchanged (topics/why-i/unimodular/; the physics is known, the reading is ours), and the why-i paper widened to a small bridge toward gravity (topics/why-i/bridge-paper/).
- The quantum–gravity bridges as one study (topics/bridges-study/), with its three candidate bridges run and four follow-ups, ending with the clock slowdown and the helix as the real and imaginary parts of one composition.
- The spinning river, rounds 23 to 34 (topics/gravity-from-theory/spinning-river/), now 152 checks; the spin-i-tie study and its sub-studies (topics/spin-i-tie/); Navier–Stokes as crossing flat flows, rounds 1 to 5, ending in a blow-up hunt that finds no blow-up in real flow (topics/ns-two-planes/); every force as curvature (topics/forces-as-curvature/); Riemann studies on the cusp tail, the cusp size and the k dial (topics/riemann/); new-ties rounds 14 to 22; torsion as the law's order gap (topics/torsion-dial/); and the physics bridges survey (topics/physics-bridges/).
- A naming review (topics/naming-review/) and a running known-vs-yours tally, [KNOWN-VS-YOURS.md](KNOWN-VS-YOURS.md): 65 known, 31 yours, 31 unsure.
- The Klein-bottle blow-up check is closed and kept as a one-time record, out of the automatic checks.

### The library and the build

- The shared Lean 4 library grew from 3,875 theorems in 181 modules to 3,957 in 183. Still core Lean only, no `sorry`, only Lean's three standard axioms.
- Per-push CI now runs the Lean builds and the quick, unit-test-style checks in three parallel jobs; the slow numeric runs moved to a separate workflow started by hand.

## v0.2.0

Second release, two days after v0.1.0 (Nov 25 to Nov 27, 2026). Most of what is new here is marked **DRAFT**: a first write-up, checked by its own scripts and by CI, but not reviewed by anyone outside the project. Nothing in this release is claimed as a discovery unless it is proved, and even the proved results below are unreviewed.

### Headline: shared-centre fractals (topics/new-ties/, DRAFT)

- When several maps share a centre, the size (dimension) of the fractal they build solves one closed equation. On the line, Bárány and Verma (arXiv:2507.05835, 2025) already proved this formula for almost every choice of parameters; that part is theirs. The plane version, and exact sizes for whole families of specific examples on the line and in the plane, are ours.
- Proved here: the size is exactly min(d, D*) for an explicit dense family with any pattern of class sizes; any accidental collision drops the size below D*, with a concrete break; and on the collision surface, under three stated conditions, the size of a break is exactly min(1, s_F) (resting on Hochman 2014).
- Rounds 10 to 13 attack the original break: no second collision for whole ranges of bad pairs covering 72.7% of the windows (computer-assisted, exact), and for every bad pair with a ≤ 10⁶. The original break's exact size is still open; the write-up explains why this route cannot finish it.
- 58 Lean theorems in the library (Topics.NewTies, Topics.NewTies11) and 178 numeric checks, all run in CI.

### The library and the build

- The shared Lean 4 library grew from 3,003 theorems in 119 modules to 3,875 theorems in 181 modules. It is still core Lean only, with no `sorry` and only Lean's three standard axioms.
- A Mathlib version of the library, built in CI, plus a Mathlib map (sixteen modules, 71 theorems) linking our results to Mathlib's.
- CI now also compiles every topic's own Lean file and runs every topic's numeric checks (every `checks*.py` plus 80 other scripts) against their saved output.

### Names

- The project now has two layers: Fused Hyperoperation Algebra (the algebra: a ∘ₖ b = a + b + k·ab, the plane law, the hyperoperation rungs) and Affine Composition Theory (the theory built on it). Mult-add is now called the circle operation, its name in ring theory. The repository and citation keep the algebra's name until the 1.0 name question is settled.

### New write-ups (89 new topic folders, almost all DRAFT)

- Mathematics: the union-find and order dials, decrement paths, the q-deformation dial, the Riemann studies, Lee–Yang on the ℓp ball (partly proved), the Tate curve drafted for Mathlib, hypersphere and division algebras, hat primes, and dial surveys across shapes, ratios and odds.
- Physics bridges, all DRAFT and none claimed as new physics: the spinning river (rebuilding general relativity's slow-orbit field from three rules plus the field's own energy; it matches the known 2PN results, and no prior-art check has been done), the new-force study, the Conley and Conley–Zehnder indices on a three-body return map, Newton's G, three-body orders, the Klein-bottle flow, quantum-to-classical, and others.
- [INDEX.md](INDEX.md) lists every page with its status; [IN-PROGRESS.md](IN-PROGRESS.md) logs what changed and when.

## v0.1.0

First public release of Fused Hyperoperation Algebra.

- The mult-add framework (a + b + k·ab) and the one-picture idea behind it.
- A shared Lean 4 library of 3,003 theorems in 119 modules, core Lean only, no `sorry`, standard axioms only.
- Lens results for general relativity and electromagnetism: where mult-add sits in gravity, and how light sets the dial q = −1/c². No new predictions are claimed.
- Licenses: Apache 2.0 for code (LICENSE), CC BY 4.0 for write-ups, figures and generated data (LICENSE-docs).
