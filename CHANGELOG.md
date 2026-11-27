# Changelog

Each `## vX.Y.Z` section below becomes the notes of that GitHub release (and its Zenodo version).

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
