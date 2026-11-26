# Two layers: Fused Hyperoperation Algebra under Affine Composition Theory

Decided by Brennan, Nov 26, 2026 ("Set it up that way", answering whether Fused Hyperoperation Algebra can be its own thing that Affine Composition Theory relies on).

## The split

| Layer | What it holds | Depends on |
|---|---|---|
| **Fused Hyperoperation Algebra** (the algebra) | The fused operation a ∘ₖ b = a + b + k·ab and the plane law (a + b + k·ab)/(1 − d·ab); their laws (associativity, identity, inverses, walls, chains); the hyperoperation rungs (fused add, multiply, power, tetration) with the hypernatural logarithm and exponential; the Lean library (`FusedAlgebra`). | Nothing in the theory. It stands alone. |
| **Affine Composition Theory** (the theory) | The shapes (ℓp exponent p from cross-polytope to cube, cubic form C to the simplex, the cube as one setting); the dials as settings (k, d, p, C, t, r, β) and the readings (Δ, curvature −Δ/4, D, holonomy ω); conics and affine spheres; order, which isn't required but a dial turns on (q, or ħ over the star product; the checkerboard's −i); the physics bridges. | The algebra. |

The dependency runs one way only. Nothing in the algebra mentions shapes, dials beyond k and d, or physics. That is the same relationship linear algebra has with geometry.

## Names

- The repo keeps the name **fused-hyperoperation-algebra**, which is the algebra's name.
- The README presents both layers: the algebra first, then the theory built on it.
- CITATION.cff, .zenodo.json and the version number stay unchanged until the 1.0 name question is settled.
- Earlier names of the theory: Hypercube Composition Theory (Nov 25), and Fused Mult-Add Algebra / Fused Hyperoperation Algebra for the whole project before that.

## Which topic folders sit where (first pass, by folder name; to be confirmed as the README is rebuilt)

- **Algebra:** axioms, abstract-algebra, lean-library, operations-catalog, named-functions, fractional-rungs, fractional-steps, rung-shape, height-flow, e, fused-constants, boolean-lattice, semiring-audit, mathlib-map, mathlib-plan, mathlib-port-check, proof-engine, proof-gallery, simple-laws.
- **Theory:** everything else: shapes (shape-dial-survey, simplex, hypersphere, tesseract, exceptional-shapes, beyond-the-cube), dials (affine-names, dial-survey, order-dial, repetition-dial, zoom-dial), geometry (conics, affine-sphere, curved-space, differential-geometry), and the physics, chemistry and biology passes.
