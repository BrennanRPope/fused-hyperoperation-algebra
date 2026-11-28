# DiffGeometry.lean

> **Names:** some words or symbols on this page (mult-add, wall) also mean something else in standard math, physics or computing, or inside this project. [NAMING-COLLISIONS.md](../../../NAMING-COLLISIONS.md) says which meaning is in play.

- 32 theorems. Compiles clean with core Lean 4.37.0-rc1 against the shared library build
  (`/tmp/lib/build`): no errors, no warnings, no `sorry`. `#print axioms` on cr_base,
  dial_flat_iff, rmul_form, central_chart and cr_num showed only propext, Classical.choice,
  Quot.sound (the `#print` lines were removed afterwards).
- Imports only `FusedAlgebra.Base`. Namespace `DiffGeometry`. Defines `rmul`, `form`, `dN`, `dD`,
  `curv4`, `CR`.
- Intended module path: `FusedAlgebra/Topics/DiffGeometry.lean`
  (add `FusedAlgebra.Topics.DiffGeometry` to build.sh and an import line to FusedAlgebra.lean).

Compile:

    export PATH=/opt/lean-4.37.0-rc1-linux/bin:$PATH
    LEAN_PATH=/tmp/lib/build lean DiffGeometry.lean

Sections: (1) turns of t² + Kx² = 1 keep the form, and the central chart composes by F(0, K);
flat case adds; (2) doubling 2a/(1 − da²) (Poincaré → Klein, stereographic → gnomonic) respects the
law, in cleared numerator/denominator form; (3) 4d − p² is chart-free and scales by λ²; values for
the circle operation ∘₁ (called mult-add before 1.0), the dial, tangent, velocity, parity; real walls give negative curvature; (4) cross-ratio
with the walls is the lens; (5) circle operation as dilations of the upper half plane, isometries of
(dx² + dy²)/y², wall on the horizon; (6) Möbius and Einstein addition agree on a line.
Left to checks.py (needs real numbers): Gaussian curvatures, arc lengths, Cayley–Klein logs,
Ungar's 2-D doubling, Gauss–Bonnet.
