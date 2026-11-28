# Lean: DialPointProcesses

> **Names:** some words or symbols on this page (mult-add, wall) also mean something else in standard math, physics or computing, or inside this project. [NAMING-COLLISIONS.md](../../../NAMING-COLLISIONS.md) says which meaning is in play.

- File: `DialPointProcesses.lean`, namespace `DialPointProcesses`, imports `FusedAlgebra.Base` only.
- 40 theorems (38 on the topic, 2 small helpers: `prodL_zero_of_mem`, `sq_nonneg'`).
- Compiles clean with Lean 4.37.0-rc1, core only (no mathlib), no `sorry`.
  `#print axioms` on tilt_lens, wall_mode_fills, fermion_subpoisson, var_superpose, tilt_in_ruler,
  modeVar_wall: only propext, Classical.choice, Quot.sound (lines removed afterwards).
- Intended module path: `FusedAlgebra/Topics/DialPointProcesses.lean`.

Compile:

    export PATH=/opt/lean-4.37.0-rc1-linux/bin:$PATH
    LEAN_PATH=/tmp/lib/build lean lean/DialPointProcesses.lean

Parts: (1) hole lift 1 + k·fold = ∏(1 + kμ), named dial values, a wall eigenvalue fills the region;
(2) mode variance μ(1 + kμ), Var = Σμ + kΣμ² (Fano line), sub/super-Poisson, second factorial moment;
(3) superposition (k, K) × c = (k/c, cK) and its variance; (4) tilt ruler and the lens law at
p = 1 + 2kμ, d = kμ(1 + kμ), parabola p² − 4d = 1, named points (the circle operation ∘₁ (called mult-add before 1.0), union of failures,
velocity addition with speed limit 2), particle-hole symmetry.
