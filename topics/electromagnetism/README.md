# Electromagnetism through the circle operation (Nov 25, 2026)
Page: https://claude.ai/artifact/A71iTBsgUHZSkcbNUrEacY (slug 49b926a7-cb25-47e8-9c1f-da8a404822e9, version 1795626602-cf6b;
source index.html, built by build_page.py from figures/*.svg and results.json)

Task (Brennan, Nov 25, voice): "let's look at electromagnetism" — show how Maxwell fits under the circle-operation lens, in historical
order (EM, then special relativity, then gravity: "general relativity... stood on the shoulders of electromagnetism").

## Answer
- Same structure (exact, known physics): (1) the wave equation survives a dial boost only at d = -1/c^2, and c^2 eps0 mu0 = 1:
  EM sets the dial; (2) boosts scale the light-cone parts of E, B by D and 1/D, ratio by D^2 = 1 + d(v) (the lift); D's multiply,
  d's combine by the circle operation ∘₁ (called mult-add before 1.0); light is a null field on the wall, energy density x 1/lift; (3) loop phase = charge x flux, same "loop turn =
  curvature x area" as gravity; lattice plaquettes are the circle operation of four links, neighbours combine by the circle operation (discrete Stokes).
- Renamings: phases as the circle operation of e^{i theta} - 1 (any multiplicative group); F as the wedge part of the "corner" grad A
  (content = Clifford's dot + wedge split); grad F = J is Hestenes 1966.
- Small frame of ours: gauge phase (d = +1, circle) and boost (d = -1, hyperbola) both obey a * conj(a) = 0.
- EM linear vs gravity nonlinear = commuting phases (the plain circle operation) vs non-commuting boosts (sandwich law, order gap = commutator).
- Stops: source law grad F = J from outside (as Einstein's equations); vacuum superposition is plain addition (k = 0);
  Born-Infeld shares the wall shape but not the composition law (ATLAS: scale > ~100 GeV, arXiv:1703.08450).
- No new prediction, no lead toward new physics. Kept: the bridge to gravity as a gauge field (teleparallel river).

## Files (for the repo thread to push; nothing pushed from here)
- checks.py -> checks_output.txt, results.json: 49/49 pass (sympy, numpy, scipy). Uses ../gr-bridges/ga.py by path.
- lean/Electromagnetism.lean: 41 theorems, core Lean 4.37.0-rc1, no sorry, standard axioms only (checked with #print axioms).
  Compiled against the shared library's build/ (imports Base, Topics.Thermodynamics, Topics.Relativity, Topics.TwistedBoosts).
  NOT yet merged into lean-library: to fold in, copy to FusedAlgebra/Topics/Electromagnetism.lean, add
  FusedAlgebra.Topics.Electromagnetism to build.sh after GRBridges (needs TwistedBoosts earlier, which it is) and an import line in
  FusedAlgebra.lean, then ./build.sh. It defines its own local sq_nn (Electromagnetism.sq_nn), so no name clash.
- figures.py -> figures/*.svg (history, dial, field, circle, holonomy, commute, corner, one, walls).
- build_page.py -> index.html.
