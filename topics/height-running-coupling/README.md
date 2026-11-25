# DRAFT: the height flow read as a running coupling

> **DRAFT, not finished.** Pushed Nov 25 so it isn't lost. The work stopped partway when Brennan merged it into the height-flow stream, now [topics/height-flow/](../height-flow/README.md). No figures yet. [HANDOFF.md](HANDOFF.md) holds everything found so far.

**Early verdict:** the height flow of E(x) = eˣ − 1 behaves exactly like a marginal running coupling, with one-loop and two-loop running, scheme-free b0 and ρ = −b1/b0² = 1/3, and renormalon-like divergence whose singularity exponent is set by ρ. But any map with a double fixed point does the same, so nothing here is special to E. Its singularities sit on the imaginary Borel axis, not the real one as in quantum field theory. No physics test was found. It doesn't open a lead toward physics; it is worth keeping as a correct dictionary entry between iteration theory and the renormalization group.

**Known vs ours:** the pieces are known on both sides (Julia, Écalle, Dudko–Sauzin, 't Hooft, Beneke, Curtright–Zachos 2011). The side-by-side dictionary for E is ours. Curtright–Zachos should be read before claiming anything new.

Run `python3 running_coupling.py` (14 of 14 checks pass; output in `running_coupling_output.txt`, numbers in `running_coupling.json`).
