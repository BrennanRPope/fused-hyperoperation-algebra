# Special relativity's velocity addition through mult-add (Nov 25, 2026)
Page: https://claude.ai/artifact/S2RVDVVQWR4UAQVtmDQQMq (source index.html)
- checks.py: 37/37 pass (checks_output.txt)
- lean/Relativity.lean: 34 theorems, compiled with core Lean 4.37.0-rc1 against the shared library (imports Base, Topics.Thermodynamics). No sorry, standard axioms only. Fold in as FusedAlgebra/Topics/Relativity.lean after Thermodynamics.
Findings: velocity addition is mult-add in d = 2u/(1-u) (Doppler factor squared minus 1); rapidity is half the log of the lift; landscape point k = 0, q = -1/c^2 (no corner), walls ±c; rung s has shape s, and in natural units rungs 0 to 1 run along q = k^2/4 - 1 from velocity addition to mult-add (third layer only).
Next lead: boosts in different directions (Thomas–Wigner rotation), Lorentz group as Möbius maps.
