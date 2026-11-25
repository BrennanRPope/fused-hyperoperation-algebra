# Twisted boosts: non-collinear boosts in the mult-add picture (Nov 25, 2026)
Page: https://claude.ai/artifact/3ozW93quvmikMu5TGaKuw4 (source index.html)
- checks.py: 36/36 pass (checks_output.txt). Exact sympy algebra; numpy against 4x4 Lorentz matrices (2000 random pairs), SO(3) for the sphere end, triangle areas, Thomas precession by composing 20,000 small boosts.
- lean/TwistedBoosts.lean: 32 theorems, core Lean, imports Base, Thermodynamics, Relativity. Compiled in the shared library as Topics/TwistedBoosts (clean build 1,616 theorems, 54 modules); one proof step fixed there. Statements checked true on 300 random rational cases.
Findings:
- Dial boost B_q(u) = [[1,u],[-q ū,1]] on half-velocities (t = tanh(φ/2)); q is the landscape q (q=-1 relativity, 0 Galilean, +1 rotations).
- B(u)B(v) diagonal = X ± iY, X = 1 - q(u·v), Y = q(u×v): the corner ū·v splits into dot + i cross. Twist e^{iω} = (X+iY)/(X-iY), tan(ω/2) = q(u×v)/(1-q(u·v)). No twist iff q·cross = 0.
- u⊕v = (u+v)/(1-q ū v) = twist·(v⊕u). Light wall in 2D: |1+ūv|²-|u+v|² = (1-|u|²)(1-|v|²).
- ω = q × signed area of the velocity triangle (rest, u, u⊕v): Gauss–Bonnet with curvature K = q. Thomas precession per lap = 2π(γ-1) = area of the velocity circle.
- Half-velocity = half a kick: d(t⊕t) = d(t) ✳ d(t).
Known: Thomas 1926, Wigner 1939, Borel 1913, Ungar gyrogroups, κ-stereographic model. Ours: landscape q = curvature, twist = phase of split corner.
Next lead offered: Gravity Probe B geodetic precession (~6.6 arcsec/yr) with twist = curvature × area bookkeeping on spacetime.

## Gravity Probe B follow-up (gpb/)
- gpb/gpb_checks.py: 21/21 pass. Exact Schwarzschild parallel transport gives -2π√(1-3M/r) per orbit vs radial; space slice K = -M/r³, curvature outside orbit = 2π(1-√(1-2M/r)); v² = M/(r-2M), γ² = (r-2M)/(r-3M); √(1-3M/r) = √(1-2M/r)/γ, so −p = (−p_space) ✳ (−p_velocity).
- GP-B at r = 7027.4 km: velocity twist 2201.3 + space twist 4402.6 = 6603.9 mas/yr vs measured 6601.8 ± 18.3 (GR with Earth's shape 6606.1). Either alone off by >100 sigma. Frame dragging ≈ 40.8 vs 37.2 ± 7.2 (GR 39.2), as gravitomagnetic circulation per area.
- Caveat: exact velocity factor is 1/γ, Thomas's flat-space twist is γ−1; equal to first order only.
- gpb/lean/GravityProbeB.lean: 16 theorems, sent to library for first compile.
- Page updated (version 1795594639-04ee), section #gravity-probe-b.
