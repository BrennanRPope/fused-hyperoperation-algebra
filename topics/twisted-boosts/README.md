# Twisted boosts: non-collinear boosts in the circle-operation picture (Nov 25, 2026)
Page: https://claude.ai/artifact/3ozW93quvmikMu5TGaKuw4 (source index.html)
- checks.py: 36/36 pass (checks_output.txt). Exact sympy algebra; numpy against 4x4 Lorentz matrices (2000 random pairs), SO(3) for the sphere end, triangle areas, Thomas precession by composing 20,000 small boosts.
- lean/TwistedBoosts.lean: 32 theorems, core Lean, imports Base, Thermodynamics, Relativity. Compiled in the shared library as Topics/TwistedBoosts (clean build 1,616 theorems, 54 modules); one proof step fixed there. Statements checked true on 300 random rational cases.
Findings:
- Dial boost B_d(u) = [[1,u],[-d ū,1]] on half-velocities (t = tanh(φ/2)); d is the landscape d (d=-1 relativity, 0 Galilean, +1 rotations).
- B(u)B(v) diagonal = X ± iY, X = 1 - q(u·v), Y = q(u×v): the corner ū·v splits into dot + i cross. Twist e^{iω} = (X+iY)/(X-iY), tan(ω/2) = q(u×v)/(1-q(u·v)). No twist iff d·cross = 0.
- u⊕v = (u+v)/(1-d ū v) = twist·(v⊕u). Light wall in 2D: |1+ūv|²-|u+v|² = (1-|u|²)(1-|v|²).
- ω = d × signed area of the velocity triangle (rest, u, u⊕v): Gauss–Bonnet with curvature K = d. Thomas precession per lap = 2π(γ-1) = area of the velocity circle.
- Half-velocity = half a kick: d(t⊕t) = d(t) ∘₁ d(t).
Known: Thomas 1926, Wigner 1939, Borel 1913, Ungar gyrogroups, κ-stereographic model. Ours: landscape d = curvature, twist = phase of split corner.
Next lead offered: Gravity Probe B geodetic precession (~6.6 arcsec/yr) with twist = curvature × area bookkeeping on spacetime.

## Gravity Probe B follow-up (gpb/)
- gpb/gpb_checks.py: 21/21 pass. Exact Schwarzschild parallel transport gives -2π√(1-3M/r) per orbit vs radial; space slice K = -M/r³, curvature outside orbit = 2π(1-√(1-2M/r)); v² = M/(r-2M), γ² = (r-2M)/(r-3M); √(1-3M/r) = √(1-2M/r)/γ, so −p = (−p_space) ∘₁ (−p_velocity).
- GP-B at r = 7027.4 km: velocity twist 2201.3 + space twist 4402.6 = 6603.9 mas/yr vs measured 6601.8 ± 18.3 (GR with Earth's shape 6606.1). Either alone off by >100 sigma. Frame dragging ≈ 40.8 vs 37.2 ± 7.2 (GR 39.2), as gravitomagnetic circulation per area.
- Caveat: exact velocity factor is 1/γ, Thomas's flat-space twist is γ−1; equal to first order only.
- gpb/lean/GravityProbeB.lean: 16 theorems, sent to library for first compile.
- Page updated (version 1795594639-04ee), section #gravity-probe-b.

## Binary pulsar follow-up

Checks: `pulsar/pulsar_checks.py` (13/13 pass). Lean: `pulsar/lean/BinaryPulsar.lean` (10 theorems, merged unchanged as Topics/BinaryPulsar; library builds clean at 2,263 theorems).

- Three twists, each weighted by the companion's mass m₂: velocity ½·m₂/M, space m₂/M, companion's drag 2·m₁/M. Their sum times m₂ is m₂(4m₁ + 3m₂)/(2M), the standard spin-orbit coefficient. Light-gyroscope limit gives Gravity Probe B's 3/2.
- Double pulsar, pulsar B: 0.75 + 1.51 + 2.81 = 5.07°/yr (one turn per 71 years), against 4.77 +0.66 −0.65 measured (Breton et al. 2008, arXiv:0807.2644). Within 1σ.
- Without the drag twist: 2.26°/yr, 3.9σ off, so the data needs it. Without the velocity twist: 4.32°/yr, still within 1σ, so the pulsar can't isolate it.
- Correction to the Gravity Probe B notes: the exact ∘₁ fusion can't be tested here either. M/r ≈ 4.35e-6 (about 7,000× Gravity Probe B's), and the corner is still about 1e-6 of the effect, far below the 14% error.
- Hulse–Taylor: same bookkeeping gives 1.21°/yr.
- Verdict: the twist law reproduces a second, independent measured precession, including the drag twist that a heavy partner adds. No further lead with data was judged fairly promising (strong-field gyroscope data doesn't exist; periastron advance as a twist is speculative), so the work stops here.
