# DRAFT · Deriving the Kerr helix from the law (on the d dial)

Status: DRAFT, not checked by anyone outside this project.

**Correction (Nov 28):** the dial below is not new. It is the repo's d dial, from the two-dial law of Coleman and McGuinness 1991 (see ../../boolean-lattice/, ../../action/, ../../affine-sphere/); at k = 0 it is the landscape q of ../../twisted-boosts/. It is called d below. The missing 1/γ² is now derived in ../winding-slowdown/. Follows ../c2-one-body/, where the helix came from the horizon, not from the law.

## Result

The helix comes out of the law on the law's existing d dial. One piece is still taken from Kerr, not derived.

1. **The two-dial law.** a∘b = (a + b + k·ab)/(1 − d·ab). It's still associative and commutative, 0 is still the identity, and d = 0 gives back the law exactly (K1, K2). It's the law conjugated by f(x) = (1+px)/(1+qx), with k = p+q and d = pq. This is the Coleman–McGuinness two-dial law, already in the repo. Reading the helix on it is ours.
2. **What the d dial does.** With k = 0, d = −1 is relativistic velocity addition (boosts) and d = +1 is tangent addition (rotations) (K3). Multiplying by i flips the sign of d, which is the project's "two legs meet at i" shape showing up inside the law.
3. **A Kerr horizon is a spun Schwarzschild horizon.** Its equatorial radius is 2M, the same as Schwarzschild. Its rim moves at v = a/r₊ (K4). Its surface gravity is Schwarzschild's 1/4M times the rim's 1/γ², so κ = (1 − v²)/4M (K5).
4. **The helix is the rim composed with itself.**
   - Real time: Ω_H/κ = v∘v at k = 0, d = +1, which is tan χ (K6).
   - Imaginary time (a = iα): the turn per loop is tanh u = t∘t at k = 0, d = −1, with t = α/r₊ = tanh(u/2) (K7).
5. **The law's multiplier sits at the horizon pole.**
   - Real time: (r₊ − M + ia)/M = e^{iχ}, on the unit circle.
   - Imaginary time: (r₊ + α − M)/M = e^u, on the real line. This is the distance from the moved mass to the horizon pole.
   - The helix is tan χ or tanh u of that multiplier (K8).
6. **Side finding (Kerr–Newman).** κ·2R = 1 − v² − Q²/r₊². The slowdowns from spin and charge add with no cross term (K8).

## What's derived and what isn't

- From the law (two-dial form): the helix has the form of a self-composition, the dial settings are d = ±1, and i swaps them.
- Taken from Kerr: the rim's 1/γ² factor in κ. Part 8's "clock rates multiply" gives one factor of the rim's clock rate (1/γ). Why κ takes it squared is the open step.

## Does it open a lead?

Yes, a modest one. The d dial is where signature lives (rotation versus boost), and i flips it. That's a concrete place for the shared shape inside the law itself. The 1/γ² slowdown is now derived from the shifted source in ../winding-slowdown/.

## Checks (checks.py, 8/8)

| Check | What it tests | Result |
|---|---|---|
| K1 | two-dial law is multiplication conjugated by a Möbius map | pass |
| K2 | associative for all k, d; d = 0 is the original law | pass |
| K3 | d = −1 velocity addition, d = +1 tangent addition | pass |
| K4 | Kerr horizon radius 2M; rim speed a/r₊ | pass |
| K5 | κ = (1 − v²)/4M | pass |
| K6 | real-time helix = v∘v at d = +1 | pass |
| K7 | imaginary-time helix tanh u = t∘t at d = −1 | pass |
| K8 | multiplier at the pole; Kerr–Newman slowdowns add | pass |

Figure: figures/c2f_helix_from_law.png. Panel a is the dial plane, panel b the helix as rim∘rim with actual Kerr horizons as dots, panel c the multiplier on the circle and on the line.
