# DRAFT · Phase winding at E/ħ and gravity's clock slowdown

Status: DRAFT, not checked by anyone outside this project. Brennan's picture: each clock's line is the centre of a helix; its quantum phase winds around it at E/ħ.

## Result

Gravity's clock slowdown comes out as slower winding. The same composition also gives the helix.

- Take the law as it is (k = 1) and give it an imaginary argument, x = i·v. Then 1 + x∘x = (1 + iv)².
  - Its real part, 1 − v², is the clock slowdown squared.
  - Its imaginary part over its real part, 2v/(1 − v²), is the helix's turns per loop.
- Real time uses x = i·v. Imaginary time turns this into velocity addition and tanh u.

This also closes the gap left in ../helix-from-law/: the squared 1/γ² in the horizon's surface gravity is now derived, not taken from Kerr.

1. **No spin (W1).** Seen from far, a clock winds at (E/ħ)·N, with N² = Re(1 + x∘x) and x = i·(river speed). Winding stops where the river reaches light speed, at r = 2M.
2. **Kerr (W2).** The shifted source 1/(r − ia) carries both pieces. Its real part gives the shells' clock factor Δ/(r²+a²), which is zero at the horizon. Its imaginary part is the shells' turn rate a/(r²+a²).
3. **Why squared (W3, W4).** Surface gravity is the slope (the pull) of the potential. The slope of 1/(r − ia) is 1/(r − ia)², and that square turns the rim's 1 − v into Re(1 + iv)² = 1 − v². So κ = (1 − v²)/4M follows from the shifted source.
4. **The helix (W5).** Ω_H/κ = Im/Re of the same (1 + iv)². The slowdown and the helix are the real and imaginary parts of one composition. Continuing to a = iα gives tanh u.
5. **Bigger mass winds faster (W6, known).** A mass's own phase per imaginary-time loop is M·β = 2S without spin, and β(M − 2Ω_H J) = 2S with spin (Smarr). The winding rate grows with M and the loop gets longer, so the phase per loop grows as M².
6. **Where G enters (W7).** The source's own winding rate gives the river speed only once the Planck length (so G) is put in. So the slowdown follows from winding plus the river speed, but the river speed still needs G. That's the same open item as ../../newton-g/.

## Correction to ../helix-from-law/

The "added dial c" there is not new. It is the repo's d dial, from the two-dial law (x + y + kxy)/(1 − d·xy) of Coleman and McGuinness 1991 (see ../../boolean-lattice/, ../../action/, ../../affine-sphere/). At k = 0 it is the landscape q of ../../twisted-boosts/ (d = −1 relativity, 0 Galilean, +1 rotations, on half-velocities). The write-up now calls it d. The rim speed t = tanh(u/2) found there is that project's half-velocity, and helix = t⊕t is its "half a kick, doubled". Seen this way, no extension of the law was needed: an imaginary argument does the same job.

## Checks (checks.py, 7/7)

| Check | What it tests | Result |
|---|---|---|
| W1 | Schwarzschild winding factor N² = Re(1 + x∘x), x = i·v_river; stops at 2M | pass |
| W2 | Kerr shifted source: real part = shell clock factor, imaginary part = shell turn rate | pass |
| W3 | κ = M·Re[1/(r₊−ia)²] = (1−v²)/4M | pass |
| W4 | 1 − v² = Re(1 + x∘x) with x = i·(rim speed) | pass |
| W5 | helix = Im/Re of (1 + x∘x); tanh u in imaginary time | pass |
| W6 | phase per loop = 2S (Smarr) | pass |
| W7 | source winding gives river speed only with l_P² (G) | pass |

Where exact simplification stalls, checks fall back to 40-digit evaluation at 9 points. The fallback was confirmed to reject false identities.

Figure: figures/c2g_winding_slowdown.png. Panel a shows one clock's helix at four distances, panel b the composition (1 + iv)² in the complex plane, and panel c phase per loop against 2S.
