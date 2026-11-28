# Running tally: what the bridges study derives and what it assumes

Kept up to date with each round. Tags follow ../dependency-map/assumptions.md: **known** (published), **ours** (a step of ours), **input** (put in by hand). That ledger covers the main path, and ../KNOWN-VS-YOURS.md covers known-versus-ours across the repo; this one covers derived-versus-assumed in bridges-study/ only.

## Derived (from the law or the river, with checks)

| Result | Folder | Rests on |
|---|---|---|
| Two legs (entropy = area/4Għ, phase = action/ħ) meet at i | README | known pieces; the reading is ours |
| k = 0 is Dong's tensionless point (Rényi–Tsallis lift) | README | known (Dong 2016); the law's k = 1−n is ours |
| Area exponent δ is the existing roundness dial p = 1/δ | candidates/ | ours |
| Kerr T and S = A/4 from the source shifted by i·a | candidates/ | shift known (Newman–Janis); river reading ours |
| Moved mass with all moments from the exact shift | c2-one-body/ | known (Hansen 1974); river reading ours |
| River horizon shell rate = Kerr horizon rate, every angle | c2-one-body/ | ours |
| Helix = horizon rim composed with itself (d dial = ±1) | helix-from-law/ | two-dial law known (Coleman–McGuinness 1991, in the repo); helix reading ours |
| Surface gravity's (1 − v²) from the slope of 1/(r − ia) | winding-slowdown/ | ours |
| Clock slowdown and helix as Re and Im/Re of (1 + iv)² | winding-slowdown/ | ours |

## Assumed (not yet derived)

| Assumption | Where used | Tag |
|---|---|---|
| G's value; the river speed v² = 2GM/rc² needs it | winding-slowdown/, all horizons | input (open, see ../newton-g/) |
| ħ's value and the rule that phase winds at E/ħ | winding-slowdown/ | input |
| The shift by i·a itself (why spin is an imaginary displacement) | candidates/, c2-*, helix-from-law/ | known (Newman–Janis), not derived from the law |
| Horizon smoothness fixes β = 2π/κ | c2-one-body/, helix-from-law/ | known (Gibbons–Hawking) |
| Phase per loop = 2S | winding-slowdown/ | known (Smarr) |

## Failed or ruled out

| Claim | Folder |
|---|---|
| Spin = two half-masses at imaginary distance (fails at the dipole) | c2-two-body/ |
| Part 8's 1PN rules alone give the helix (no horizon, so no) | c2-one-body/ |
| Rényi-of-black-hole entropies with linear mixing (would need a brane) | candidates/ |
