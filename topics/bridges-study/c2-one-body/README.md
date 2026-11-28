# DRAFT · One spinning body: do the river's rules give the moved mass and the helix?

Status: DRAFT, not checked by anyone outside this project. Follows ../c2-two-body/ (which found the "two half-masses" reading fails).

## Result: partly yes

1. **The moved mass:** the river's *exact* shift of the source by i·a gives it, with every moment M(ia)^l. Part 8's moving-body rules (1PN) give only the first step: the spin's own drag term, continued to a = iα, is exactly the dipole of a mass moved by −α (H1). The second step, 2Ma²cos²θ/r³, is needed and is missing from those rules (H2).
2. **The helix:** the river's rigid shells turn at a/(r²+a²). At the horizon this equals Kerr's frame-drag rate at every angle (H3). Run in imaginary time, that shell rate times the loop length β gives the Euclidean Kerr turn per loop 2πα/√(M²+α²) = 2π tanh u (H4).
3. **What the 1PN rules can't do:** they give the far swirl 2J/r³, which continues to a real 2Mα/r³, but they have no horizon, so they can't set β and can't give the helix (H5).

So the picture holds together for one body: the moved mass comes from the exact shift, and the helix comes from the horizon shell. The moving-body rules carry the first step only.

## Known vs ours

- Known: Hansen 1974 (Kerr's moments M(ia)^l); Gibbons–Hawking 1977 (Euclidean Kerr and its helical identification); Newman–Janis shift.
- Ours: reading these through the river's shell rate and Part 8's rules, and finding where the 1PN rules stop.

## Checks (checks.py, 5/5)

| Check | What it tests | Result |
|---|---|---|
| H1 | spin drag term at a = iα is the dipole of a mass moved by −α | pass |
| H2 | second-order term needed; missing from 1PN rules; exact shift has it | pass |
| H3 | Kerr horizon frame-drag rate = river shell rate a/(r₊²+a²), every angle | pass |
| H4 | β × river shell rate at r₊ = 2π tanh u | pass |
| H5 | far swirl continues to real 2Mα/r³; β only from the horizon | pass |

Figure: figures/c2e_one_body_test.png (a: moments, exact shift vs 1PN rules; b: turn per loop for each shell, helix at the horizon).
