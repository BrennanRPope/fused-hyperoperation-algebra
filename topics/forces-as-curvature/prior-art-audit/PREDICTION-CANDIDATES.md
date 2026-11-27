# DRAFT · Prediction candidates across the whole theory

*DRAFT, Nov 27 2026, from the prior-art audit (`README.md` in this folder). This is a list, not a test. The forces thread's
rounds own the testing: round 27 (`../round27/`) worked gravity candidates G1–G6, round 28
(`../round28/`) worked C2 and C4, and round 29 takes C6.*

**What counts as a candidate:** a number or result that the law fixes, where standard physics or mathematics leaves it
free, fits it from data, or gives a different value. Each row says where it sits (file and line), whether the law
forces it or it is only chosen, whether it's already known or measured, and what would decide it.

Status key: **forced** = follows from the law with no extra choice; **chosen** = needs an input the law doesn't supply;
**ruled out**, **open**, **agrees (textbook)**, **math check** (a yes/no a calculation settles, not a measurement).

## A. Physics

| # | Candidate | Where | Forced or chosen | Known or measured | What would decide it |
|---|---|---|---|---|---|
| C1 | Λ = 3q with q the curvature of today's space, so Ω_k = −Ω_Λ ≈ −0.69 | `forces-as-curvature-round20/README.md:24`; `-round16/README.md:36` | Forced if read literally | **Ruled out** (round 27 G1): Planck + BAO give Ω_k = 0.001 ± 0.002 | Done; only holds at a moment of rest (round 27 G2) |
| C2 | Space is closed (Ω_k < 0) if matter's rest slices are the law's 3-spheres | round 27 G3 | Chosen (de Sitter can be sliced flat or open) | **Open, but not a test of the law.** DESI DR2 + CMB lean open at 2.1σ. Round 28 S1–S3: the law doesn't force the sign; closed space needs an added premise, and a closed result would fit ordinary closed ΛCDM too | **Worked in round 28: no lead.** Data still decide the premise |
| C3 | Λ's size from the number of 600-cell refinements (~202 halvings from Planck size) | `forces-as-curvature-round16/README.md:36` | Chosen: the count is an input. 202 ≈ log₂ 10⁶¹ restates Λ in Planck units | Not a prediction | Only if something in the law fixed a smallest kick |
| C4 | Spin–twist coupling fixed by the law at a value different from Einstein–Cartan's | `forces-as-curvature/README.md:24` | Nothing to fix (round 28 S4) | The law's connection is torsion-free, so its twist is curvature, not torsion, and it adds no number to Einstein–Cartan's coupling; the twisted members give ~10⁻⁴³ GeV against bounds near 10⁻³¹ GeV | **Worked in round 28: no lead** |
| C5 | κ-Poincaré: dial k = −1/E_Planck, wall at the Planck energy | `why-i/kappa-poincare/README.md:7` | Chosen: E_Planck is put in | Known (Lukierski et al. 1991; Judes–Visser). Whether it gives photon time delays depends on the momentum basis | Fermi GRB 090510 and LHAASO GRB 221009A time-delay bounds already exceed E_Planck for a linear delay |
| C6 | The dial as rolling dark energy; rescaling forces an exponential potential in ln k; beats Λ by Δχ² ≈ 4.85 | `dial-potential/README.md:11` | The exponential form is forced by the rescaling assumption; the slope is fitted | Exponential quintessence is known (Wetterich 1988; Ratra–Peebles 1988). The hint matches DESI 2025 | DESI DR3 and Euclid w₀–w_a; a cosmological constant at > 3σ rules it out. **Taken up as round 29** |
| C7 | Higgs mass 128.7 ± 1.0 GeV from a dial flat at the Planck scale | `dial-potential/README.md:10` | Forced within that reading | Known idea (Froggatt–Nielsen 1995; Shaposhnikov–Wetterich 2009). Measured 125.2 GeV, 3.3σ off | **Disfavoured** |
| C8 | Corrected 2PN term: two comparable one-charge dilaton black holes precess at fourth order in G, unlike half-BPS pairs in N = 8 supergravity | `hypersphere/README.md:363` | Forced by the corrected Lagrangian (not by the law) | Open in the literature, as far as the page found | A theory check: an amplitude or post-Newtonian calculation in N = 8 supergravity |
| C9 | Newton-G ratio ρ: 1/3 against gravity's 3/5 | `newton-g/README.md:78` | Chosen: needs a clock (base dial r) the law doesn't fix | Sakharov induced gravity is known | Not decidable until something fixes the clock |
| C10 | Gluon colour factor C_A = 3 from the hidden space's curvature | `forces-as-curvature-round5/README.md:9` | Forced once SU(3) is chosen | **Agrees, textbook** (LEP: 2.93 ± 0.60, 3.02 ± 0.55) | Already decided |
| C11 | Law's f = 1 twist filling a closed universe; lowest gravitational wave and couplings (3, −6) | round 27 G4, G5 | Forced | Fixed but untestable: ≤ 4e-126 of critical density; wavelength ≥ 41 Hubble radii | Out of reach |
| C12 | A qubit's two filtering measurements leave a Bloch turn equal to their triangle's area, sign set by d = −1 | `bridge-survey/README.md:126` (P3) | Forced within standard quantum mechanics | Unchecked; would test known geometric-phase physics, not the law | Calculation first, then a lab test with weak measurements |
| C13 | The dial boost covers barriers and mixers | `tunnelling/README.md:206` | Forced within standard wave optics | Unchecked against optics literature | Multilayer or interferometer measurement; likely standard transfer-matrix physics |

## B. Mathematics (sharp yes/no checks, not measurements)

| # | Candidate | Where | Status | What would decide it |
|---|---|---|---|---|
| C14 | Critical Poisson(1) survival gets a next term −(2/3)(log n)/n² (general: −(2·residue/σ²)(log n)/n²) | `parabolic-fingerprints/README.md:90` | Forced. The order is in print as an error bound; the constant wasn't found | Read Zolotarev 1957 and Jurayev 2020; exact iteration to large n |
| C15 | Hirota soliton corner weights pass the holographic monogamy test | `bridge-survey/README.md:124` (P1) | Open | Check I₃ ≤ 0 for S(A) = log det(I + M_A), all positive Cauchy matrices |
| C16 | Kerr's ISCO factors like Schwarzschild's | `bridge-survey/README.md:125` (P2) | Open | Factor r² − 6Mr ± 8a√(Mr) − 3a² |
| C17 | Free-fermion entanglement glues by ∘₋₁ at every cross-ratio | `bridge-survey/README.md:127` (P4) | Open | Two-interval free-fermion entropy formula |
| C18 | κ-Minkowski spectral dimension is a two-wall curve | `bridge-survey/README.md:128` (P5) | Open | Fit the published curve |
| C19 | The AdS₃ cross-ratio clock fails in excited states | `bridge-survey/README.md:66` (N4) | Open; falsifiable | CFT excited-state entanglement computation |
| C20 | Cooperative-binding zeros on a Lee–Yang circle, radius ω^(−(n−1)/2); concerted (MWC) binding on another circle | `bridge-survey/README.md:101` (G5) | Forced for the pairwise model | Fit published hemoglobin Adair constants; read Briggs 1985 first |
| C21 | Three-body period from a five-number run-price law, 0.28% typical error | `three-body-word-formula/README.md:39` | An empirical fit, not from the law | Score it on catalogues published after 2018 |
| C22 | Landscape laws at the colour-count dial values k_n = 1 + 2cos(2π/n) | `lee-yang/README.md:17` | Tested: 0 hits at the 7 non-integer values | Already decided: no |

## Short read

- **No candidate yet is both forced by the law alone and testable.** C1 was forced and is ruled out. C11 is forced and
  out of reach. Everything testable needs an extra choice (C2, C5, C6) or comes from standard physics rather than from
  the law (C10, C12, C13).
- **The two flagged physics items are now worked, with no lead.** Round 28 found the law doesn't force the closed sign
  (C2, S1–S3) and has no torsion, so no spin–twist coupling to compare (C4, S4). The one physics item still open is C6,
  the dial as rolling dark energy, which round 29 takes up.
- **The best mathematical prediction is C14.** It is forced and checkable, but it is a statement about branching
  processes, not physics.
