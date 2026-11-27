# Do rotations force quantum mechanics to use complex numbers? (Nov 26, 2026)

**Verdict.** Rotations on their own don't force it, because real-number quantum mechanics also evolves by rotations. Add two textbook facts, energy = ħ × turning rate and "only energy differences matter", and a five-line proof follows: every Hamiltonian must turn states with one and the same quarter turn J (J² = −1). That J does the job of i in Schrödinger's equation, and it is shared by the whole world, which rules out real quantum mechanics with independent parts.

Asked by Brennan (Nov 26): "I really want to prove quantum mechanics requires complex numbers because the steps along the way are rotations."

Method: published work first. Anything not found in the sources is labeled *ours*.

Checks: `checks.py`, 22/22 pass (`checks_output.txt`). Figures: `figures.py` writes `figures/r1`–`r6`.

---

## 1. What rotations give you, and what they don't

![Three panels. Left: a single real energy level is a line, where the state can only stay, stretch or flip, so it can't turn steadily. Middle: a level that is a plane turns at rate E/ħ. Right: the turn of that plane is a quarter turn J with J squared = −1.](figures/r1-a-rotation-needs-a-plane.png)

A step that keeps the size of the state is a rotation. Its generator A is antisymmetric (Aᵀ = −A). Three facts follow:

- **In odd dimension every rotation has a frozen direction** (det A = 0; check T1). Steady turning needs each energy level to be a plane.
- **Each turning plane has its own quarter turn** with frequencies in ± pairs (T1b).
- **For one dynamics, rotations always supply a J.** Take J = A|A|⁻¹ and H = ħ√(−A²); then A = JH/ħ (T2). This is the positive-frequency complex structure of Segal and Ashtekar–Magnon, recently re-derived for this question by Aste (2019) and Vaughan (2026).

That last point is why rotations alone are not a proof. **Every** rotation dynamics, including real quantum mechanics, comes with a J of its own. The catch is the words "of its own": in real quantum mechanics different Hamiltonians get different J's (counterexample C2, two Hamiltonians on ℝ⁴ whose quarter turns differ). Then there is no single i in the theory.

So the real question is: **what forces every Hamiltonian to use the same J?**

## 2. The energy-zero lever (the proof)

![Left: three boxes. A Hamiltonian H turns with its own J_H; shifting the zero, H + c, adds a common turn cJ; steady turning for every c forces J_H = J. Right: a plot of how badly steady turning fails against the shift c. With a shared turn the failure is zero for every c; with two different quarter turns it grows linearly as soon as c moves off zero.](figures/r2-the-energy-zero-lever.png)

### Setup and assumptions

States are real unit vectors in ℝⁿ. Observables are symmetric matrices, and any observable can serve as the energy of some arrangement (A0). Each Hamiltonian H drives the state by ψ' = A(H)ψ, where:

| | Assumption | Where it comes from |
|---|---|---|
| A1 | A(H) is antisymmetric: **every step is a rotation** | Brennan's premise; continuous reversibility (Hardy; Masanes–Müller) |
| A2 | A(H) commutes with H: **energy is conserved** | Noether; Stone. **Not needed: it follows from A3 and A4 (step 3)** |
| A3 | A(H)² = −H²/ħ²: **a state of energy E turns at rate E/ħ** | Planck–Einstein |
| A4 | A(H + c) = A(H) + c·A(1): **only energy differences matter** | the zero of energy is a convention |

### Proof

1. **A3 at H = 1:** J := ħA(1) satisfies J² = −1, and it is antisymmetric by A1. So J is a quarter turn, and n is even.
2. **A3 for H + c:** with B := ħA(H), A4 gives ħA(H + c) = B + cJ, and A3 gives (B + cJ)² = −(H + c)². Matching powers of c gives B² = −H² and BJ + JB = −2H.
3. **Energy conservation comes free.** Because J² = −1, J commutes with BJ + JB, so **J commutes with H**. Then B commutes with BJ + JB too, because B²J − JB² = −H²J + JH² = 0, so **B commutes with H**: that is A2, derived rather than assumed (check T3d; Lean `EnergyShift.A2_redundant`).
4. **Write B = JH + E.** JH is antisymmetric, because J and H commute. So E is antisymmetric, anticommutes with J and commutes with H. Then B² = −H² + E², because the cross terms cancel (EJH = −JEH = −JHE). So **E² = 0**.
5. **An antisymmetric E with E² = 0 is zero**, because |E|² = tr(EᵀE) = −tr(E²) = 0 (T3b).

So **A(H) = JH/ħ for every Hamiltonian, with one and the same J.** Schrödinger's equation reads ħψ' = JHψ, which is iħψ' = Hψ with i = −J.

Checks:
- **T3:** on ℝ⁴, with all six entries of B left unknown, B = JH is the only solution.
- **T3a:** the degenerate case H = λ·1 leaves a two-parameter family after the linear steps. The Planck–Einstein condition then gives −(b₃² + b₄²) = 0, which kills it over the reals.
- **T3c:** a perturbation E that keeps A1 and A2 breaks A3. Two states of the same energy would then turn at different rates.

### What it costs real quantum mechanics

Since J commutes with every Hamiltonian (step 3), A0 makes it commute with every observable. In real quantum mechanics every symmetric matrix is an observable, and only multiples of 1 commute with all of them (T4). No quarter turn can do that, so **real quantum mechanics fails A4**: shifting the zero of energy would change the physics.

## 3. Real, complex or quaternion

![Three boxes. Real numbers: 0 shared quarter turns, fails A4. Complex numbers: 1 shared quarter turn J, meets A1 to A6, ordinary quantum mechanics. Quaternions: 3 quarter turns i, j, k, meets A1 to A4 with one of them, but hides observables and composites.](figures/r4-real-complex-quaternion.png)

Once one shared J exists, the observables are symmetric matrices that commute with J, which are exactly the complex Hermitian matrices (T5). Frobenius (1878) leaves three candidate number systems, and quaternions fail in one of two ways:

- **Genuine quaternionic QM (Adler 1995) fails A4.** Its dynamics must respect all three quarter turns i, j, k. The only operators commuting with every quaternionic observable are built from those three turns, and none of them respects the other two. So there's no J to carry the energy shift, and the energy zero isn't a convention there (C3b; caught in the rigor pass on Nov 27).
- **A variant that picks J = R_i meets A1–A4 but hides observables.** It allows 6 observables on a quaternionic two-level system, against 16 symmetric operators that commute with J (C3). That's where one more assumption comes in:
- **A5, nothing hidden:** every symmetric operator commuting with J is an observable. This is Masanes–Müller's "all effects are allowed" in this setting. With A5, the observables are the complex Hermitian matrices, (n/2)² of them (T5), and the theory is complex quantum mechanics.

Vaughan (2026) reaches the shared-structure point from normalisation and energy conservation: complex structure "at least", with quaternions not excluded. A4 and A5 close that gap here.

## 4. Two parts share one i

![Left: part A turns with J_A, part B with J_B. Add c to the total energy: booked to A it turns by J_A, booked to B by J_B, and since it is the same physics, J_A ⊗ 1 = 1 ⊗ J_B. Right: bars of joint parameters against products of local observables. Real product: 10 against 9, one joint turn hidden. Shared turn: 16 and 16. Quaternionic: 28 against 36, too many, no composite fits.](figures/r3-composites-share-one-i.png)

The lever works for composites too:

- **A6, a total shift has no owner.** Adding c to the total energy of A and B can be booked to A, which turns the joint state by J ⊗ 1, or to B, which turns it by 1 ⊗ J. It's the same physical change, so the two must agree on every joint state: J ⊗ 1 = 1 ⊗ J.

That condition cuts the plain real product ℝ^{n₁} ⊗ ℝ^{n₂} down to half its size (T6). On what remains:
- J ⊗ J = (J ⊗ 1)² = −1, so the whole system has one i.
- The space is the complex tensor product, which is Bang–Cho–Baek's "balanced tensor product" (2026).
- Products of local observables span every joint observable: 1/1, 16/16 and 36/36 for 1×1, 2×2 and 2×3 (T6). This is **local tomography**, the property that Hardy and Masanes–Müller use to single out complex numbers.

Without A6 (C4), on the plain real product the two bookings differ and J ⊗ J is a free hidden parameter. That is the extra joint quarter turn on the e-i-pi page (σ_y ⊗ σ_y). The 2021 Renou test and the 2022 experiments rule out exactly that version of real quantum mechanics.

## 5. What each assumption does

![Table of assumptions A1 to A6 with where each comes from and what goes wrong without it. A1 to A3 (blue) give a quarter turn; A4 to A6 (orange) make it one shared i.](figures/r5-assumptions.png)

- **A1–A3 give a J for each dynamics.** This part is known and is not a proof: real quantum mechanics passes it.
- **A4 makes the J shared by every Hamiltonian.** This is the step Brennan's idea needed. Without it you get real quantum mechanics (C2, T4).
- **A4 also rules out genuine quaternionic QM** (C3b), and **A5 rules out the quaternionic variant** (C3).
- **A6 makes one J for the whole world** and gives local tomography (T6, C4).

## 6. Where this sits in the 2021–2026 debate

![Timeline. 2021 Renou et al.: real QM with the plain real product can be falsified. 2022 experiments side with complex QM. 2025 Hoffreumon–Woods: real QM with a new combination rule matches everything. 2026 Barrios Hita et al.: same claim, quotient-space construction. 2026 Bang–Cho–Baek: those real theories carry a hidden J; Vaughan: evolution needs at least a J. 2026 this page: rotations, E = ħω and only energy differences matter force that one shared J.](figures/r6-the-debate.png)

- **Renou et al. (2021), Chen et al. and Li et al. (2022):** real quantum mechanics with the plain real product is falsified.
- **Hoffreumon–Woods (2025) and Barrios Hita et al. (2026):** real theories with a modified composition rule match every experiment.
- **Bang–Cho–Baek (2026):** the Barrios Hita construction needs a hidden J with J² = −1 and a balanced tensor product, so it is complex quantum mechanics in real notation.
- **Feng–Ren–Vedral (2025):** a real theory with independent sources needs a nonlocal map.

**Tested against the theorem (Nov 27, checks R1–R4).** All three real formulations that match experiment, Hoffreumon–Woods, Barrios Hita et al. and Maioli–Curado–Gazeau (arXiv:2604.19482), meet every assumption once their dynamics is taken as the image of the complex one. None of them defines a Hamiltonian flow, so that step is our reading.
- **Hoffreumon–Woods:** their observables are exactly the symmetric matrices commuting with their J (their Lemma B.6), and their product rule equals projection onto W (R1, R2).
- **Barrios Hita et al.:** their quotient throws away exactly the J⊗J = +1 half, so their composite is W (R3).
- **Maioli–Curado–Gazeau:** their product rule is the real form of the complex tensor product (R4).

So the shared J isn't something these authors chose. The theorem says any theory meeting the assumptions must have it, and the one real theory without it is the one the experiments ruled out.

This proof says that whatever notation you use, a theory with rotation steps, E = ħω and a free energy zero has one shared quarter turn. Whether a real notation that carries that J "counts" as real quantum mechanics is a question about names, not physics. The physics is the shared J, and it is forced.

## 7. In the framework's terms

The dial's step is 1 + (t/N)·k·H. A real k stretches the state at first order. An imaginary k keeps its size to second order and builds exp(tJH/ħ) (T7; see the chain-resolves page). The theorem adds one thing: **the dial's turning direction is one J, the same for every Hamiltonian and every part of the world.** There is a single imaginary dial, k = −i/ħ, not one per system.

## 8. What's known and what's ours

| Piece | Status |
|---|---|
| A rotation dynamics carries its own J (positive frequency) | known: Segal; Ashtekar–Magnon; Aste 2019; Vaughan 2026 |
| A universal J commuting with all observables | postulated by Stueckelberg (1960); derived from Poincaré symmetry by Moretti–Oppio (arXiv:1611.09029) |
| Complex structure required for energy-conserving, normalised evolution | known: Vaughan 2026 |
| Balanced (complex) tensor product inside real constructions | known: Bang–Cho–Baek 2026 |
| Local tomography picks ℂ | known: Hardy 2001; Masanes–Müller 2011 |
| **Energy-zero lever: A4 plus Planck–Einstein forces one J for all Hamiltonians (§2)** | *ours*, as far as I found |
| **"A total shift has no owner" gives the balanced product (§4)** | *ours*, as far as I found |

I didn't find the energy-zero argument in what I read. It is short enough that it may exist elsewhere, perhaps in older work on real Hilbert spaces or in the Moretti–Oppio line. An expert should be asked before it's called new.

## 9. Is it a lead?

**Yes, as a clean proof. It is not new physics.** It turns Brennan's premise into a theorem with one honest correction: rotations get you a quarter turn, and "only energy differences matter" makes it one shared i. The argument is five lines, uses only textbook assumptions, and lands on the experimentally tested property (local tomography). It is worth writing up as a short note once someone who knows the real-Hilbert-space literature has checked that the energy-zero step isn't already published.

## Checks

`python3 checks.py` gives 21/21 PASS (sympy and numpy):

- Q0: complex QM satisfies A1–A4.
- T1, T1b: odd dimension; turning planes.
- T2: positive-frequency J.
- C2: two Hamiltonians with different J's.
- T3, T3a, T3b, T3c: the proof, generic and degenerate.
- T4: real QM fails A4.
- C3, C3b, T5: quaternions (variant and genuine); complex Hermitian count.
- T6, C4: composites and local tomography.
- T7: the dial chain.
- R1–R5: the three recent real formulations tested against every assumption (§6), and why A4 has content.

## Assumptions

- The setting is finite-dimensional, with states as real unit vectors and dynamics linear and generated by A(H).
- A0, that any observable can be some Hamiltonian, is needed for §2's conclusion about observables. Without it, J only commutes with the Hamiltonians.
- A4 is read as A(H + c) = A(H) + c·A(1), so a shift of the energy zero adds a fixed turn and nothing else.
- The symbolic check T3 is on ℝ⁴ only. The proof in §2 holds in any finite dimension, and its last step is checked in general by T3b.
- The algebraic core is in `lean/EnergyShift.lean` (18 theorems in core Lean, any ring, no matrices needed; a copy sits in `arxiv/anc/`). It compiled clean on Nov 27 in the shared library as `FusedAlgebra.Topics.EnergyShift`. Step (v), E antisymmetric with E² = 0 ⇒ E = 0, uses the transpose and trace, so it stays with the paper and the sympy checks.
- The literature placement comes from abstracts and the papers I opened. Vaughan's date, and whether the energy-zero lever is already published, aren't confirmed.

## Files

- `checks.py`, `checks_output.txt`
- `figures.py`, `figures/r1`–`r6`
- `NOTE-DRAFT.md`: a short note on the proof for an expert to check (draft, not sent)
- `arxiv/`: the proof as an arXiv paper (`main.tex`, figures, ancillary checks, `SUBMIT-STEPS.md`); draft, not submitted

## Sources

- Stueckelberg, "Quantum theory in real Hilbert space", *Helv. Phys. Acta* 33, 727 (1960); Moretti and Oppio, [arXiv:1611.09029](https://arxiv.org/abs/1611.09029)
- Aste, "Origin of the complex structure of quantum mechanics", [arXiv:1905.12894](https://arxiv.org/abs/1905.12894) (2019)
- Vaughan, "The requirement of (at least) complex structure for quantum mechanics", [arXiv:2609.23041](https://arxiv.org/abs/2609.23041) (2026)
- Bang, Cho and Baek, "Hidden complex structure in quotient-space real quantum mechanics", [arXiv:2607.05865](https://arxiv.org/abs/2607.05865) (2026)
- Barrios Hita et al., *PRL* 136, 240202 (2026), [arXiv:2503.17307](https://arxiv.org/abs/2503.17307)
- Hoffreumon and Woods, [arXiv:2504.02808](https://arxiv.org/abs/2504.02808) (2025)
- Feng, Ren and Vedral, "Locality implies complex numbers in quantum mechanics", [arXiv:2504.07808](https://arxiv.org/abs/2504.07808) (2025)
- Renou et al., *Nature* 600, 625 (2021); Chen et al., *PRL* 128, 040403 (2022); Li et al., *PRL* 128, 040402 (2022)
- Hardy, quant-ph/0101012 (2001); Masanes and Müller, *NJP* 13, 063001 (2011)

Textbook, not re-fetched: Frobenius (1878); Segal (1950s, positive-frequency complex structure); Ashtekar and Magnon (1975); Adler (1995, quaternionic QM).
