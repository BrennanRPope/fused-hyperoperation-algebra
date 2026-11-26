# DRAFT · Q#'s body, Adjoint and Controlled on our shapes

**DRAFT, Nov 26 2026.** Part of Hypercube Composition Theory (repo: fused-hyperoperation-algebra).
Asked for by Brennan (his "Yes." to mapping Q#'s operation forms onto our shapes). Ideas from the voice call, tested
here as hypotheses: single-qubit operations live on the 3-sphere S³, Adjoint is the "mirror point", and Controlled
acts on one Hopf fiber.

Builds on `hypersphere/` (i, π and the sphere of i's on S³, plane dial = 1/R², the swap rule as curvature) and
`polar-recast/` (corner weights, the right shape per family). Checks: `checks.py` → `checks_output.txt`,
**37 of 37 pass** (numpy only). Pictures: `figures.py` → `figures/`. Colour key: **orange = known**, **blue = ours** (a
reading, not a new theorem), **red = does not hold**. "Proved" below means a short algebraic proof written here and
checked numerically; no Lean was compiled (no Lean on this machine).

## The short answer

![Sorting board](figures/fig9-sorting-board.png)

*Figure 9. Every claim tested, where it comes from, and whether it holds.*

- **The body is our core law on a sphere.** Single-qubit gates are the sphere |1 + kx| = 1 (centre −1/k, radius
  1/|k|), and composing gates is exactly a ∘ₖ b there, for every k ≠ 0. The dial k is that sphere's curvature
  (d = k²); as k → 0 the sphere flattens onto the plane of Hamiltonians, where ∘₀ is plain adding.
- **Adjoint is our ∘ₖ inverse, and it is a mirror.** On the gate sphere the ∘ₖ inverse of x is its conjugate x̄, the
  reflection through the real axis. It is **not** the antipode: −U is the same gate up to a phase of −1. The mirror
  reverses order, which is why Q# builds an Adjoint by running the steps backwards.
- **Controlled is a cube corner.** In our coordinates Controlled is x ↦ P1 ⊗ x. It survives ∘ₖ only because
  P1·P1 = P1, which is why Q# can push Controlled into every step. Several controls light one corner of the control
  cube; "control on a OR b" is our union rule ∘₋₁.
- **They commute because the corner is a mirror-symmetric projector** (P1† = P1). The two properties do separate jobs:
  P1·P1 = P1 gives "distribute", P1† = P1 gives "commutes with Adjoint" (both shown failing one at a time).
- **"Controlled acts on one Hopf fiber": no, as stated.** It holds only when the control is surely 1. What is true and
  standard is the reverse direction: Controlled turns a move along the target's phase circle into a visible turn of
  the control (phase kickback).
- **The Clifford gates are corners of our shapes on S³:** the Paulis are the cross-polytope's 8 corners, the
  axis-cycling Cliffords the tesseract's 16, and H, S and friends the 24 corners of the dual 24-cell. T sits on none.

Most of this is standard quantum computing in our language. What is ours is the dictionary: body = ∘ₖ, Adjoint = ∘ₖ
inverse, Controlled = multiply by a corner projector, the dial k = curvature between Hamiltonians and gates.

---

## 1. The body: gates are the sphere |1 + kx| = 1

![Gate sphere](figures/fig1-gates-are-the-sphere-1-plus-kx.png)

*Figure 1. A slice through the real axis and one imaginary direction. Each round sphere is the set of gates for one k;
all pass through the identity at x = 0. As k → 0 they flatten onto the orange plane Re x = 0, the Hamiltonians.*

Match quaternion units to qubit gates by i ↔ −iσx, j ↔ −iσy, k ↔ −iσz. Then unit quaternions multiply exactly as
SU(2) gates (check B1). A Bloch turn by θ about axis n is the quaternion cos(θ/2) + sin(θ/2) n, at distance θ/2 from 1
on S³; a full 2π turn lands on −1, distance π (B4). This is the hypersphere page's "π and spin" again.

**In our coordinates.** The lift u = 1 + kx turns ∘ₖ into multiplication. A gate is a lift with |u| = 1, so in x the
gates fill the sphere |x + 1/k| = 1/|k|. Composing two gates is a ∘ₖ b, and the result stays on the sphere because
sizes multiply (B2, for k = ±1, 0.5, 2, 0.1). So **the body is ∘ₖ, restricted to one sphere.**

**The dial is curvature.** That sphere has radius 1/|k|, so by the hypersphere page's rule its plane dial is
**d = k²**. As k → 0 it flattens onto Re x = 0, and x = (e^{kh} − 1)/k tends to the generator h itself. There ∘ₖ
becomes plain adding: the error falls exactly in half each time k halves (B3). So one dial runs from **Hamiltonians
that add** (k = 0, flat) to **gates that multiply** (k ≠ 0, round).

**Known:** SU(2) = unit quaternions = S³; the Lie algebra is the flat tangent space at 1 (İnönü–Wigner contraction is
the same flattening). **Ours:** reading that flattening as our dial k going to 0. **Status: proved.**

## 2. The Clifford gates are corners of our shapes

![Clifford corners](figures/fig2-clifford-gates-are-corners-of-our-shapes.png)

*Figure 2. The same 48 points of S³ in one fixed view, split into three shapes. Orange: the shape named above each
panel. Grey: the other two.*

Generating the single-qubit Clifford group from H and S gives 48 points of S³ (C1):

| Gates | Points on S³ | Our shape |
|---|---|---|
| ±1, X, Y, Z (the Paulis) | ±1, ±i, ±j, ±k | the 8 corners of the **cross-polytope** |
| the gates that cycle X → Y → Z | (±1 ± i ± j ± k)/2 | the 16 corners of the **tesseract** |
| H, S, √X and the rest | two of ±1, ±i, ±j, ±k, over √2 | the 24 corners of the **dual 24-cell** |

The cross-polytope and the tesseract together are the 24-cell (the Hurwitz units). The roundness dial on the
polar-recast page ran cross-polytope → sphere → cube; here those exact shapes sit inside the sphere as the gates that
are easy to make. T is 0.39 from the nearest corner (C4), and H with T keeps producing new points (1,398 different
ones from words up to length 14): leaving the corners is what makes a gate set universal.

**Known:** the Clifford group mod phase is the binary octahedral group, whose 48 elements are the vertices of a
24-cell and its dual. **Ours:** naming the pieces with our shapes. **Status: proved (counted).**

## 3. Adjoint: our ∘ₖ inverse, which is a mirror

![Adjoint is the mirror](figures/fig3-adjoint-is-the-mirror.png)

*Figure 3. A: U, its Adjoint (the mirror image through the real axis) and its antipode (red). B: the mirror reverses
order. C: the gates that are their own Adjoint lie on the sphere of i's.*

**Adjoint = ∘ₖ inverse = x̄.** Solve x ∘ₖ y = 0: y = −(1 + kx)⁻¹x. On the gate sphere, (1 + kx)⁻¹ is the conjugate of
1 + kx, and the size condition 2·Re x = −k|x|² turns this into **y = x̄** for every k. Its lift is U† (A1). So the
voice-call idea "Adjoint is the inverse" is our ∘ₖ inverse exactly.

**It is a mirror, not the opposite point.** x ↦ x̄ keeps the real axis and reverses the three imaginary directions:
a reflection of R⁴ (determinant −1) whose only fixed gates are 1 and −1 (A2). The antipode −U is something else: it
moves every state exactly as U does, because −1 is only a phase (A3). So if "mirror point" meant the opposite point,
that part is **no**.

**The mirror reverses order.** conj(x ∘ₖ y) = ȳ ∘ₖ x̄, and the same order misses (A4). So the Adjoint of a
five-step body is the five Adjoints run backwards (A5): this is Q#'s `adjoint invert` (the default), derived from
the core law.

**Which gates are their own Adjoint?** Up to phase, exactly the points with Re x = 0: the sphere of i's from the
hypersphere page, the π-turns. There the mirror and the antipode agree. X, Y, Z and H are on it; S and T are not
(A6). These are the gates Q# marks `adjoint self`.

**Known:** U† = U⁻¹, the quaternion conjugate, the order reversal. **Ours:** that this is our ∘ₖ inverse for every k,
and the shape reading. **Status: proved.**

## 4. Controlled: multiply by a corner of the control cube

![Controlled lights one corner](figures/fig4-controlled-lights-one-corner.png)

*Figure 4. A: with three controls, U acts only at corner 111. B: the controlled gate as blocks. C: "a OR b" lights
three corners, and the projector for it is our union rule.*

Controlled U = 1 + P1 ⊗ (U − 1), with P1 = |1⟩⟨1|. Subtract 1 and divide by k: **in our coordinates Controlled is
x ↦ P1 ⊗ x**, the same for every k (K2). That keeps ∘ₖ, because

  P1x ∘ₖ P1y = P1x + P1y + k·P1x·P1y = P1(x ∘ₖ y), using P1·P1 = P1.

So Controlled(UV) = Controlled(U)·Controlled(V) (K1), and Q# can build the Controlled form of a body by controlling
each step: `controlled distribute`, derived. Replace P1 by something with e·e ≠ e (half of P1, or a nilpotent) and it
breaks, by 1.8 and 7.3 (K2).

**Several controls are a cube corner.** Two controls are 1 + P1 ⊗ P1 ⊗ (U − 1) (K3). The 2ⁿ control patterns are the
corners of the control cube: their projectors never overlap and add to 1, which is the cube's own expansion
(P0 + P1)ⁿ (K4). "Control on zero" picks another corner. "Control on a OR b" uses P_a ∘₋₁ P_b =
P_a + P_b − P_aP_b = 1 − P0 ⊗ P0, our union rule, and it is again a projector (K5).

**Known:** the block form, controlled-on-zero, OR built from projectors. **Ours:** Controlled as "multiply x by a
corner projector", and OR as ∘₋₁. **Status: proved.**

## 5. Why Adjoint and Controlled commute

![Commuting square](figures/fig5-adjoint-and-controlled-commute.png)

*Figure 5. A: the square closes. B: which property of the corner does which job.*

(1 + P1 ⊗ (U − 1))† = 1 + P1 ⊗ (U† − 1), because P1† = P1. So Controlled Adjoint = Adjoint Controlled (D1). It does
follow from our algebra, and it needs two separate facts about the corner:

- **P1·P1 = P1** gives "distributes over a body".
- **P1† = P1** gives "commutes with Adjoint".

Tested one at a time: an oblique projector E (E·E = E, E† ≠ E) still distributes but stops commuting with Adjoint and
stops being a gate (D2); half of P1 (Hermitian, but e·e ≠ e) commutes with Adjoint but does not distribute (D3).
Q#'s control needs both, and that is exactly an orthogonal projector.

**Known:** the commuting is standard. **Ours:** splitting it into the two corner properties. **Status: proved.**

## 6. Controlled and the Hopf circle: testing the hypothesis

![Controlled and the Hopf circle](figures/fig6-controlled-and-the-hopf-circle.png)

*Figure 6. A: an eigenstate of U only slides round its phase circle. B: the control turns by that angle. C: for
random targets the control's Bloch vector shrinks, so this is not "one fiber".*

The Hopf circle over a qubit state is its phase. A global phase moves a state only round that circle and is invisible
(H1). Controlled makes it visible: Controlled(e^{iφ}·1) is a turn of the **control** about Z by φ (H2). If the target
is an eigenstate of U, the target only slides round its own circle and that angle is handed to the control (H3):
phase kickback, the engine of phase estimation.

**The hypothesis as stated fails.** For a random target, branch 1's Bloch point moves (not just its phase) in all
2,000 trials, and the control ends up entangled (mean Bloch length 0.66 instead of 1) (H4).

![Quaternionic Hopf base](figures/fig7-two-qubits-on-the-quaternionic-hopf-base.png)

*Figure 7. Two qubits on S⁷, drawn on the quaternionic Hopf base S⁴. Height is the control's population difference.*

**Where it does hold.** Write two qubits as a pair of quaternions (q0, q1) on S⁷ (with C² = H, every gate is left
multiplication by a unit quaternion, H5). The quaternionic Hopf map sends S⁷ to S⁴ with S³ fibers. Controlled-U
sends each fiber to a fiber, keeps the height |q0|² − |q1|², and turns the rest of S⁴ about that axis (H6). At the
pole where the control is surely 1, it acts entirely inside one S³ fiber (H7). So "acts on one Hopf fiber" is exactly
true there and nowhere else.

**Known:** kickback (Cleve, Ekert, Macchiavello and Mosca 1998); the quaternionic Hopf fibration. **Ours:** the S⁴
picture of Controlled as a turn about the control axis. **Status: "one fiber" no as stated; kickback and the S⁴ turn
proved.**

One structural point falls out: the body and Adjoint only see S³ (phase does not matter), but Controlled needs the
full U(2), the S³ plus its phase circle. That is why Q# keeps track of global phase at all.

## 7. The swap rule for gates

![Swap rule](figures/fig8-swap-rule-for-gates.png)

*Figure 8. A: turn about x, about y, undo both: a turn about z is left. B: leftover ÷ area goes to 1.*

Rx(a) Ry(b) Rx(−a) Ry(−b) leaves a turn about Z of size a·b (W1: ratio 0.99999 at small sizes). That is the swap rule
"turn = dial × area" on the Bloch sphere. **Known** (the Lie bracket; it is the trick inside the Solovay–Kitaev
algorithm). **Status: proved.**

**In our law the order cost is exact:** x ∘ₖ y − y ∘ₖ x = k[x, y] (W2). And splitting a Hamiltonian A + B into a
step of A then a step of B, each of size k, misses the true gate by (k²/2)·|[A, B]| (W3: ratio 1.0000). Read on the
gate sphere, the step size k is the dial, k² is the sphere's curvature d, and the miss is d times the swap-rule
commutator. **Known:** first-order Trotter error. **Ours:** that it is the swap rule with d = k². **Status: proved.**

---

## Leads, judged on their own merits

1. **The dial k as the knob between Hamiltonians and gates (sections 1 and 7).** Worth following. One parameter
   turns "add generators" into "multiply gates"; it is the time step, its square is the gate sphere's curvature, and
   the first-order cost of splitting a Hamiltonian (Trotter error) is that curvature times the swap-rule commutator.
   The next test has a definite answer either way: do the higher Trotter orders (the k³ and k⁴ terms of
   Baker–Campbell–Hausdorff) come out as higher rungs or corner terms of our law?
2. **Controlled as "multiply by a corner projector" (sections 4 and 5).** Worth following for its reach, not for new
   physics. It says which algebraic facts any "controlled" construction needs (e·e = e for distributing, e† = e for
   Adjoint), and it ties multi-controls and OR-controls to the cube and the union rule we already have. It would make
   a clean Lean file in our library.
3. **Clifford gates as corners of our shapes (section 2).** A good picture and a good teaching hook; known, no new
   result.
4. **The Hopf-fiber reading of Controlled (section 6).** Not worth pursuing as stated. The true parts are standard.

## Next (offered, not started)

- Test lead 1 further: the k³ and k⁴ Trotter terms against our rungs and corner terms.
- Write the Lean file for sections 3 to 5 (∘ₖ inverse on the gate sphere, order reversal, idempotent homomorphism,
  the commuting square) on the library's ring structures.

## Files

- `README.md`: this write-up
- `checks.py` → `checks_output.txt` (37/37)
- `figures.py` → `figures/fig1` … `fig9`
- `page_build.py` → `index.html`

## References

- W. R. Hamilton, quaternions, 1843; the SU(2) ≅ S³ correspondence in any quantum mechanics text.
- H. S. M. Coxeter, *Regular Polytopes*, 1973 (the 24-cell and the binary polyhedral groups).
- D. Gottesman, "The Heisenberg representation of quantum computers", 1998 (the Clifford group).
- R. Cleve, A. Ekert, C. Macchiavello, M. Mosca, "Quantum algorithms revisited", Proc. R. Soc. A 1998 (phase kickback).
- E. İnönü, E. P. Wigner, "On the contraction of groups and their representations", PNAS 1953.
- C. M. Dawson, M. A. Nielsen, "The Solovay–Kitaev algorithm", 2005.
- Microsoft, Q# language documentation: functors Adjoint and Controlled, specializations `adjoint invert`, `adjoint self`, `controlled distribute`.

Page: https://claude.ai/artifact/VeVmLnshwkAUSvqFcvXof4 (slug e7fc20f4-1f9e-4a5e-bdfa-fffe89815177, version 1795656952-2106).
Not yet pushed to the repo (BrennanRPope/fused-hyperoperation-algebra): the repo thread does that.
