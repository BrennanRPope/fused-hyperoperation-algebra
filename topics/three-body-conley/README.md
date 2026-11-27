# DRAFT · The Conley index on a real three-body return map

**DRAFT, Nov 27 2026.** Part of Affine Composition Theory, built on Fused Hyperoperation Algebra. Follows
[entropy-topology](../entropy-topology/README.md), which found that the Conley index tells a stable orbit from its
unstable twin where knot type can't, in the logistic, standard and Hénon maps, and graded "test it on a three-body
return map" as worth following (high). Asked for by Brennan on Nov 27 ("Yes" to that follow-up).

Checks: [`checks.py`](checks.py) → [`checks_output.txt`](checks_output.txt), **21 of 21 pass** (numpy, sympy, gcc;
a few seconds; it reads the saved JSON and recomputes the key numbers live). No Lean module. Pictures:
[`figures.py`](figures.py) → [`figures/`](figures/). Colour key where it applies: **orange = known, blue = our reading
of known pieces, green = not found in print** (nothing on this page is green).

## The answer first

**On a real three-body return map, the index tells twins apart but cannot tell stable from unstable.**
That is a clean negative result with a one-line reason, and it marks the road that could work.

![Figure 2](figures/fig2-classes-and-index-pairs.png)

*Figure 2. Left: the seven ways a fixed point of the reduced three-body map can sit, and the fixed-point index of the
map F and of F² for each. The four in the blue box share the index pair (+1, +1); only the top one (E-E, both pairs on
the unit circle) is stable. Right: the multipliers of 21 refined orbits, coloured by class.*

1. **The map is real and checked.** The reduced planar three-body problem (equal masses, zero angular momentum, fixed
   energy) has a four-dimensional section: the line-ups with body 3 in the middle. Its return map F was built and
   checked two independent ways (figure 1; checks A1–A4). The figure-eight is a fixed point with all four multipliers
   on the unit circle (E-E).
2. **In four dimensions the index forgets stability.** For a nondegenerate fixed point the index is
   (−1)^(number of real multipliers above 1). With two Floquet pairs that only counts hyperbolic pairs mod 2, so a
   stable point (E-E) and a fully unstable one (H-H, IH-IH, or a complex quartet) get index +1 **for every iterate**
   (B1–B2). In the plane there is one pair, so the same count sorts everything; that is why it worked last round.
   On the Li–Liao catalog: of the 52 orbits whose type is safe from rounding and whose index is +1 for every iterate,
   only 4 are stable (figure 3; B6).
3. **A value-only count works on the real map.** The index of F^k was computed from map values alone (a simplicial
   degree count in four dimensions) and agrees with the Jacobian on every orbit where it settled, 14 of 15 (figure 4;
   B4–B5). So the index really is computable without derivatives here, the step interval arithmetic would certify.
4. **Twins: the index tells them apart, not which is stable** (figure 5). The figure-eight (index +1) and its twin,
   Simó's H orbit (index −1), straddle the wall trace = 2 in one Floquet pair: a fold, seen in the numbers. In 5 of the
   7 same-word pairs the indices are opposite; in those 5, the +1 member is the less unstable one only 3 times (C1–C4).
5. **The pointwise Conley index exists only for fully unstable orbits.** A fully hyperbolic orbit is an isolated
   invariant set with Conley index a pointed 2-sphere and the map ±1 on it. A point with a generic elliptic pair is
   wrapped in invariant tori or circles, so no small box isolates it and it has no Conley index of its own (figure 6;
   D1–D2). So the Conley index proper can certify full instability; it never certifies stability.
6. **What the dials add: each Floquet pair carries its own curvature dial** (figure 7). Iterating an elliptic pair is
   the tangent law (k = 0, d = +1) on x = tan(θ/2); a hyperbolic pair is velocity addition (k = 0, d = −1) on
   x = tanh(μ/2); squaring the map sends each pair's residue R to 4R(1 − R) (E1–E3). The index is the **product** of
   the pairs' signs of d. In the plane the product is the pair; in the three-body map (+1)(+1) = (−1)(−1). That is the
   whole failure in one line, and it says what a fix needs: the dials pair by pair, not their product.
7. **A correction to an earlier page.** The finite-difference instability numbers in
   [three-body-order](../three-body-order/README.md) (`stability.json`) are too high by more than 1 (in log) for 78 of
   239 orbits (A6); the exact variational numbers here replace them. That page's stability check (H1) should be rerun.

**Is it a lead?** As a way to read stability in three bodies: **no**, and now we know exactly why (item 2 is a
theorem, not a sampling accident). Grades two ways:

- **Worth following: low** for stability by the Conley index or its index shadow; **moderate** for two narrower uses:
  certifying that an orbit is fully unstable (an isolating block with a 2-sphere index), and twin bookkeeping (opposite
  indices mean a pair that can be born or die together).
- **Could reach a proof: medium** for existence and index of a given orbit by interval arithmetic (the value-only count
  in item 3 is the piece that would be certified); **low** for stability, which needs KAM (done for the figure-eight by
  Kapela and Simó).

**What would make it a lead:** a topological invariant that keeps the pairs separate. The known candidate is the
Conley–Zehnder index and its iterates (Bott's formula; Long's index iteration theory): an elliptic pair adds a
rotation that grows irregularly with the iterate, a hyperbolic pair a fixed even or odd amount, so E-E and H-H part
ways. It is computable rigorously (Joung and van Koert 2024 do it for the restricted three-body problem). A test on
the Li–Liao catalog, and a check of whether the per-pair dial sign can be read from it, would be the follow-up.

## 1. The return map

![Figure 1](figures/fig1-section-and-map.png)

*Figure 1. Left: the figure-eight; all three bodies run along the one grey curve. The start (Li–Liao's form: bodies at
(−1, 0), (1, 0), (0, 0)) is a line-up with body 3 in the middle. Centre: the four coordinates. Right: 1,500 returns from
three starts near the eight, plotted in two of the four coordinates; they stay on rings, slices of invariant tori.*

**Section.** States with centre of mass 0, momentum 0, angular momentum L = 0 and energy E, at the moments the
triangle's oriented area passes through zero (shrinking), with body 3 between bodies 1 and 2. Turn the frame so body 2
is to the right of body 1. Then four numbers fix a point: p = |r₂ − r₁|, q = where body 3 sits from the midpoint,
u₁ and u₂ = the speeds of bodies 1 and 2 along the line. The sideways speeds are a multiple of
(x₃ − x₂, x₁ − x₃, x₂ − x₁) (from zero momentum and L = 0), sized by the energy and signed by the crossing direction
(`section.to_state`). F sends a point to the next such line-up. Every Li–Liao orbit starts on this section.

Count: 12 phase-space numbers, minus 4 for the centre of mass, 2 for rotation and L, 2 for time and energy, leaves 4.
F preserves volume (and the symplectic form).

**Two independent checks.** (A1–A4)

- The figure-eight returns to itself after one hit per period, to 2 × 10⁻¹¹ before polishing.
- The Jacobian of F (finite differences of the section map) and the monodromy of the variational equations reduced
  to the same 4 × 4 block (`cb3.reduced_block`: Φ on {dCOM = dP = dH = dL = 0} modulo the flow and rotation
  directions) give the same Floquet traces to 10⁻⁶.
- Figure-eight multipliers: e^(±0.0529 i) and e^(±1.8730 i), all on the unit circle (E-E). The small angle is the pair
  that meets Simó's H orbit at the fold (section 4).

The integrator is Dormand–Prince 5(4) at tolerance 10⁻¹³ with the 12 × 12 variational equations carried along the same
steps (`tbvar.c`); section hits are found by bisection inside the step.

## 2. Why the index can't see stability here

A 4 × 4 symplectic matrix has two pairs of multipliers (λ, 1/λ), each summarised by its trace τ = λ + 1/λ. A pair is
elliptic (E) if |τ| < 2, hyperbolic (H) if τ > 2, inverse hyperbolic (IH) if τ < −2; if the traces are complex the four
multipliers form a quartet off the circle (CU, "complex unstable"). Since det(I − M) = (2 − τ₁)(2 − τ₂), and
τ for Mⁿ is a Chebyshev polynomial of τ,

> **i(Fⁿ) = (−1)^#H · ((−1)^(n+1))^#IH** for a nondegenerate fixed point (check B1 on 108 orbits, n = 1..4).

This is the known formula index = (−1)^(number of real multipliers > 1), read for symplectic spectra. It depends only
on the number of hyperbolic and inverse hyperbolic pairs mod 2. So for every n:

| Class | i(F) | i(F²) | all iterates | stable? |
|---|---|---|---|---|
| E-E | +1 | +1 | +1 | linearly stable |
| H-H | +1 | +1 | +1 | no |
| IH-IH | +1 | +1 | +1 | no |
| CU | +1 | +1 | +1 | no |
| E-H | −1 | −1 | −1 | no |
| E-IH | +1 | −1 | alternates | no |
| H-IH | −1 | +1 | alternates | no |

In the plane the same formula has one pair, and the index pair (i(F), i(F²)) = (sign R, sign R(1 − R)) is exactly the
residue class, which is what entropy-topology found. One more dimension and three unstable classes join the stable one.

![Figure 3](figures/fig3-catalog-tally.png)

*Figure 3. Left: Floquet classes of the 692 Li–Liao orbits that integrate (grey) and of the 108 whose type is safe from
rounding (coloured). Right: the 52 safe orbits with index +1 for every iterate; 4 are stable.*

**The catalog.** Li and Liao's 695 orbits (equal masses, zero angular momentum, collinear start) give their start to
10 digits. Floquet pairs near τ = ±2 are sensitive to that rounding: refining four orbits as fixed points of F changed
their class (I.A14 went from E-E to E-H). So the tally keeps only orbits whose traces sit at least 0.01 from ±2 (and
from each other): 108 orbits, 4 E-E, 21 E-H, 2 E-IH, 26 H-H, 33 H-IH, 7 IH-IH, 15 CU. Of the 52 with index +1 at every
iterate, 4 are stable (B6). For comparison, Dmitrašinović, Hudomal, Shibayama and Sugita (2018) found 21 linearly
stable orbits among the more than 200 they studied (another catalog).

## 3. The index from map values only

![Figure 4](figures/fig4-value-only-index.png)

*Figure 4. Left: for 15 orbits refined to a residual below 10⁻⁸, the index of F^k from the Jacobian (orange squares)
and counted from map values only (blue dots). One count (II.C9) did not settle. Right: the method, drawn in 2D.*

The index of an isolated fixed point is the degree of g = F^k − id on a small sphere round it. In four dimensions:
cut the boundary of a small box into 3-simplices (Kuhn's triangulation of each of its eight 3-cube faces), evaluate g
at the corners only, and count the simplices whose image cone contains a fixed direction, each with the sign of its
image's orientation times its own (`degree.py`). Linear maps give the sign of the determinant exactly; (z₁², z₂) gives 2
(B3). On the three-body map the count is repeated for two box sizes, two subdivisions and five directions, and accepted
only when all agree.

- 14 of 15 well-refined orbits settle, and every settled count equals sign det(I − DF^k) (B4). The 15th (II.C9, an
  H-H orbit with 12 hits per period and residual 3 × 10⁻⁹) never settled.
- Live in the checks: figure-eight +1, Simó's H orbit −1 (B5).
- The box is shaped by (DF^k − I)⁻¹ so the count is well conditioned. That uses the Jacobian to decide where to sample;
  the count itself uses only values of F. A plain cube gives the same answers where it settles.
- Counts for F² (two periods) settled less often: the orbits are unstable enough that a box small enough to be linear
  is close to the refinement error. Those are not used for any claim.

## 4. Twins

![Figure 5](figures/fig5-twins.png)

*Figure 5. Left: the 7 same-word pairs found earlier ([three-body-word-formula](../three-body-word-formula/README.md)),
each orbit placed at its index, coloured by class, dot size growing with instability. Right: the figure-eight and its
twin on either side of the wall τ = 2 in the pair that folds.*

- **The figure-eight and Simó's H orbit** (C1–C2): one pair has τ = 1.9972 (elliptic) for the eight and 2.0084
  (hyperbolic) for the twin; the other pair is −0.595 against −0.589, elliptic in both. One pair straddling τ = 2 with
  the other unchanged is what a fold looks like (Galán and coauthors 2002 found this fold as a mass changes). Indices
  +1 and −1: here +1 is the stable one.
- **All 7 pairs** (C3–C4): 5 have opposite indices, as fold partners must. II.C66 and II.C75 do not; that fits an
  eigenvalue-1 bifurcation in one of them since their birth, or a birth that is not a simple fold (not traced here).
  In the 5 opposite pairs, the +1 member is the less unstable in 3 (I.A1, II.C73, II.C192) and the more unstable in 2
  (II.B8, II.C246). Away from the figure-eight both twins are unstable, and +1 often marks the fully hyperbolic one.

What it contributes: the index separates twins wherever they were born at a fold and has the right sign rule, which is
why the pair can appear and disappear together. It does not say which twin is stable unless one of them is E-E.

## 5. Isolation: where a Conley index exists

![Figure 6](figures/fig6-isolation.png)

*Figure 6. 24 random starts at each distance from four fixed points, followed for 300 returns. A start "stays" if it
never gets farther than 50 times its starting distance.*

The Conley index belongs to an isolated invariant set: one that is the largest invariant set inside some box and does
not touch the box's edge.

- **Fully hyperbolic points** (H-H, IH-IH, H-IH, CU) are isolated; their Conley index is the pointed 2-sphere (the
  unstable directions) with the map acting by +1 or −1 (the sign is the orientation of the unstable plane; +1 for H-H,
  IH-IH and CU, −1 for H-IH). Every start near I.A4 (H-H) and II.C24 (CU) left, at all four distances (D2).
- **The figure-eight** keeps every start within 10⁻³ (D1), on rings like those in figure 1. Kapela and Simó proved KAM
  tori around it. Tori that fit inside a box but reach its edge make every small box fail to isolate, so the eight has
  no pointwise Conley index. The same goes for any elliptic pair that satisfies the usual twist condition (invariant
  circles in its centre plane), so E-H points are not isolated either, though their staying set has zero volume and
  the test in figure 6 can't see it.
- Exceptions exist: Fayad and Saprykina (2016) built smooth Hamiltonians with isolated elliptic fixed points. The
  statement above is for the generic, twist case.

So in three bodies the Conley index proper certifies full instability. Its shadow, the index, separates twins.
Neither certifies stability.

## 6. The dial reading of one Floquet pair

![Figure 7](figures/fig7-dial-reading.png)

*Figure 7. Left: squaring the map sends a pair's Greene residue R = (2 − τ)/4 to 4R(1 − R) (blue dots: pairs of the
refined orbits). Centre: an elliptic pair, written as x = tan(θ/2), iterates by the tangent law; a hyperbolic pair,
x = tanh(μ/2), by velocity addition. Right: the index is the product of the pairs' dial signs.*

- **E1.** R(F²) = 4R(1 − R) for each pair (Chebyshev doubling), on every refined orbit.
- **E2.** For an elliptic pair with angle θ, x = tan(θ/2) and the n-th iterate is x ∘ x ∘ … ∘ x in the law
  (x + y)/(1 − xy), our two-dial law at k = 0, d = +1. For a hyperbolic pair with log-multiplier μ, x = tanh(μ/2)
  composes by (x + y)/(1 + xy), d = −1 (velocity addition, rapidities adding). This is the jet: the multipliers are the
  derivative of the return map, as Brennan's frame-dragging note asks for in place of plain velocity addition.
- **E3.** The index factor of a pair is sign(2 − τ) = sign(4R), and R = x²/(1 + x²) > 0 (elliptic) or
  −x²/(1 − x²) < 0 (hyperbolic): the factor is the sign of the pair's curvature dial d.

What the dials add: one language for the whole failure. The stable/unstable wall of a pair is the Δ wall of the dial
plane (d changes sign), the same wall that switched entropy on and off last round. The index multiplies the pairs'
d-signs and so loses them in four dimensions. No new dial: k = 0 throughout and only the sign of d matters.

## 7. Neighbours

- **Dial map** ([dial-map](../dial-map/README.md)): only the curvature dial d is used, at ±1. New here: a Floquet pair
  as a point on the line k = 0, and the index as the product of d-signs.
- **Mathlib** ([mathlib-map](../mathlib-map/README.md)): nothing to port yet. The parity law is linear algebra
  (det(I − M) = Π(2 − τᵢ) for a block-diagonal symplectic M) and could be a small Lean lemma; I don't know of a Brouwer
  degree in Mathlib, so the value-only count has no Lean home.
- **three-body-order:** item 7 above; its stability numbers should be replaced by `floquet.json` here.

## Assumptions

- **Equal masses, G = 1, zero angular momentum**, planar. Energy fixed per orbit at its Li–Liao start.
- **Floquet types** come from double-precision variational equations (tolerance 10⁻¹³). Pairs within 0.01 of τ = ±2 are
  left out of the tally because catalog rounding can move them across; refinement of 21 orbits changed 4 classes.
- **The parity law** holds for nondegenerate fixed points (no multiplier equal to 1). It is exact; the catalog checks
  only confirm it was applied right.
- **Value-only counts** are numerical: an outer box, a fine enough triangulation and the absence of other fixed points
  in the box are checked by agreement across sizes and directions, not proved. Interval arithmetic would make them
  proofs; that was not done.
- **Isolation** (section 5) uses KAM for elliptic pairs with a twist condition; the twist was not checked for each orbit,
  only for the figure-eight through Kapela and Simó's proof. The staying test is a picture, not a proof.
- **Twins:** pairing by identical word, from the earlier page; why II.C66 and II.C75 break the sign rule was not traced.
- **"Stable" means linearly stable** (E-E) throughout. Long-time stability needs KAM on top.

## Known vs new

| Item | Status | Source |
|---|---|---|
| Figure-eight: found, existence proved | Known | Moore 1993; Chenciner and Montgomery 2000 |
| Figure-eight linearly stable; KAM stable | Known | Simó 2002 (numerical); Roberts 2007 (ETDS, computer-assisted linear stability); Kapela and Simó (computer-assisted, linear and KAM) |
| Figure-eight's fold with Simó's H orbit | Known | Galán, Muñoz-Almaraz, Freire, Doedel and Vanderbauwhede 2002 (PRL 88, 241101) |
| 695 orbits, equal masses, L = 0 | Known | Li and Liao 2017 |
| Linear stability counts of zero-L orbits | Known (other catalog) | Dmitrašinović, Hudomal, Shibayama and Sugita 2018 |
| Index of a nondegenerate fixed point = (−1)^(# real multipliers > 1); index of iterates | Known | Lefschetz, Dold 1983; standard |
| The parity law for symplectic spectra and the table of section 2 | Elementary consequence; our statement | this page |
| Index pair = residue class in the plane | Elementary | entropy-topology |
| Conley index; isolating blocks; Conley index of a hyperbolic point = pointed sphere | Known | Conley 1978; Conley and Easton 1971; Easton 1975 |
| Isolating blocks in the three-body problem | Known | Easton 1971; Moeckel 2004 (collinear relative equilibria) |
| Generic elliptic points accumulated by KAM tori; smooth isolated elliptic exceptions | Known | KAM; Fayad and Saprykina 2016 |
| Conley–Zehnder index and its iteration | Known | Conley and Zehnder 1984; Bott 1956; Long 2002 |
| Rigorous Conley–Zehnder indices in the restricted three-body problem | Known | Joung and van Koert 2024 |
| Reduced 4D section map of the planar problem at L = 0, value-only 4D degree counts on it | Our computation | this page |
| Catalog tally: 4 E-E among 52 orbits with index +1 at every iterate | Our computation | checks A5, B6 |
| Floquet pair as a dial point (tangent law / velocity addition); index = product of d-signs | Our reading of known pieces (Chebyshev, half-angle addition) | this page |

**Searches (Nov 27):** "fixed point index iterates symplectic map elliptic-elliptic hyperbolic-hyperbolic"; "figure-eight
three-body linear stability Floquet multipliers Roberts Kapela Simó"; "Conley index isolating block three-body problem
Moeckel Easton"; "Galán Muñoz-Almaraz Freire Doedel Vanderbauwhede figure-8 fold"; "Conley index Hamiltonian isolated
invariant set elliptic KAM". None stated the parity table as a remark about stability; it is too elementary to call new.
The figure-eight's multiplier angles (0.0529, 1.8730) were not compared with a published value: the papers found give
the method, not the numbers, in what I could read.

## Files

- `tbvar.c` → `libtbvar.so` (built on first import): the equations of motion with the 12 × 12 variational equations,
  and the section-hit finder.
- `cb3.py`: integration, conserved quantities, the reduced 4 × 4 block, Floquet traces, classes, index from traces.
- `section.py`: the four-dimensional return map F, its Jacobian, Newton refinement, hits per period.
- `degree.py`: the value-only degree count in four dimensions.
- `floquet_catalog.py` → `floquet.json`: Floquet classes of the 695 orbits.
- `curated.py` → `curated.json`: 21 orbits refined, typed and counted.
- `twins.py` → `twins.json`: the 7 same-word pairs.
- `analysis.py` → `analysis.json`: the tally and the isolation test.
- `checks.py` → `checks_output.txt`: 21 checks (A the map, B the index, C twins, D isolation, E dials).
- `figures.py` → `figures/fig1`–`fig7`.
