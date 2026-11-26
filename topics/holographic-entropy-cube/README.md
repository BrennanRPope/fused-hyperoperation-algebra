# DRAFT · Free exploration: spacetime is the entanglement pattern that passes a hypercube test

Nov 25, 2026. Open exploration, asked for by Brennan ("come back with something you find interesting... the goal is
finding this North Star"). Checks: `checks.py` → `checks_output.txt` (23/23 pass). Pictures: `figures.py` →
`figures/`. Second pass: `checks2.py` → `checks2_output.txt`, pictures from `figures2.py`. Pushed Nov 25 as a draft (working folder `exploration/`). Marked DRAFT: the computations reproduce known results; the "ours" items below are
readings and proposals, not theorems yet.

## The finding in two spoken sentences

When a quantum system has a geometric spacetime behind it, the entanglement of its pieces sits on the corners of a
hypercube, and our corner sum at k = −1 is exactly the quantity that tells "has a geometry" apart from "doesn't".
And as of 2025 there is a proof that an entanglement rule holds for every geometry exactly when one hypercube can be
mapped onto another without stretching a single edge, so in this corner of quantum gravity, proofs are maps between
cubes.

## Why I picked it

Three things made me stop.

1. **It is the North Star's best-tested bridge.** Ryu–Takayanagi (2006) says that in a quantum system with a gravity
   dual, the entanglement entropy of a region equals the area of the smallest surface in the curved spacetime that
   hangs off that region, divided by 4G. Quantum information on one side, geometry on the other, one equation.
   From it, people derived Einstein's equations themselves (linearized) from the first law of entanglement
   (Faulkner, Guica, Hartman, Myers, Van Raamsdonk 2013), and in 2019 its quantum version (the "island" formula)
   solved the black-hole information paradox's Page curve. If quantum mechanics and general relativity meet anywhere
   we can compute, it is here.
2. **Its natural coordinates are our cube.** With n regions, the entropies S_T of all the unions T are a function on
   the corners of the n-cube. The quantities that physicists actually use (mutual information, tripartite
   information I₃) are our n-fold cube sum with corner weight k^(rank−1), at k = −1.
3. **The proofs are pictures of hypercubes.** Which entanglement patterns can come from a geometry is an open
   problem (the "holographic entropy cone", fully known only up to five regions). The standard proof method maps one
   Boolean cube into another while never increasing Hamming distance. In 2025 Bao, Furuya and Naskar argued that
   this cube-map test is not just sufficient but necessary. That is about as close to "Hypercube Composition Theory
   applied to quantum gravity" as the literature gets, and nobody has said it in those words.

## The pictures

### Fig 1. Entropies on the cube; I₃ is our corner sum at k = −1

![fig1](figures/fig1_entropy_cube.png)

Each corner of the 3-cube is a union of regions; the empty corner has S = 0. Weight each corner by k^(rank−1) and
add: at k = −1 you get the tripartite information I₃ = S_A+S_B+S_C − S_AB−S_AC−S_BC + S_ABC (check A2). This is the
same sum that, fed product weights, gives a ∘ₖ b ∘ₖ c (check A1, all k).

### Fig 2. Where the entropies come from: shortest curves in the hyperbolic disk

![fig2](figures/fig2_rt_disk.png)

A time slice of three-dimensional anti-de Sitter space is the hyperbolic disk; its boundary circle is where the
quantum system lives. The entropy of a boundary arc is the length of the shortest curve (a geodesic, a circle meeting
the edge at right angles) that ends on the arc's two ends. For two arcs, two pairings compete: (b) far apart, each arc
keeps its own curve and they share no information; (c) close together, the curves reconnect, the shaded region (the
entanglement wedge) becomes one piece, and mutual information turns on. That jump is a geometric phase transition.

### Fig 3. Proof by hypercube

![fig3](figures/fig3_contraction_map.png)

The monogamy of mutual information, MMI: S_AB + S_BC + S_AC ≥ S_A + S_B + S_C + S_ABC, equivalently **I₃ ≤ 0**
(Hayden, Headrick, Maloney 2011). It holds for every geometry and fails for many quantum states.

The proof (Bao et al. 2015): give each term on the left a bit and each term on the right a bit. Each region (A, B, C
and the "rest of the universe" O) marks a corner on both cubes. If you can map the left 3-cube into the right 4-cube,
sending each region's corner to its own corner, without ever stretching a distance, the inequality holds for every
minimal-cut geometry. `checks.py` searched all maps: **there is exactly one**, and it is strikingly simple (panel a):
the four region corners go where they must, and the four other corners all collapse to the single point 0001. All 12
edges land on edges.

Panel (b) shows why the reverse inequality is false: one corner of its domain is a single step from A, B, C and O at
once, and no corner of the target cube is within one step of all four. Some edge must stretch, and Bao, Furuya and
Naskar (2025) argue that a forced stretch can always be turned into a geometry that breaks the inequality (see the
corrections below for the exact scope).

The same search finds the map for the five-region cyclic inequality (a 5-cube into a 6-cube) and finds none for its
reverse (check B).

### Fig 4. The test separates geometric states from quantum ones

![fig4](figures/fig4_who_passes.png)

- 1,500 random weighted graphs, entropy = minimal cut (the discrete version of Ryu–Takayanagi): I₃ ≤ 0 every time.
- 3,000 random triples of arcs on the AdS₃ boundary, entropy = geodesic length: I₃ ≤ 0 every time.
- Random 4-qubit states: 0.2% of typical (Haar) states and 13% of sparse states have I₃ > 0. The GHZ state has
  I₃ = +1 bit exactly; the W state +0.245. A perfect tensor (the 4-qutrit absolutely maximally entangled state, the
  building block of holographic codes) passes with I₃ = −2.
- A state with I₃ > 0 cannot be the boundary of any classical geometry. GHZ-type entanglement (everyone correlated
  through one shared bit) is the wrong shape for spacetime; perfect-tensor entanglement is the right one.

### Fig 5. The k = −1 sums are finite where single entropies are not

![fig5](figures/fig5_cutoff.png)

A single region's entropy diverges as you resolve the boundary more finely (it grows like log(1/cutoff)). The k = −1
sums over a face (mutual information) and over the whole cube (I₃) don't move at all (check E: identical to 1e-14 at
two cutoffs a thousand times apart).

### Fig 6. Why k = −1, and no other dial value

![fig6](figures/fig6_why_k_minus_1.png)

The divergent part of each entropy belongs to one region at a time (it is "one-body", a sum over the regions'
endpoints). For a one-body function the cube sum is n(1+k)^(n−1) times its average, which vanishes only at k = −1
(checks A3, A4). So k = −1 is forced: it is the only dial setting at which the cube sum forgets the cutoff and keeps
only the shared, geometric information.

## How it might touch the North Star

The open question "which quantum states have a spacetime?" has a partial, checkable answer: their entropy cube must
pass every hypercube-contraction test. MMI is the first test (three regions), and five regions are fully mapped;
for six regions, Hernández-Cuenca, Hubeny and Jia (2024) found more than 1,800 new inequalities. Each one is a
statement about one Boolean cube mapping into another. A framework whose home is the hypercube, with a Lean library
that already proves finite cube facts, is suited to this work in a concrete way: a contraction map is a finite object,
so each inequality's proof is decidable by computation.

This does not give a quantum theory of gravity. It gives a sharp necessary condition that any such theory's states
must meet when they look classical, phrased in our cube's own terms, inside the one framework (holography) where
quantum states and curved geometry are already known to be the same thing.

## Known vs ours

| Piece | Status | Source |
|---|---|---|
| Entropy = minimal area / 4G | Known | Ryu & Takayanagi 2006 (hep-th/0603001) |
| MMI, I₃ ≤ 0 for holographic states | Known | Hayden, Headrick, Maloney 2011 (arXiv:1107.2940) |
| Holographic entropy cone; contraction-map proof method | Known | Bao, Nezami, Ooguri, Stoica, Sully, Walter 2015 (arXiv:1505.07839) |
| Five-region cone complete | Known | Hernández-Cuenca 2019 (arXiv:1903.09148) |
| Over 1,800 new six-region inequalities | Known | Hernández-Cuenca, Hubeny, Jia, JHEP 08 (2024) 238 (arXiv:2309.06296) |
| Deterministic filling of contraction maps; partial-cube reformulation | Known | Bao, Furuya, Naskar, JHEP 06 (2024) 039 (arXiv:2403.13283); JHEP 03 (2025) 117 (arXiv:2409.17317) |
| Contraction map necessary as well as sufficient | Known (2025 claim, read in full on the second pass; scope and caveats under "Corrections") | Bao, Furuya, Naskar, JHEP 12 (2025) 140 (arXiv:2506.18086) |
| Einstein's equations from entanglement | Known | Faulkner, Guica, Hartman, Myers, Van Raamsdonk 2013 (arXiv:1312.7856) |
| Multipartite information as alternating (Möbius) sums over subsets | Known | McGill 1954; standard |
| I₃ and the other multipartite informations are our cube sum at k = −1 | **Ours (reading)** | this note, checks A1–A2 |
| k = −1 is the only dial value that removes the one-body (cutoff) part; cube sum of a one-body function is n(1+k)^(n−1)·mean | **Ours (elementary, but stated in our terms)**; the cutoff-freedom of "balanced" inequalities is known | checks A3–A4, Fig 6 |
| MMI's contraction map is unique, and collapses all four non-region corners to one point | Computed here; may be noted in the literature, I didn't find it said | check B |
| Lean-verified contraction maps as a route into the six-region cone | **Ours (proposal)** | next steps |

## Next steps (proposals, nothing started)

1. Put MMI and the five-region cyclic inequality into the core-Lean library as decidable statements (a map on
   `Fin 2^m → Fin 2^n`, checked by `decide`). No mathlib needed.
2. Ask whether other dial values mean something here: Rényi entropies in holography are areas of "cosmic branes"
   (Dong 2016), and the Tsallis form already composes by our ∘ₖ. A k-deformed MMI could be a lead or a dead end.
3. Run our proof engine on a handful of the six-region inequalities and compare with the published contraction maps.

## Runners-up

- **Gravity as the square of electromagnetism (Kerr–Schild double copy).** Schwarzschild's metric is flat space plus
  φ·kk with φ = 2M/r, which is the square of our river speed v² = 2M/r, and its "single copy" φk is the Coulomb field;
  it would join our general-relativity and electromagnetism pages.
- **Curvature from counting (causal sets).** The Benincasa–Dowker action gets spacetime curvature by counting small
  causal intervals with fixed alternating integer weights, close kin to our corner weights.

---

# Second pass (Nov 25, evening): following it deeper

Brennan: "Whoa that is a great find!!! Follow it deeper." Checks: `checks2.py` → `checks2_output.txt` (32/32 pass),
helpers in `checks_lib.py`. Pictures: `figures2.py` → figures 7–12. Still DRAFT.

## The result in plain words

The project's own operations show up inside holographic entanglement. In the simplest curved space (AdS₃), when you
glue two neighbouring regions together, their cross-ratios combine by our union rule ∘₋₁; the flip φ₋₁ turns that
into ∘₁, the height-flow law E(a+b) = E(a) ∘₁ E(b); and holographic mutual information is just the log of that
flipped number. For three neighbouring regions, monogamy of mutual information then comes down to one line:
X ∘₁ Y ≥ X·Y, and the slack in the inequality is exactly the "+X+Y" of the circle operation. The formulas are
known; reading them through ∘₋₁, φ₋₁, E and ∘₁ is ours.

## Corrections from reading the 2025 paper in full

Bao, Furuya and Naskar, "On the completeness of contraction map proof method for holographic entropy inequalities"
(arXiv:2506.18086, JHEP 12 (2025) 140):

- The theorem is for **linear rules with whole-number (so also rational) coefficients**. It uses plain Hamming
  distance, with a term that appears twice counted as two bits. (Bao et al. 2015 used a weighted version; our search
  handles both.)
- "Valid" means: holds for every geometry obeying Ryu–Takayanagi. The proof **assumes every such geometry has a graph
  model**, which the authors say they expect but don't prove. Nonlinear rules are not ruled out.
- The necessity step builds a violating graph by giving a stretched edge infinite (or very large) weight. Some steps
  are informal and a couple of equations look mistyped; the idea is clear, the rigour is not complete.
- Their worked examples, a "deformed MMI" and a five-region candidate, have no contraction map in our search either
  (checks F). Our search also finds maps for the published five-region facets "MMI uplift" and "cyclic" (checks F).

So in the first-pass text, "show" is now "argue", and the necessity result should be cited with these caveats.

## Fig 7. The ladder: union rule → flip → height flow

![fig7](figures/fig7_gluing_ladder.png)

Take a boundary region A and two side-by-side regions B₁, B₂. The cross-ratio x(A,B) is the one Möbius-invariant
number of four boundary points. Exact identity (sympy, check H1): **x(A, B₁B₂) = x(A,B₁) ∘₋₁ x(A,B₂)**. That is the
union rule ("at least one"), the project's k = −1. The flip φ₋₁(x) = x/(1−x) (the project's odds) carries ∘₋₁ to ∘₁
(H2). Writing the odds as X = E(t) = eᵗ − 1 with t = −log(1−x) (H4) makes gluing plain addition of t, and
E(t₁+t₂) = E(t₁) ∘₁ E(t₂) is the height-flow law (H3). Ryu–Takayanagi gives I(A:B) = (c/3)·max(0, log X), checked
against direct geodesic computations in 3,000 configurations, in vacuum and at finite temperature (BTZ black hole),
error below 10⁻¹¹ (H5).

## Fig 8. Monogamy is ∘₁ beating multiplication

![fig8](figures/fig8_mmi_from_circle.png)

I₃(A:B₁:B₂) = log X₁⁺ + log X₂⁺ − log(X₁ ∘₁ X₂)⁺ (⁺ means "or 0 if negative"). Since X₁ ∘₁ X₂ = X₁X₂ + X₁ + X₂ ≥ X₁X₂,
and also ≥ each of X₁, X₂, I₃ ≤ 0 in every phase. When all three pairs are connected, I₃ = −(c/3)·log(1 + 1/X₁ + 1/X₂)
exactly (H7, 274 configurations): the "+X+Y" of the circle operation is MMI's slack. This covers one family (three
neighbouring intervals in AdS₃), not MMI in general.

## Fig 9. Turning the dials

![fig9](figures/fig9_dial_sweeps.png)

- **(a) The corner-sum dial k.** On a grid from −2 to 0.5, the only k at which the cube sum is ≤ 0 for every
  min-cut geometry is k = −1. For k ≥ −½ it is ≥ 0 for everything (it follows from subadditivity). No other dial
  value carves out a known cone.
- **(b) Softening the min cut** (log-sum-exp at softness h, the Maslov / tropical dial the operations catalog already
  lists under "rescale"). h = 0 is the min cut. For h > 0 these are Ising free energies, which is what the averaged
  Rényi-2 entropy of a random tensor network with bond dimension D = e^(1/h) becomes (Hayden et al. 2016; I did not
  check our normalization against theirs exactly). **Numerical observation:** with positive weights, no violation of
  subadditivity, MMI, the MMI uplift or the five-region cyclic rule at any h tested (up to h = 100 for MMI and h = 10 for the five-region rules, in side probes beyond checks.py); with
  mixed-sign weights 30% break MMI. I did not find this stated in the literature. Subadditivity here looks like the
  FKG inequality for ferromagnets; I have no proof for MMI or the five-region rules. **This is a lead worth
  following**: if true in general, the holographic cone would hold for a whole family of "fuzzy" geometries, not
  just the sharp classical limit.
- **(c) Rényi index.** Generic random states break MMI at every α from ¼ to ∞; sparse states about 10–16% of the time.
  Changing the entropy doesn't rescue them.

## Fig 10. Which basis makes the rules short?

![fig10](figures/fig10_rule_length.png)

Rewriting each rule in the dial-k corner basis: SA and MMI become a single term only at k = −1. The MMI uplift is
shortest at k = 0 (7 terms vs 9 at k = −1), because it is I₃ of merged parties (A, BC, DE), and the cyclic rule ties
at k = −1 and 0. Every other k spreads all of them out to every corner. So k = −1 is the natural basis for the
"atoms" (Hernández-Cuenca, Hubeny and Jia work in this basis), but not uniformly the shortest.

## Fig 11. The concrete North Star link

![fig11](figures/fig11_first_law.png)

- **First law of entanglement, checked from Ryu–Takayanagi in d = 2, 3, 4, 5, 6** (sympy, check J1): for a ball of
  radius R in a slightly heated state, the change in the minimal surface's area gives δS = m R^d / (8G(d+1)) per unit
  solid angle, which equals the change in modular energy δ⟨H⟩ exactly.
- Faulkner, Guica, Hartman, Myers and Van Raamsdonk (2013) proved the converse: if δS = δ⟨H⟩ holds for every ball,
  the bulk obeys the linearized Einstein equations. So Einstein's equations (to first order) are the statement
  that entanglement obeys a first law.
- At the next order the two sides separate, and the gap is relative entropy, which is never negative (BTZ, exact,
  J2). Lashkari, Lin, Ooguri, Stoica and Van Raamsdonk (2016) read that positivity as gravitational positive energy.
- **Where our operation enters, honestly:** the first law itself is linear, so ∘ₖ adds nothing there. What survives
  at finite temperature is the gluing structure: the ∘₋₁ law for cross-ratios holds in the BTZ black hole too (the
  distances become sinh's, and the identity still holds; check H5, BTZ). That is checkable but not yet a new
  prediction.

## Fig 12. GHZ, the state with no geometry

![fig12](figures/fig12_ghz.png)

Every entropy of an n-party GHZ state is 1 bit, so its entropy function is the all-ones product and its cube sum is
the n-fold circle product 1 ∘ₖ 1 ∘ₖ … ∘ₖ 1 = ((1+k)ⁿ − 1)/k (the q-number [n] with q = 1 + k). Every n passes
through the same point at k = −1, value 1: GHZ breaks the k = −1 test by exactly one bit at every size (check I4).

## Did the cube-map proofs match our structure?

Partly. What matched: the k = −1 corner basis (the multipartite informations), the union rule ∘₋₁ on cross-ratios,
the flip φ₋₁, the height-flow law E, and the tropical "rescale" dial. What did not: the contraction maps themselves
use the cube as a **metric space** (Hamming distance; image graphs are "partial cubes"), which the project has not
used. Our corner weights k^rank are distances from the empty corner, and a contraction map that fixes the purifier
corner never increases rank, but I found nothing deeper than that. The curvature dial (the plane dial) d and the soliton pair weights did
not show up.

## Known vs ours (second pass)

| Piece | Status | Source |
|---|---|---|
| Holographic mutual information in AdS₃ as a function of the cross-ratio, with the phase transition at x = ½ | Known | Headrick 2010 (arXiv:1006.0047) |
| Cross-ratio product identity under gluing | Known (standard cross-ratio algebra) | — |
| Reading it as x₁ ∘₋₁ x₂, flipping to ∘₁, and MI = log E(t) with E the height-flow generator | **Ours (reading)** | checks H1–H5, Fig 7 |
| MMI for neighbouring intervals as X ∘₁ Y ≥ XY, slack −log(1 + 1/X + 1/Y) | **Ours (short proof for this family)**; MMI itself known | check H7, Fig 8 |
| Only k = −1 makes the cube sum sign-definite on geometries | **Ours (numeric, grid)** | check I1 |
| Soft min cut (Ising free energy, ferromagnetic) obeys MMI and the five-region rules at every h tested | **Ours (numeric observation, unproven, literature not found)** | checks I2, I2c |
| Random tensor networks: Rényi-2 ↔ Ising free energies | Known | Hayden, Nezami, Qi, Thomas, Walter, Yang 2016 (arXiv:1601.01694) |
| First law from RT in d dimensions; converse gives linearized Einstein | Known | Blanco, Casini, Hung, Myers 2013; Faulkner et al. 2013 (arXiv:1312.7856) |
| Relative entropy positivity ↔ positive energy | Known | Lashkari, Lin, Ooguri, Stoica, Van Raamsdonk 2016 (arXiv:1605.01075) |
| GHZ cube sum = n-fold ∘ₖ of 1 | **Ours (one line)** | check I4 |

## What looks worth following next

1. **The soft-min observation (strongest lead).** Try to prove MMI for ferromagnetic Ising free energies with boundary
   conditions, or find a counterexample on larger graphs. Check the literature on correlation inequalities (GKS,
   FKG) first. If it holds, it says the holographic cone survives "1/D corrections" in this model.
2. **Lift the ladder beyond three intervals.** For n neighbouring intervals the cross-ratios still glue by ∘₋₁; test
   whether the five-region rules reduce to ∘₁ inequalities in the same way.
3. **Lean.** Both the contraction maps (finite) and the X ∘₁ Y ≥ XY step are small enough for the core-Lean library.

---

# Third pass (Nov 25, night): two proofs (DRAFT)

**The result in plain words.** Both claims are now proved, not just argued. First, softened ("fuzzy") minimal cuts
can never break the monogamy rule, at any softness and on any network with positive weights, so this piece of the
geometry survives being blurred. Second, in the simplest curved space, gluing two neighbouring regions combines
their cross-ratios by our union rule ∘₋₁ because the regions' "shadows" on a minimal surface simply add end to end.

Files: `checks3.py` (24/24 pass, output in `checks3_output.txt`), `figures3.py` (Figs 13–17), `lean/` (three Lean
files, **compiled into the shared library on Nov 26**: see part 3).

## Part 1. The gluing law, proved by shadows (ours)

**Assumptions.** Vacuum AdS₃ (or the planar BTZ black hole, which is the same space in other coordinates), and the
Ryu–Takayanagi formula: the entropy of a boundary region is (c/6) × the length of its minimal geodesic, and for two
intervals one takes the shorter of the two ways of pairing the endpoints. Nothing else.

**Statement.** Let A = (a₁, a₂) and let B = (b₁, b₂) be cut at m into neighbours B₁ = (b₁, m), B₂ = (m, b₂). With
the cross-ratio x(A,B) = (a₂−a₁)(b₂−b₁) / ((b₁−a₁)(b₂−a₂)),

  **x(A, B₁B₂) = x(A,B₁) ∘₋₁ x(A,B₂) = x₁ + x₂ − x₁x₂.**

**Proof (Fig 13).** Drop each endpoint of B onto A's minimal geodesic γ_A along the geodesic that meets γ_A at a right
angle. Call the stretch of γ_A between the two feet the *shadow* of B, with hyperbolic length t(A,B).

1. *Shadow lemma: 1 − x(A,B) = e^(−t(A,B)).* Move a₁ to 0 and a₂ to ∞ with a Möbius map f (this changes no
   lengths and no cross-ratios). γ_A becomes the vertical axis. The right-angle geodesic from a boundary point at
   distance h from the origin is the quarter circle of radius h, landing at height h. On the vertical axis length
   is log(height ratio), so t = log h(b₁)/h(b₂). Ptolemy's relation for four points on a line,
   (a₂−a₁)(b₂−b₁) + (b₁−a₂)(b₂−a₁) = (b₁−a₁)(b₂−a₂), turns 1 − x into exactly h(b₂)/h(b₁).
2. *Gluing.* The shadows of B₁ and B₂ lie on the same geodesic and meet at m's foot, so t(B₁B₂) = t(B₁) + t(B₂).
   Exponentiate: 1 − x₁₂ = (1 − x₁)(1 − x₂). That is the definition of ∘₋₁ read backwards
   (1 − (x ∘₋₁ y) = (1 − x)(1 − y): "complements multiply", like the union of independent events). ∎

**What follows (Fig 14).** Everything is one number in three coordinates. If δ is the distance between the minimal
geodesics of A and B, then x = sech²(δ/2), the odds X = x/(1−x) = csch²(δ/2), and the shadow t = −2 log tanh(δ/2),
with X = eᵗ − 1. So the shadows add (+), the odds glue by ∘₁ (the height-flow law, via the flip φ₋₁) and the
cross-ratios glue by ∘₋₁. With RT, the mutual information is I(A:B) = (c/3)·max(0, log X); A and B share
entanglement exactly when δ < 2 arcsinh 1 ≈ 1.76. Monogamy for A next to B₁, B₂ is then one line:
(e^s − 1)(e^u − 1) ≤ e^(s+u) − 1, the gap being X + Y, so I₃ = −(c/3) log(1 + 1/X + 1/Y) ≤ 0 when all three
connect, and I₃ ≤ 0 in the other phases because log(eᵗ − 1) only grows. The same holds for the BTZ black hole
because conformal factors cancel from every cross-ratio (check L10).

**Checks.** L1–L3 are the algebra (sympy). L4 finds the feet of the perpendiculars by direct minimisation (no
Möbius map) and confirms e^(−t) = 1 − x. L5: shadows add. L6: x = sech²(δ/2). L7: RT mutual information equals
(c/3) max(0, log csch²(δ/2)). L8–L9: monogamy. L10: conformal factors cancel.

## Part 2. Fuzzy min cuts never break monogamy (ours)

**Setting.** Take any network with positive edge weights, boundary regions A, B, C and the rest O. Replace the
minimal cut by a softened one at softness h: S_h(T) = −h log Σ exp(−cut/h), summing over all ways to label the
inside nodes, with T on one side and the rest on the other (minus the same with T empty). As h → 0 this is the
min cut (the holographic entropy of a graph model). For h > 0 it is the free-energy cost of a domain wall in a
ferromagnetic Ising magnet at temperature h, and (in Hayden et al.'s 2016 dictionary, up to normalisation) the
"annealed" second Rényi entropy of a random tensor network whose bond dimension is D = e^(weight/h).

**Theorem (ours as far as we could find).** For every such network and every h > 0,
S_AB + S_BC + S_AC ≥ S_A + S_B + S_C + S_ABC (monogamy of mutual information). Because regions can be any sets of
nodes, the "uplifted" five-party rule (S_ABC + S_ADE + S_BCDE ≥ S_A + S_BC + S_DE + S_ABCDE) follows too.

**Proof (Fig 16).**
1. *Reformulate.* Glue each region into one node; treat A, B, C, O as four Ising spins and sum out everything else.
   The result is the four spins' joint weight, exp(Σ J_U s^U) over even sets U (the magnet has no field, so odd
   terms vanish). Monogamy turns out to be exactly J_ABCO ≤ 0: the four-body coupling the network induces among
   its four terminals is never positive (check K6: I₃ = 8h·J_ABCO). Written out, it is an inequality between 8
   boundary pinnings: Z(++++)Z(++−−)Z(+−+−)Z(+−−+) ≤ Z(+++−)Z(++−+)Z(+−++)Z(−+++).
2. *Four copies.* A product of four Z's is one Z of four independent copies. Stacking the four "even" pinnings
   pins each terminal to a column of the 4×4 Hadamard matrix; stacking the "odd" ones pins them to the odd frame
   v₁ = (+,+,+,−), v₂ = (+,+,−,+), v₃ = (+,−,+,+), v₄ = (−,+,+,+). Negating copy 4 everywhere is a symmetry and
   turns the Hadamard frame into (v₁, v₂, v₃, −v₄) (check K4). So the claim is Z(v₁,v₂,v₃,+v₄) ≥ Z(v₁,v₂,v₃,−v₄).
3. *Rotate.* Write each 4-copy spin in the orthogonal frame: αᵢ = vᵢ·σ/2. The energy σ·σ′ = Σ αᵢαᵢ′ is unchanged,
   terminal j now pushes only coordinate αⱼ, and the 16 spin values become the 8 even sign vectors
   (α₁α₂α₃α₄ = +1) plus the 8 points ±2eᵢ. Every moment of that set is ≥ 0 (checks K2, K3, K5, K9).
4. *Griffiths.* Z(+v₄) − Z(−v₄) = 2Z′·⟨sinh(push on α₄)⟩, a series of odd monomials with positive coefficients.
   With ferromagnetic couplings, positive pushes on α₁, α₂, α₃ and nonnegative single-site moments, Griffiths'
   first inequality (Ginibre's general form) makes each such average ≥ 0. In pictures: A, B and C push α₁, α₂, α₃
   up, and the parity rule α₄ = α₁α₂α₃ makes α₄ lean up too. ∎

**What is known and what is ours.** The four-copy-and-rotate technique is Ellis–Monroe–Newman's (Commun. Math. Phys.
46, 1976), from their proof of the GHS inequality; Griffiths' first inequality is from 1967 (Ginibre 1970 for
general single-site measures). Strong subadditivity of these fuzzy entropies is standard (it is the FKG lattice
property of the terminals' marginal). Applying the technique to the 8-pinning inequality, and so to holographic
monogamy for Ising free energies and annealed tensor-network Rényi-2 entropies at any bond dimension, I did not find
in the literature (searched: random tensor networks at finite bond dimension, holographic entropy cone, annealed
Rényi entropies and domain walls, effective multi-spin couplings of Ising marginals). Treat "ours" as "not found",
not as "certainly new".

**The two-state symmetry is essential.** With 3 or more spin states (Potts magnets, the natural analogue for higher
Rényi entropies of tensor networks) the same fuzzy cut does break monogamy: check K12b gives an explicit network with
positive weights where 3 states give I₃ = +3.3·10⁻⁵ and 2 states give −1.5·10⁻⁴; on a single hub, 8 and 20 states
break it (Fig 15b). Single-hub Ising also reproduces Kadanoff's decimation sign: the generated four-spin coupling
(1/8)[log cosh 4K − 4 log cosh 2K] is ≤ 0 because cosh⁴2K − cosh 4K = (cosh²2K − 1)² (check K11).

**Five regions (evidence only, Fig 17).** Optimising weights to break each of the five-party facets from
Hernández-Cuenca 2019 (items 3–7 of his table; I could not retrieve item 8) found nothing positive (check K13).
Item 3 is proved by the theorem above; items 4–7 remain open for fuzzy cuts.

**What it means physically.** The rule that entanglement can't be shared freely (monogamy) is the signature of a
geometric, Einstein-like bulk. It survives turning the sharp minimal surface into a thermal ensemble of domain walls,
that is, away from the strict large-D (classical-gravity) limit of a random tensor network, as long as the replica
symmetry is the two-state one that Rényi-2 has. It is a small, exact statement about "fuzzy geometry" keeping one
of its defining rules.

**Our dials.** The softness h is the tropical (Maslov) dial: h → 0 sends log-sum-exp to min, which is the rescale
dial already in the project. The proof itself lives at the k = −1 corner: I₃ is the k = −1 corner sum, and in spin
language it is the single top coefficient J_ABCO.

## Part 3. Lean (compiled Nov 26)

`lean/HoloGluing.lean` (core Lean 4.37, the library's style, `grind`, no `sorry`) states with denominators cleared:
Ptolemy's relation, the shadow numerator, the gluing law, "complements multiply" for ∘₋₁, the odds gluing by ∘₁,
the height-flow identity, cancellation of conformal factors, and the monogamy step
max(1,X)·max(1,Y) ≤ max(1, X ∘₁ Y) on Rat. `lean/HoloGluingField.lean` has the same laws with division.
`lean/FuzzyFrame.lean` checks the finite facts of the fuzzy proof by `decide` (orthogonal frame, the copy-4 flip,
nonnegative moments up to exponent 3, where the pinned terminals sit).

**Compiled.** The three files were written on a machine that couldn't download Lean. On Nov 26 they were ported to the 1.0
names and compiled into the shared library as `FusedAlgebra/Topics/HoloGluing`, `HoloGluingField` and `FuzzyFrame`
(19 theorems, no `sorry`, standard axioms only). Two proofs needed help: `HoloGluingField.gluing` splits the divisions
into inverses before `grind`, and `FuzzyFrame.moments_nonneg` uses `decide +kernel`.

**Which steps can't be formalized in this library.** The Ryu–Takayanagi formula (a physical assumption), the
hyperbolic-geometry part of the shadow lemma (core Lean has no hyperbolic plane, logarithms or exponentials), and
Griffiths' first inequality for general graphs (a real but sizeable formalization; Mathlib would make it feasible).

## Known vs ours (third pass)

| Piece | Status | Source |
|---|---|---|
| Ptolemy's relation, cross-ratio invariance, perpendicular feet on a geodesic | Known (classical) | — |
| Holographic MI in AdS₃ = (c/3) max(0, log odds) | Known | Headrick 2010 (arXiv:1006.0047) |
| The shadow lemma 1 − x = e^(−t) and the proof of gluing by adding shadows; x = sech²(δ/2) as the ∘₋₁ rung | **Ours (proof; the facts are elementary hyperbolic geometry)** | checks L1–L10, Figs 13–14 |
| Four-copy rotation, Griffiths/Ginibre inequalities | Known | Ellis–Monroe–Newman 1976; Griffiths 1967; Ginibre 1970 |
| Fuzzy entropies obey SSA | Known in substance (FKG) | — |
| Fuzzy entropies obey MMI for every graph and softness (J_ABCO ≤ 0) | **Ours (proof; not found in the literature)** | checks K1–K10, Figs 15–16 |
| 3-state Potts version breaks MMI | **Ours (explicit counterexample)** | check K12b |
| Five-party facets for fuzzy cuts | **Open (numerical evidence only)** | check K13, Fig 17 |
| Rényi-2 of random tensor networks ↔ Ising free energies | Known | Hayden et al. 2016 (arXiv:1601.01694) |

## Next, if wanted

1. Try the same four-copy idea on the cyclic five-party rule (it needs a 2⁵-copy frame; whether one with the right
   parity structure exists is the question).
2. Write the fuzzy-monogamy theorem as a short note (nothing goes out without your yes).
