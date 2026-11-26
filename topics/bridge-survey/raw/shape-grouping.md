# Shape grouping: shapes and applied folders

Raw notes for the bridge survey. The companion table is `shapes-applied.tsv`: rows S01–S149 cover shapes and A01–A125 cover applied folders.

Sources: the README, or else the page text or checks output, of every shapes and applied folder in both copies. The shared copy was used wherever it is newer, and polar-recast's shared copy adds a literature section. No new literature searches were made, and every verdict is the write-up's own.

Coverage notes:
- **shape-dial-survey** has only `checks_output.txt` (35/38 pass; H1, H2 and J1 fail). It has no write-up, so its rows are provisional and lean on the check labels.
- **simplex/shape_dial_checks_output.txt** finds that Vinberg curvature runs from −1 to 0 with temperature (26/26), but no write-up covers it yet.
- **photosynthesis** has only `photosynthesis_calc.py`, with no write-up. The script splits the Hall & Rao loss chain three ways (waterfall, Shapley and log-proportional). It gets no row, but it feeds item 7 of part (iii).
- **proof-gallery** names no bridges.

---

## The map in one picture

```
                       roundness dial p (polar-coordinates §7, shape-dial-survey)
   cross-polytope  ----------------->  sphere / ball  ----------------->  cube / polydisc
      (p = 1)                            (p = 2)                           (p = inf)
                                           |  square the coordinates,          |  stick-breaking, k = -1
                                           |  forget the phases                |  (simplex §1)
                                           v                                   v
                                        simplex  <------- same object, three charts: v = cos^2(theta) -------
                                           |
     corner weights (polar §6, Fig 7):  simplex -> associahedron -> permutohedron -> cube
     budget region (beyond §1):         cube (k=-1) -> simplex (k=0) -> axes (k -> inf)
     Hilbert metric (simplex §5):       ball (Klein, curved -1)  --temperature-->  simplex (flat Bayes, 0)
     landscape (k, d):                  circle  |  parabola  |  hyperbola      by the sign of k^2 - 4d
```

---

## (i) Where a write-up says two shapes are one object, or one slides into another

| # | Statement | Folder(s) |
|---|---|---|
| 1 | Stick-breaking at k = −1 sends the cube onto the simplex. Cube, simplex and sphere angles are three charts of one thing, with v = cos²θ. | simplex §1–3 |
| 2 | Sphere and simplex are one object: square the coordinates and forget the phases (the Born-rule map). | simplex §3, hypersphere §3, polar-coordinates §3 |
| 3 | The cube is the polydisc, i.e. a product of discs read through radii. | hypersphere |
| 4 | Budget region slides from cube (k = −1) to simplex (k = 0) to the coordinate axes (k → ∞). | beyond-the-cube §1 |
| 5 | In log-radius coordinates every budget region is one simplex of size log(1+k)/k, so the cube is the simplex pushed to infinity. | polar-coordinates §8 |
| 6 | The p dial runs cross-polytope (1), then sphere (2), then cube (∞), with the Maslov/tropical limit at p → ∞. | polar-coordinates §7, shape-dial-survey |
| 7 | Corner weights slide simplex → associahedron → permutohedron → cube (Postnikov generalized permutohedra, y_I = tesseract corner coordinates). | polar-coordinates §6, Fig 7 |
| 8 | The type-B family holds both cube and cross-polytope. | polar-coordinates §7 |
| 9 | Hilbert geometry: the ball is relativity's velocity space (Klein model, curved, precesses) and the simplex is flat Bayes geometry. The two coincide in 1D, where log odds = 2 × rapidity. | simplex §5 |
| 10 | The temperature free-energy body slides from ball to simplex, with curvature running from −1 to 0. | simplex fig 8 / R5; shape_dial checks (no write-up) |
| 11 | Cluster coupling deforms the cube into the associahedron. | number-theory-geometries §3, beyond-the-cube §3 |
| 12 | The permutohedron grows from a point at k = 0. Cube paths are the permutohedron's vertices, and the permutohedron is the shadow of the Birkhoff polytope. | number-theory §2, boolean-lattice, beyond-the-cube §7 |
| 13 | The cut polytope sits inside the elliptope, and the quantum boundary is straight in angle coordinates. | beyond-the-cube §2, polar-recast 6c |
| 14 | The entropy polymatroid is a shape, and the holographic cone is "the net carved". | polar-recast 6b |
| 15 | Pinching takes torus to cylinder to plane. | shape-dependencies, polar-recast §3, shape-directions |
| 16 | The conics split by the sign of k² − 4d: circle / parabola / hyperbola, i.e. Euclid / Galileo / Einstein. | polar-coordinates §2, landscape, simple-laws |
| 17 | Sphere / plane / saddle as d = 1/R², where an imaginary radius gives the saddle. | hypersphere §5 |
| 18 | Spin sphere contracts to the Heisenberg plane (Inönü–Wigner). | polar-recast §4 |
| 19 | Tropical subdivisions with ∘ₖ as the uncut point, and the cube as the q → 1 subspace lattice. | number-theory §1, §5 |
| 20 | Rung sheet: square corner → line → hyperbola → axes. Super log is the change between shapes. | rung-shape, polar-recast 10d |
| 21 | Tesseract dial runs x² = 0, x² = x, e² = ±1 (dual, idempotent, split/complex). | tesseract |
| 22 | Polygons' CHSH value slides from 4 to 2√2 as the polygon rounds. | shape-dial-survey C5 |
| 23 | Orthants are the cube's corners inside the sphere. | hypersphere §7 |
| 24 | The cube splits into n! simplices, one for each ordering. | simplex E2, hypersphere |
| 25 | Fermat curves are the line shrunk by n. | polar-recast 8 |

---

## (ii) Which fields' laws were sorted into each shape family

| Shape family | Laws filed there (folder) |
|---|---|
| **Cube** (and its corners) | probability and logic, inclusion–exclusion, binomial (statistics, boolean-lattice); KdV solitons, ring theory, Rota Möbius functions, path signatures (patterns, boolean-lattice); CISD, Mayer clusters, isotope patterns, spin multiplicities (chemistry 1, 5, 6); Bliss independence, G×E interaction (biology 3, 4); independent risks, claim counts, decrements, annuities (actuarial); rate paths and return cube (econ-forecasting); fan layers, PUE corner (datacenter-energy); committees (consensus); spillover corners, factorial "society tensor" (rulecraft); holographic corner weights (polar-recast) |
| **Simplex** | Bayes, Dirichlet, softmax/Gibbs, Born rule, arcsine law, replicator dynamics, Fisher metric, conservation laws, e, categorical models (simplex, polar-recast) |
| **Cross-polytope** | p = 1 end of the dial, type-B, PR box, e² (polar-coordinates, shape-dial-survey). **No science laws.** |
| **Sphere / ball** | i, π, spin, Coulomb, relativity (celestial sphere, Klein velocity ball, Thomas–Wigner), Kepler/Fock S³, hydrogen n², de Sitter, Born cap, Berry phase, swap rule, d = 1/R², Tsirelson, Gaussian (hypersphere, polar-recast, simplex) |
| **Torus / helix / circle of angles** | roots of unity / Fermat, elliptic curves, Lee–Yang circle, q on the unit circle, rung helix (polar-recast, number-theory, rung-shape) |
| **Conic: circle (turning)** | trigonometry, waves above the barrier (Prüfer), Hückel band, bound states (landscape, simple-laws, chemistry 4) |
| **Conic: parabola (border)** | Galilean addition, Euler characteristic, cusp, Lavoisier at c = 1, second-order kinetics, one-compartment synapse, wave turning point (landscape, chemistry 2, 9, biology 7) |
| **Conic: hyperbola (two walls)** | relativity, logistic, Kalman, epidemics, ee/Horeau, Bayes pooling, telegraph line at (0, −1), Kosambi, Hückel gaps, free particles (simple-laws, chemistry 4, 12, biology 4, statistics) |
| **Permutohedron** | orderings, Plackett–Luce, flag-variety genus prediction, base polytope / core with Shapley at the centre (number-theory, boolean-lattice, beyond-the-cube); actuarial decrement paths implicitly |
| **Associahedron** | φ³ amplitudes, cluster algebras, AdS₃ linking, Y-system / Lee–Yang pentagon, double copy (beyond-the-cube, number-theory) |
| **Elliptope / cut polytope** | quantum correlations only (beyond-the-cube, polar-recast 6c) |
| **Generalized permutohedra / polymatroids** | entropy cone, cores of games, Smets belief functions (polar-coordinates §6, polar-recast 6b) |
| **Others** | half-line for rungs (rung-shape); strain triangle for Navier–Stokes (shape-directions); **log axis**, where polar-recast's inventory files almost all applied sciences (chemistry, biology, actuarial, econ) |

Two lopsided spots stand out in this table:
- Most applied laws are filed only by their composition law (cube corners or log axis), not by a shape.
- Three families hold one field each, or none: cross-polytope, elliptope and associahedron.

Part (iii) is mostly about those gaps.

---

## (iii) Coincidences that no write-up names

Every item below is **unnamed, needs checking**. Each one gives the concrete identity that would have to hold. Items 1–5 are exact algebra I checked by hand, item 5 also numerically for n = 3, 4, 5. For those, what needs checking is whether the cross-link is new and useful, not whether the algebra holds.

1. **One law at (−2, 1) shared by kinetics, neurons and copulas.** *Unnamed, needs checking.*
   - Three write-ups use the same law, x₁₂ = (x₁ + x₂ − 2x₁x₂)/(1 − x₁x₂):
     - second-order kinetics (chemistry 9);
     - one-compartment synaptic summation (biology 7);
     - the AMH/Hamacher θ = 1 end (statistics).
   - In each case the odds x/(1−x) add: k[A]₀t in kinetics, conductance in units of the leak in the neuron.
   - x → −x carries it onto the Euler-characteristic law (a + b + 2ab)/(1 − ab) at (2, 1) (landscape, polar-recast 10a, shape-dial-survey F4).
   - Lavoisier's c = 1 sits at (2, 1) too, via the law of cosines with k = 2c.
   - None of these folders cites another.

2. **Genetics map-function line = the copula dependence line.** *Unnamed, needs checking.*
   - Map biology 4's line 2k + q = −4 onto statistics' AMH/Hamacher line k + q = −1 by rescaling x = 2r, i.e. (k, q) → (k/2, q/4).
   - The endpoints match:
     - Haldane (−2, 0) → independence, AMH θ = 0, (−1, 0);
     - Kosambi (0, −4) → velocity addition, θ = −1, (0, −1);
     - θ = +1 → (−4, 4) in r units, which is the item-1 law in x = 2r units.
   - So map functions are a one-parameter copula family on 2r.

3. **Haldane's addition rule is the piling-up lemma.** *Unnamed, needs checking.*
   - The identity: 1 − 2(r₁ + r₂ − 2r₁r₂) = (1 − 2r₁)(1 − 2r₂).
   - That is statistics' XOR/piling-up at k = −2 applied to recombination fractions (biology 4).
   - It holds trivially. What is new is only the cross-link.

4. **Hückel chain = Lucas corner weight on the q = 1 line.** *Unnamed, needs checking.*
   - The identity: the N-atom Hückel chain polynomial (α = 0, β = 1) equals Lucas U_{N+1}(P = x, Q = 1), the boolean-lattice rank-(N+1) corner weight at (k, d) = (x, 1).
   - boolean-lattice calls its formula "ours, probably folklore". Hückel/Chebyshev (chemistry 4, 12) would be its literature anchor.
   - The levels 2cos(jπ/(N+1)) are:
     - patterns' Beraha points, where the (x, 1) step has period N+1;
     - simple-laws' tilted-axes rotation (2cosθ, 1) on the law-of-cosines line.
   - Rings with all-integer levels exist only for N ∈ {1, 2, 3, 4, 6}, the crystallographic restriction.
   - Frost-circle ring levels are roots of unity, which moves Hückel rings into the torus-of-angles family rather than the log axis.

5. **Cooperative binding zeros are a Lee–Yang circle.** *Unnamed, needs checking.*
   - Two sites (chemistry 8): zeros of 1 + 2x + ωx² sit on |x| = 1/√ω exactly when ω > 1. With y = √ω·x the polynomial is 1 + (2/√ω)y + y², a ferromagnetic two-spin gas.
   - n sites with pairwise ω (Ising-type, Σ C(n,j) ω^{C(j,2)} x^j): zeros lie on |x| = ω^{−(n−1)/2} (checked numerically for n = 3, 4, 5).
   - The MWC fold, (1 + kx)ⁿ = 1 − k, instead puts its zeros on a circle centred at −1/k.
   - The zero locus therefore tells pairwise cooperativity apart from concerted cooperativity, the same distinction chemistry 11 draws by fitting. It links to biology 1's "Lee–Yang dial gas" helix–coil row and to shape-dial-survey's failed Lee–Yang check J1.

6. **Horizon n is the roundness dial p in discounting.** *Unnamed, needs checking.*
   - Econ-forecasting's certainty-equivalent factor is (Σ wᵢ dᵢⁿ)^{1/n}, with dᵢ = 1/(1+rᵢ) per world. That is a weighted n-norm, so the horizon n plays the role of p.
   - "Lowest corner wins" (Weitzman) is then the p → ∞ Maslov/max limit of polar-coordinates §7 and shape-dial-survey M1/M2.
   - Numerically, with rates 1/4/7%: the equivalent rate at n = 1, 10, 100, 1000 is 3.9%, 3.7%, 2.1%, 1.1%, tending to 1%.

7. **Decrement splits are base-polytope vertices and centroid.** *Unnamed, needs checking.*
   - Actuarial's UDD = Shapley should be the centroid of the vertices of the k = −1 base polytope / permutohedron (number-theory B5).
   - The forward waterfall in order σ gives each cause the stick-breaking piece q′ᵢ ∏_{j before i}(1 − q′ⱼ), which is a vertex (simplex §1).
   - Constant force is the log split, i.e. datacenter-energy's "nepers".
   - photosynthesis_calc.py already computes all three splits for one chain without naming the geometry.

8. **Spillover corners are Postnikov y_I.** *Unnamed, needs checking.*
   - Rulecraft's Harsanyi dividends d_I are polar-coordinates §6's y_I.
   - Identity to check: for a supermodular game (all d_I ≥ 0), core = Σ d_I Δ_I.
   - Benedict's synergy A > 0 would be number-theory's k > 0 supermodular / convex game.

9. **Econ's effective dimension lives on the elliptope.** *Unnamed, needs checking.*
   - ED = n²/‖C‖_F² = n²/(n + Σ_{i≠j} C_ij²) is a function on the elliptope (all correlation matrices).
   - Its two ends:
     - ED = 1 exactly at the rank-one ±1 points, the cut-polytope vertices (cube corners);
     - ED = n at the identity.
   - That would give the elliptope a second, non-quantum field. Statistics never uses it.

10. **Synthesis routes are associahedron vertices.** *Unnamed, needs checking.*
    - Convergent vs linear routes in chemistry 5 and 12 are bracketings of a product, i.e. associahedron vertices.
    - Overall yield y^{max depth} is depth-sensitive, so routes would form an associahedron with a non-uniform weight.
    - Today the associahedron hosts only amplitudes and cluster algebras.

11. **Selection is a Hilbert isometry of the simplex.** *Unnamed, needs checking.*
    - Biology 8's discrete selection on odds is simplex's k = +1 Bayes move (Aitchison perturbation).
    - By simplex H4, selection steps should be Hilbert-metric isometries that commute, with no precession, unlike relativity's ball.

12. **Chemistry's n² shells belong in the hypersphere family.** *Unnamed, needs checking.*
    - Chemistry 3 (n², Fock SO(4)) and chemistry 2 (hydrogen) cite Fock 1935, as hypersphere's Kepler/S³ section does.
    - polar-recast's inventory files chemistry under the log axis, and the folders never cross-link.

13. **Two rewrites of Dirac hydrogen.** *Unnamed, needs checking.*
    - Chemistry 2's binding angle (sin θ = Zα, E = mc² cos θ) and hypersphere/orbits' velocity disc of radius αc/κ both rewrite the Dirac ground state.
    - One is a turning angle and the other an area on the velocity sphere.
    - To check: do they agree for excited states (κ ≠ 1)?

14. **One q runs through Rényi, Hill numbers, reaction order and Clayton.** *Unnamed, needs checking. Weak: a shared parameterization.*
    - shape-dial-survey K4 has Rényi q = p-norm p, with k = 1 − q.
    - The same index appears as ecology's Hill-number order q (biology 10, D_q = ‖p‖_q^{q/(1−q)}), chemistry 9's reaction order n = q, and statistics' Clayton generator.
    - All would then index the roundness dial: p = 1 Shannon / first order; p = 2 sphere / Simpson / second order; p = ∞ cube / Berger–Parker.

15. **Duckweed mixtures sit where no binding scheme can go.** *Unnamed, needs checking.*
    - Chemistry 13 (Yonetani–Theorell) says two-inhibitor binding can't pass Loewe (k ≥ 0).
    - Biology 3's measured mixtures sit at α ≈ −0.94 and −0.98, near HSA (k < 0).
    - If α_CT = 1/α_YT on the same variables, those mixtures lie in a region no two-inhibitor binding scheme reaches. That would point to a mechanism other than binding, or to a mismatch of variables.

16. *(Minor)* **Hill's muscle hyperbola is Boyle's hyperbola.** *Unnamed, needs checking.*
    - Biology 6's force–velocity law (F + a)(v + b) = const is simple-laws' Boyle hyperbola, the mult-add law's unit circle.
