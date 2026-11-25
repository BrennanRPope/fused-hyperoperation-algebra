# Fused Hyperoperation Algebra

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24539305.svg)](https://doi.org/10.5281/zenodo.24539305)

Cite all versions with DOI [10.5281/zenodo.24539305](https://doi.org/10.5281/zenodo.24539305); Zenodo also gives each release, such as v0.1.0, its own version DOI.

Brennan Pope's math theory, with the proofs, checks, figures and write-ups gathered in one place. (It was called Fused Mult-Add Algebra until Nov 25, 2026.)

## The idea in plain words

Ordinary arithmetic throws information away. 0 × 4 = 0 forgets the 4, and 4 / 0 is simply "undefined". This project looks for operations that keep that information, and for one simple geometric picture that ties a lot of existing math together.

**Mult-add.** The core operation is

    a ⊕ b = a + b + ab = (1 + a)(1 + b) − 1

Draw a unit square, stretch it by a one way and by b the other. The new area minus the old one is a + b + ab: the two edge strips (the additions, linear) plus the corner patch (the multiplication, area). Nothing is lost. Do it with three numbers and you get a cube whose pieces come in layers of 1, 3, 3, 1, which is a row of Pascal's triangle. That "mult-add cube" is the main foothold.

**Rungs, fused.** The shift x ↦ 1 + x turns mult-add into ordinary multiplication and back. Conjugating any rung of the hyperoperation ladder (adding, multiplying, powers, and so on) by that shift gives its fused version, so every rung gets a mult-add twin, with inverses that map back to the standard ring.

**Numbers as structures.** In this view a number carries dimensions: additions are linear pieces, multiplications are areas, higher orders are higher-dimensional corners. Curves come in through i, which turns the picture into a helix.

**The big goal.** An "upgrade" of the standard ring where calculus falls out, and a common language in which different departments of math and physics can be seen through the same lens. The aim is not to solve famous open problems, but to see how each one looks from this angle.

**The rule.** An outcome only counts when it is geometric, with a geometric underpinning we are sure of before calling it good. Pictures are how the proofs are understood.

**Honesty note.** Each piece found so far matches known mathematics (formal group laws, jet and dual numbers, Clifford algebra, K-theory, the Hirota form of KdV solitons, and more). The possibly new part is the single geometric picture that unifies them. The write-ups say which is which.

## Where it unifies

It's easy to get stuck in one corner of mathematics and miss the forest for the trees. The point of the mult-add box is that one shape shows up in areas that are usually kept apart. Stretch a unit square by a one way and b the other: the new area is (1 + a)(1 + b), and the part beyond the old square is a + b + ab. With n numbers the box has one corner piece for every subset of them, and a dial k can weigh the corner: a + b + k·ab.

![The mult-add box in the middle, with counting, primes, probability, special relativity, general relativity and topology hanging off it](docs/where-it-unifies.svg)

The table adds two areas that came in after the picture was drawn: electromagnetism and Clifford algebra.

| Area | What the box says there | Already known as | Folder |
|---|---|---|---|
| Counting and power sets | Multiplying out ∏(1 + tₓ) gives one term per subset; all dials at 1 give the 2ⁿ subsets, and two random subsets joined have odds t ✳ s. | The multilinear expansion; p-biased product measures | [`fractional-powerset-dials`](topics/fractional-powerset-dials/) |
| Primes and zeta functions | Put dial p⁻ˢ on each prime: subsets of primes are the squarefree numbers, so ∏(1 + p⁻ˢ) = ζ(s)/ζ(2s). Turning the dial k gives ∏(1 − k·p⁻ˢ)^(−1/k), which is ζ(s) itself at k = 1. | Euler products; the fermionic and bosonic primon gases (Julia 1990, Spector 1990) | [`dial-survey`](topics/dial-survey/) |
| Probability and measure | k = −1 is "A or B" for independent events, p + q − pq. Sugeno's λ-measure is the dial with k = λ, and α-determinantal point processes are the dial with k = α (fermions −1, Poisson 0, bosons 1). | Sugeno 1974; Shirai and Takahashi 2003 | [`measure-probability`](topics/measure-probability/), [`dial-point-processes`](topics/dial-point-processes/) |
| Electromagnetism | The wave equation keeps its form under a dial boost only at q = −1/c², so light sets the dial (c²ε₀μ₀ = 1). A boost scales the light-line parts of E and B by D and 1/D, and D² − 1 composes by mult-add. The phase around a loop is charge times flux, the same "loop turn = curvature × area" rule gravity follows. No new prediction. | Maxwell; Lorentz invariance of light; Aharonov–Bohm; Hestenes 1966 (∇F = J) | [`electromagnetism`](topics/electromagnetism/) ([page](https://claude.ai/artifact/A71iTBsgUHZSkcbNUrEacY)) |
| Special relativity | Velocity addition is mult-add in d = 2u/(1 − u). Then 1 + d is the Doppler factor squared, and the rapidity ½·log(1 + d) simply adds. | Rapidities add; Doppler factors multiply | [`relativity`](topics/relativity/) |
| General relativity | In the "river" picture of a black hole, one boost stretches clocks and radial rulers by mult-add inverses: (−g_tt − 1) ✳ (g_rr − 1) = 0. That's why time and space bend equally and light bends twice Newton's angle. | Painlevé 1921, Gullstrand 1922; Hamilton and Lisle's river model, 2008 | [`general-relativity`](topics/general-relativity/) |
| Composition laws | Every law (a + b + k·ab)/(1 − q·ab) is a point on one plane, and k itself is a choice of units (only its sign survives rescaling). Tsallis' q-algebra is the dial with k = 1 − q. | One-parameter Möbius groups (Lie, Klein); Tsallis 1988 | [`landscape`](topics/landscape/), [`dial-survey`](topics/dial-survey/) |
| Clifford algebra | Mult-add's dial is the light-cone (null) element of the split-complex numbers. With vectors as inputs, x ✳ y = x + y + x·y + x∧y stays associative, and x ✳ y − y ✳ x is exactly twice the wedge product. Of 16 pieces of Clifford's work, 12 have exact analogs, 2 partial, 1 only by transfer, and 1 none (his 1878 theorem on special divisors). | The geometric product (Clifford, Grassmann) | [`clifford`](topics/clifford/) |
| Topology | The Todd genus composes by mult-add, the signature by velocity addition, and all of Hirzebruch's χ_y genera lie on one straight line of the plane. | Hirzebruch's genera | [`landscape`](topics/landscape/) |

Every link in this table is known mathematics. What the box adds is one shape that holds them all, so a fact in one row suggests where to look in another. No new theorem is claimed here.

## What's here that could help mathematicians

These are the pieces most worth a working mathematician's time. The bar is strict: an item is here only if it gives something usable (a formal check, a sharp open question, a measurable number, a negative result that saves effort, or a lookup tool). Pure re-derivations of known math are left out; INDEX.md has those.

Each entry says how solid it is (**proved in Lean**, **exact computation**, **numerical**, **open**) and how new it is. "Not found in the literature" means the write-up's own search, plus a spot check on Nov 25, 2026, found nothing; it is not a guarantee. Where the result is known, the entry says so and names the source.

### Questions worth a mathematician's attention

1. **A fractional power set needs the absence of choice.** With choice there is no "dial" of power sets: the iterated halves P^(1/2), P^(1/4), P^(1/8), … would have to shrink forever, and well-ordered sizes can't. Any monotone half power set must land strictly between κ and 2^κ; where 2^κ = κ⁺⁺ everywhere (consistent, Merimovich 2007) the plain successor is one, and odd gaps have none. Without choice the obstruction disappears (Karagila 2014: decreasing chains of sizes are consistent even with DC_κ), and we found no result building the dial or ruling it out. *Why care:* a crisp ZF question with a concrete line of attack: a permutation model from a group whose orbit counts grow like e^(c√n). *Status:* with choice, proved (Lean, on a well-ordered scale; 11 theorems); without choice, open. *New?* The no-dial argument and the strictly-between bound are short and elementary, but not found stated anywhere. [`topics/fractional-powerset/`](topics/fractional-powerset/)

2. **Is there a category of signed objects whose product is a + b − ab?** Mult-add itself has an honest set model: in the category of partial functions, A ✳ B is the categorical product, with (1 + a)(1 + b) − 1 elements. Negative dial values (k = −1 is "OR", k = −2 is "XOR") have only Euler-characteristic models. Whether some category of signed or virtual objects has the k = −1 operation as its categorical product is open as far as we know. *Status:* open. [`topics/category-theory/`](topics/category-theory/)

### Tools you can use now

3. **A machine-checked library of the a + b + k·ab family.** 3,025 theorems in 120 modules, all compiling in core Lean 4 (no mathlib), none using `sorry`, each using only Lean's three standard axioms. It formalizes, among others: Sugeno λ-measures, Tsallis q-algebra, the void and variance laws of α-determinantal point processes, the diamond hierarchical lattice map, Lee–Yang ring zeros, and the bookkeeping of the Fermat proof. A port to mathlib is under way on its own branch, and a side-by-side reading of all 1,750 ported theorems found 4 real changes, all small helper lemmas (a build will confirm). *Status:* proved in Lean. *New?* The formal checks are new; most of the statements are known mathematics. [`lean-library/`](lean-library/README.md), [`topics/mathlib-port-check/`](topics/mathlib-port-check/README.md)

4. **A map of one-parameter composition laws.** Every law (a + b + k·ab) / (1 − q·ab) is a point (k, q) on one plane, and the parabola k² = 4q splits it into three regions: two real walls (saturation, like velocity addition), one double wall (addition in disguise), no real wall (turning, like tan addition). The shape number s = k²/(k² − 4q) says which fractional rung a law looks like. More than 44 laws from the other topics are placed; four empty spots turn out to hold known laws (logistic growth and Kalman variance, the Riccati log-slope of a wave, and Hirzebruch's χ_y genus among them) and one spot is provably empty. *Why care:* read off a law's walls and the coordinate that turns it into plain addition at a glance. *Status:* 65 exact and numeric checks, 27 theorems in Lean. *New?* No: the plane is exactly Coleman and McGuinness's 1991 classification of rational formal group laws (split by k² − 4q), and the three kinds go back to Lie and Klein's one-parameter Möbius groups. What the map adds is the picture and the 44 placed laws. [`topics/landscape/`](topics/landscape/)

5. **A proved map of where mult-add sits in gravity.** Mult-add's own rules build one fixed curvature (relativity's velocity space), and a dial that varies only by place can give any curvature in 2D but not in 3D or 4D, where Schur's theorem forces it to be constant. The smallest addition that works is a dial on each pair of directions, which is exactly a metric. With it, all four classical routes to Einstein's equations go through: Jacobson's heat flow across horizons (where mult-add supplies the most: area bookkeeping and focusing), teleparallel gravity, Lovelock's uniqueness and Deser's bootstrap. The river picture of a black hole is exactly a teleparallel frame, its torsion is the river's velocity gradient, and moving mass makes the river swirl. A one-number dial, the sombrero included, gives only Nordström's scalar gravity, which bends no light and is ruled out. *Status:* 65 checks pass, 22 theorems in Lean. *New?* Every bridge is known physics (Jacobson 1995, TEGR, Lovelock 1971, Deser 1970, Freund–Nambu 1968); the mult-add reading of each, and the exact point where the metric has to come in, are ours. It gives no new prediction: every route with the direction dial is general relativity exactly. [`topics/gr-bridges/`](topics/gr-bridges/), [page](https://claude.ai/artifact/LhAbnuPMJnNoq9cRYBz3xc)

6. **A template for formalizing a proof's architecture.** Fermat's Last Theorem laid out as a map of where each would-be counterexample dies, with the deep inputs (Wiles, Ribet, Mazur, the genus step) as named hypotheses: 24 Lean theorems show FLT follows from exactly those. *Why care:* a pattern for formalizing famous results top-down before the hard parts exist in a library. *Status:* proved in Lean, conditional on the named inputs. *New?* The method, not the mathematics. [`topics/fermat-counterexample-map/`](topics/fermat-counterexample-map/)

7. **Closed forms for the diamond hierarchical lattice.** Writing the fixed points of the Potts diamond map as t = 1/y with Q = (y − 1)²(y + 1) gives the slope at every fixed point as (2y/(y + 1))², so the exponent is ν = ln 2 / (2 ln(2/(1 + t*))), with exact values at Q = 1 (golden ratio), Q = 2 (tribonacci), 32/27, 3, 16, 45. *Status:* exact computation (70 checks), 48 theorems in Lean. *New?* The fixed-point cubic and the closed forms at Q = 1, 3, 16 are in Chang, Roeder and Shrock (2020); the y-parametrization and the slope formula were not found written (elementary). [`topics/diamond-lattice/`](topics/diamond-lattice/)

### Negative results that save effort

8. **Elliptic near misses.** Two irrational rungs, s = (125 ± √11665)/90 ≈ 0.189 and 2.589, are copied by a genuine elliptic formal group through degree 7 and then fail at degree 8; no rung other than addition and multiplication matches to all orders. *Status:* exact computation. [`topics/elliptic/`](topics/elliptic/)

9. **Beraha numbers do not organize the laws.** Of 85 placed laws, none sits at the non-integer colour-count values k = 1 + 2cos(2π/n); the three hits are all at k = 2, and whole numbers do as well as the Beraha values. On the diamond lattice, ν(Q) is smooth through every Beraha value. *Status:* numerical, with the integer facts in Lean. [`topics/lee-yang/`](topics/lee-yang/), [`topics/diamond-lattice/`](topics/diamond-lattice/)

### Known results, usefully reframed

Known mathematics, but the one-picture reading may help someone teaching or connecting fields. No novelty is claimed.

- **Sugeno's λ-measure is the dial**, and the existence rule for α-determinantal point processes (Shirai–Takahashi; Maunoury) is the same "wall off the range" rule: P(no point) = (1 + k·(μ₁ ✳ₖ μ₂ ✳ₖ …))^(−1/k), where ✳ₖ is a + b + k·ab and μᵢ are the kernel's eigenvalues, with the variance a straight line in k. 40 Lean theorems, 101 checks. [`topics/dial-point-processes/`](topics/dial-point-processes/), [`topics/measure-probability/`](topics/measure-probability/)
- **The half power set under joining** is unique: dial √2 − 1 on every element (odds of inclusion), with a mirror twin √2 + 1 under complement; on the primes every reasonable half gives growth X/√(log X) (Selberg–Delange). 48 Lean theorems. [`topics/fractional-powerset-dials/`](topics/fractional-powerset-dials/)
- **One family of prime products from fermions to bosons:** ∏_p (1 − k·p^(−s))^(−1/k) gives ζ(s)/ζ(2s) at k = −1 and ζ(s) at k = 1, factors into zeta values with necklace-count exponents, and has its natural boundary at Re s = 0 exactly when k ≠ ±1 (Dahlquist 1952; k = 0 is Landau–Walfisz 1920). [`topics/dial-survey/`](topics/dial-survey/)
- **Wick rotation is a change of units by i:** it moves a law from (k, q) to (i·k, −q), turning thermal tanh addition into quantum tan addition. [`topics/thermodynamics/`](topics/thermodynamics/)
- **Electromagnetism sets the dial at q = −1/c².** The wave equation survives a dial boost only there; E and B ride the boost along the light lines (scaled by the Doppler factor D and 1/D, with D² − 1 composing by mult-add); the loop phase follows gravity's loop rule, charge × flux. EM is linear and gravity is not because gauge phases commute and boosts don't. No new prediction. 41 Lean theorems, 49 checks. [`topics/electromagnetism/`](topics/electromagnetism/), [page](https://claude.ai/artifact/A71iTBsgUHZSkcbNUrEacY)
- **In Schwarzschild, the clock and radial-ruler deficits are mult-add inverses,** (−g_tt − 1) ✳ (g_rr − 1) = 0, from one boost through the Painlevé–Gullstrand "river"; it shows at a glance why light bends twice the Newtonian amount. [`topics/general-relativity/`](topics/general-relativity/)

## Promising leads

The live leads, each judged on its own merits. "Known" marks what is already established mathematics or physics; anything else is our claim and is labelled as such. This list is updated as each piece of work reports.

- **A dial that depends on direction (k as a matrix or Clifford number).** Settled as a result (item 5 above): a place-only dial is blocked in 3D and 4D, and a direction dial is a metric renamed, plus the reason a metric must appear. The Clifford deep dive adds that mult-add's dial is a null multivector, and that with multivector inputs the law stays associative while its commutator is exactly the wedge product. *Known:* Schur's theorem, Clifford and Grassmann's algebra, gauge-theory gravity (Lasenby, Doran and Gull). *Still open:* whether the light-cone dial picture is worth a pure-math write-up, and a 2D Clifford–Klein form built from ✳. No new prediction. ([GR bridges](topics/gr-bridges/), [Clifford](topics/clifford/))
- **The Jacobson chain in mult-add language.** Every geometric piece of Jacobson's derivation has a mult-add form; the one missing piece is the area law. If information-keeping could give entropy proportional to area, that would be a real step. Nothing produces it yet. *Known:* Jacobson 1995. ([write-up](topics/gr-bridges/))
- **The river as a physical field (Einstein–aether).** A testable variant, but it is existing research and the 2017 neutron-star merger already removes most of the room. ([write-up](topics/gr-bridges/))
- **Ruled out: a one-number dial for gravity, sombrero included.** It gives only Nordström's scalar gravity, which bends no light. A massive dial can hide from solar-system tests, but it never supplies the rule that matter bends space. ([write-up](topics/gr-bridges/))
- **Fractal boxes: the hidden dials are found, and one small piece is missing.** A box of fractional dimension d has 2^d corners, but many fractals share dimension ½ and still count differently at the same scale, so more dials describe how a shape is put together. The follow-up named six known ones (zoom, placement, wobble, weights, topology, branching growth). Mult-add already has the two that matter for counting: dimension is a fractional mult-add power, and weights are the per-element dial t = p/(1 − p). The one missing piece is a zoom (size) dial, which would make dimension an output, d = log(pieces)/log(1/zoom). That is a small, clean addition, not a new lead. Placement and topology are geometry that doesn't change counts. Also settled: joining two half-dimension boxes equals the √2 − 1 dial on every element, and every fractal has as many points as the continuum, so fractals alone don't touch the choice-free question in item 1 above. *Known:* all of these dials are known math (box-counting and Moran dimension, multifractal weights, lacunarity, complex dimensions). 25 + 25 checks, 15 + 14 theorems in Lean. ([first write-up](topics/fractal-boxes/), [page](https://claude.ai/artifact/NbFUcWVphAGGNvhug8XvwS); [follow-up](topics/fractal-dials/))
- **Closed: the dial k setting how many copies each element may have.** k = −1 allows each element once, k = −1/m at most m times, k = 1 without limit, which runs the prime products ∏_p (1 − k·p^(−s))^(−1/k) from ζ(s)/ζ(2s) to ζ(s). The line Re s = 0 is a natural boundary for every real k except ±1, but that is known math: a special case of Dahlquist's 1952 theorem, with k = 0 done by Landau and Walfisz in 1920. The copy-limit points k = −1/m meet the Lee–Yang values only at k = −1. What's worth keeping is the picture: each singularity's strength is a necklace count, and k = ±1 are exactly where those counts switch off. 18 + 15 theorems in Lean. ([write-up](topics/dial-survey/), [wall follow-up](topics/dial-survey/wall/), [page](https://claude.ai/artifact/LoM7iRCttrFnAi9VGz2kVr))
- **The size dial: answered, with one lead left.** The 54.5% is exact: a plain ½ on every whole number weighed by size keeps −12·Li₂(−½)/π² = 0.545206 of the total. Weighing by size adds one number, the dimension β of the set (sharp count 58.50%, whole numbers 54.52%, plane 52.42%, primes exactly 50%), and it is not a new kind of dial but the twist q = 1 + 2^(−β). What mult-add leaves over is named: repetition, which is always the same set again at double size (share 2^(−β); for the 3-D case that is the Stefan–Boltzmann 7/8), and for number fields the class group. Both close with things mult-add can express. *Known:* Fermi–Dirac functions, the Stefan–Boltzmann factor, class groups. *Lead:* a continuous "how much repetition" dial with mult-add at one end, which is Haldane's exclusion statistics in physics; known, but a useful frame. ([write-up](topics/size-dial-residual/), [page](https://claude.ai/artifact/PKMMfy7vXAZRK9cAucwSTg))
- **The dial as a rolling dark-energy field.** From the dial-potential search (71 guesses, scored against measured numbers with held-out tests): rescaling forces the dial's potential to be an exponential in ln k, and fitted to supernova and BAO data it beats a cosmological constant by Δχ² ≈ 4.8 (about 2σ) and survives cross-validation. *Known:* this is the same hint DESI reported in 2025. The best slope is 0.80, with √(2/3) inside the 1σ range, but so are 13 of 43 simple slopes, so that match means nothing yet; the natural coupling would also need screening to pass solar-system tests. *Ruled out in the same search:* the dial as the Higgs (its best version misses the Higgs mass by 3.3σ), and number formulas for the Higgs and top masses, which hit no better than random targets. ([write-up](topics/dial-potential/), [page](https://claude.ai/artifact/WpQzqeezDjmRmUcbgTA7RT))
- **Repetition: settled inside one element, with one piece outside.** How many times one element may be used falls out of the dial k itself: giving an element the weight (1 + kx)^(1/k) makes 1/k its number of cells, so k = 1/m allows at most m copies, k = 1 is mult-add's "at most once", negative k allows any number, and k can be read off the average copy count (Fermi–Dirac at k = 1, Bose–Einstein at k = −1). The catch is that between k = 1, 1/2, 1/3, … some counts come out negative. Crowding *between* elements (Haldane's exclusion statistics) is not a single k and stays outside the composition law; it would need its own dial. *Known:* every piece (Box–Cox and the Tsallis q-logarithm, the 1/(eˣ + a) occupation formula, Haldane 1991 and Wu 1994). Only the packaging (1/k as cells, the sign map, the one-dial test) is possibly new, and it hasn't been checked closely against the literature. 29 checks, 22 theorems in Lean. ([write-up](topics/repetition-dial/README.md))
- **Paper candidates (review of Nov 25).** Every catalog result was checked against the literature. The best first submission is an expository article, "One box, many laws", for *Mathematics Magazine* or the *Monthly*: known mathematics, newly connected, with every row cited. The (k, q) map matches Coleman and McGuinness 1991, so it goes in as a picture of their theorem. The one result that appears new is the fractional power set (item 1); its elementary part may be folklore, so the plan is to ask set theorists publicly first. ([review](topics/paper-candidates/README.md))

## Lean

All the Lean proofs are gathered in [`lean-library/`](lean-library/README.md): one shared library built on a single base of mult-add laws. As of Nov 25, 2026 it holds 3,025 theorems in 120 modules, every one compiling with **core Lean 4 only** (no mathlib), none using `sorry`, and each using only Lean's three standard axioms. `lean-library/THEOREMS.md` lists them all; `./build.sh` rebuilds and re-checks them.

Some topic folders also keep the original Lean file they were first written in. The one file that imports mathlib, `topics/tesseract/TesseractCoords.lean`, writes mult-add as `✳` (mathlib already uses `⊞` for biproducts) and compiles cleanly against mathlib v4.37.0-rc1 (checked on Brennan's Mac Studio, Nov 25, 2026). A port of the whole library to mathlib is under way on its own branch.

## Releases

Releases are made automatically. To cut one, add a `## vX.Y.Z` section to [CHANGELOG.md](CHANGELOG.md), set the same number in [VERSION](VERSION) and in `version:` in [CITATION.cff](CITATION.cff), and push to main. The release workflow then tags it and publishes the GitHub release, and Zenodo archives that snapshot under a new permanent DOI. Pushes that leave VERSION unchanged make no release.

## License

Copyright 2026 Brennan Pope.

- **Code** (the Lean library and every `.lean` file, the Python checks and scripts, and the page-building code) is under the [Apache License 2.0](LICENSE).
- **Everything else** (the write-ups, READMEs, notes, pages, figures and the data we generated) is under [Creative Commons Attribution 4.0 International](LICENSE-docs) (CC BY 4.0): share and adapt it freely, with credit.
- Outside data copied in for checks, such as the Pantheon+ supernova tables, keeps its own source's terms.

## Where to look

- [INDEX.md](INDEX.md) lists every topic with its folder and its published page.
- [IN-PROGRESS.md](IN-PROGRESS.md) lists folders that were still being written when this repo was made.
- [DATA.md](DATA.md) lists the large files left out and where they came from.

Each folder under `topics/` usually holds the check scripts (`checks.py`), their saved output, the figure scripts and figures, any Lean proofs, and the page source (`index.html` or `page.html`).
