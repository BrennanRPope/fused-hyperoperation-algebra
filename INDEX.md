# Index

Every topic, its folder in this repo, and its published page on claude.ai (pages open for Brennan's account). See IN-PROGRESS.md for folders that were still changing at the last sync, and PAGES-INDEX.md for the raw list of page ids.

## Drafts

Work saved mid-way so nothing is held back. Each is marked DRAFT at the top of its README.

| Topic | Folder | Page |
|---|---|---|
| DRAFT: the height flow read as a running coupling. Early verdict: behaves exactly like a marginal coupling, but so does any map with a double fixed point; no physics test found. Merged into the height-flow stream (`topics/height-flow/`). | `topics/height-running-coupling/` | — |
| Height flow: entropy family, running coupling and invariants. Invariants part DRAFT: the two-loop number of the height-s map is exactly 1/(3s), and the two half-height maps differ by the horn coefficient times \|sin πs\|. Entropy family is a modest catalogue entry (inside Tempesta's group-entropy framework); running coupling has no physics test; telling maps apart by the horn map is the known Écalle–Voronin theorem. | `topics/height-flow/` | — |
| DRAFT: Lucas corner weights U_m(k, d) for the two-dial plane law. Early verdict: the weights are known Lucas sequences and the corner formula is probably folklore. Merged into the corner-weights stream. | `topics/lucas-weights/` | — |
| DRAFT: which parts of the cube can carry a weight? Square faces as a reordering cost give soliton pair weights (known vertex-operator picture; cube framing ours) | `topics/cube-weights/` | — |
| DRAFT: the holographic entropy cone as the cube corner sum at k = −1 (physics published; cube reading ours; one open lead: a softened min cut keeps monogamy) | `topics/holographic-entropy-cube/` | — |
| DRAFT: is probability entanglement? Time, small loops, the helix and the union rule. Mostly known results; one lead: the AdS₃ cross-ratio read as a clock's tick probability, which explains why ∘₋₁ appears there. A hidden loop or helix still hits Bell's wall. | `topics/probability-entanglement/` | — |
| DRAFT: how quantum paths add up to everyday laws, and path sums as database aggregates. The circle operation does the bookkeeping along a path, but the canceling that gives Newton's laws is ordinary addition, so no new physics. Two tidy pieces: one dial map of path aggregates (counting, least action, quantum, statistical), and the height flow as a critical family tree. | `topics/quantum-to-classical/` | — |

## The library and tools

| Topic | Folder | Page |
|---|---|---|
| Shared Lean library: 3,176 theorems, 133 modules | `lean-library/` | [page](https://claude.ai/artifact/WuSVpapjX9aiXUd6juT22M) |
| Mathlib version of the Lean library: the 62 modules of the original library (1,750 theorems) on Mathlib's `PreQuasiregular`, plus the proof-complexity study. Built by CI | `lean-library-mathlib/` | see `PORT-NOTES.md` |
| Proof engine | `topics/proof-engine/` | [page](https://claude.ai/artifact/CZM59HVTwKMrK5Qnd7p2g9) |
| Proof shape gallery | `topics/proof-gallery/` | [page](https://claude.ai/artifact/NWREn3Fbw1i5reUMrYaAXZ) |
| Proof complexity: proofs through the circle operation are not shorter (negative result; its Lean files are built in `lean-library-mathlib/`) | `topics/proof-complexity/` | — |
| Shape directions map | `topics/shape-directions/` | [page](https://claude.ai/artifact/XJxHp4hRkWYBJ45ZWA7pqe) |
| Mathlib port check: does the Mathlib version prove the same 1,750 theorems? (read side by side) | `topics/mathlib-port-check/` | no page |
| Common patterns pass | `topics/patterns/` | no page yet |
| The theory's shape: the Boolean lattice (n-cube); a chain of circle operations sums over its corners | `topics/boolean-lattice/` | — |
| Corner weights: every corner is a principal minor of one matrix; Cauchy matrices give soliton pair weights (lead; frame ours; parts DRAFT) | `topics/corner-weights/` | — |
| Do shapes depend on shapes? Proof bodies order yes, whole laws no | `topics/shape-dependencies/` | not yet published |
| Mathlib port: what it takes to merge (plan; no Mathlib pull request yet) | `topics/mathlib-plan/` | — |
| Dependency map: prerequisites for a G dial, how results depend on each other, and what each rests on ([DEPENDENCY-MAP.md](DEPENDENCY-MAP.md)) | `topics/dependency-map/` | — |
| Operations catalog: names for the moves we keep using; five shorthands picked for the 1.0 rename, not applied yet | `topics/operations-catalog/` | — |

## Foundations

| Topic | Folder | Page |
|---|---|---|
| Axioms | `topics/axioms/` | [page](https://claude.ai/artifact/PrF6u5bRhkq6tPBcD5K4i2) |
| Algebra atlas | `topics/abstract-algebra/` | [page](https://claude.ai/artifact/DTBbDjwUepJAUmU9AQEAhV) |
| Linear algebra | `topics/linear-algebra/` | [page](https://claude.ai/artifact/6MmBPSMK179Q7ZRFTjz3L3) |
| Calculus from the circle operation ∘₁ (called mult-add before 1.0) | `topics/calculus/` | [page](https://claude.ai/artifact/VBRCmFrQVwQiXjMxZkY5SY) |
| Fields of jets | `topics/jet-fields/` | [page](https://claude.ai/artifact/Lg84gxR34iuBqZwQNZWjF3) |
| Osgood inverse | `topics/osgood-inverse/` | [page](https://claude.ai/artifact/TrRisNESnD3ZPv6sHbK3Gy) |
| Euler's identity (Claude Doc) | (page only) | [page](https://claude.ai/artifact/PwjUtELeuJQQbcx6vesXit) |
| e from the circle operation | `topics/e/` | [page](https://claude.ai/artifact/TBpb2hxxUxQsFwASbnAx65) |
| pi and other transcendentals | `topics/pi/` | [page](https://claude.ai/artifact/SsQ4CUfuyWHmmBSwvKYxYA) |
| Natural constants, fused | `topics/fused-constants/` | [page](https://claude.ai/artifact/TW1FbQ74pbLGD9P68doQ2y) |
| Fractional and negative rungs | `topics/fractional-rungs/` | no page yet |
| e and the supernatural logarithm at every rung height: e_s = 1 + E^s(1) | `topics/fractional-rungs/e-every-height/` | — |
| The smooth shape through the rungs (rung sheet) | `topics/rung-shape/` | [page](https://claude.ai/artifact/AHr9NVcjzr4cHjQVrRaEnJ) |
| Infinities | `topics/infinities/` | [page](https://claude.ai/artifact/YEchSB4Ft97HsNEcV4q3aJ) |
| Landscape map of fused laws | `topics/landscape/` | [page](https://claude.ai/artifact/E2ri2CF3xiGEWA5B9hAWYf) |
| The landscape from the simplest laws up | `topics/simple-laws/` | [page](https://claude.ai/artifact/KjLF8S4j6DKkgUL5bUrNr9) |
| Circle-operation crawl map | (page only) | [page](https://claude.ai/artifact/PgWqA5hKzRKzWJDSr3FCBM) |
| Tesseract coordinates | `topics/tesseract/` | [page](https://claude.ai/artifact/KV9thET91iskRFkxqPne5n) |
| Toric geometry and the primes | `topics/toric/` | [page](https://claude.ai/artifact/PV2A6FTHF98FTqBvHsCazZ) |
| Topology (K-theory, Adams operations, Hopf invariant one) | `topics/topology/` | [page](https://claude.ai/artifact/4Nb6GyQBxi8gDNdHbszA26) |
| Statistics from the circle operation | `topics/statistics/` | [page](https://claude.ai/artifact/FwYdFUqfFh4DEAjHMqNEsv) |
| Information theory | `topics/information-theory/` | [page](https://claude.ai/artifact/SBsWmduwSNADpPVbeRDq2v) |
| Fractional power sets and choice grades | `topics/fractional-powerset/` | [page](https://claude.ai/artifact/94mPKskoy8QKXDumoDjRb7) |
| Power set dials: one dial per element (merging, complement, primes, other sets) | `topics/fractional-powerset-dials/` | [page](https://claude.ai/artifact/27MDnoDnLAMxUUmWGADYk9) |
| Dial survey: the five dials, and k acting on per-element dials and on the primes | `topics/dial-survey/` | [page](https://claude.ai/artifact/LoM7iRCttrFnAi9VGz2kVr) |
| Fractal boxes: does a box of fractional dimension give a fractional power set? | `topics/fractal-boxes/` | [page](https://claude.ai/artifact/NbFUcWVphAGGNvhug8XvwS) |
| Hidden dials of a half box: the circle operation has dimension and weights, lacks only a zoom dial | `topics/fractal-dials/` | page to come |
| What could become a paper: literature check of every result, ranked | `topics/paper-candidates/` | — |
| The repetition dial: copies per element fall out of k (1/k cells); crowding between elements stays outside | `topics/repetition-dial/` | — |
| The zoom dial: size as one new input; dimension comes out at the wall (Moran's equation) | `topics/zoom-dial/` | — |
| Composition, order and the collection laws: the circle operation as stretches about one centre, order cost as a wedge, the route to quantum, no functor half power set | `topics/order-dial/` | — |
| Why quantum needs i: the twist fixes i's sign; path integrals as the weighted circle chain at an imaginary dial | `topics/why-i/` | — |
| Path integrals in a stronger sense: the three published pieces that make the chain a real path integral | `topics/why-i/path-integral-strong/` | — |
| Feynman's checkerboard and Kac's telegraph as one weighted circle chain (k = −i and k = −1) | `topics/why-i/checkerboard/` | — |
| Adding a potential to the checkerboard and telegraph chain | `topics/why-i/checkerboard-potential/` | — |
| e, i and π from one dial, and where i lives in physics | `topics/why-i/e-i-pi/` | — |
| Does resolving the chain force the turn? Only reversibility and the classical path at every length force k = −i | `topics/why-i/e-i-pi/chain-resolves/` | — |
| Dials for G: three published dials with stable ends, and G as the square of the stability-forced length | `topics/why-i/e-i-pi/g-dials/` | — |
| Units at different levels: asymptotic safety, dimension flow and G (the RG flow of G as a circle-operation chain) | `topics/why-i/e-i-pi/asymptotic-safety/` | — |
| Beyond one loop: G's running matches the circle operation only at one loop or large N (limiting result) | `topics/why-i/e-i-pi/asymptotic-safety/beyond-one-loop/` | — |
| Height dial: e at rung height 1.10 fits 23/11 but misses every other published number (negative result) | `topics/why-i/e-i-pi/asymptotic-safety/beyond-one-loop/height-dial/` | — |
| k as the Cayley–Klein curvature dial (κ = −k²/4) | `topics/why-i/cayley-klein/` | — |
| The circle operation as κ-Poincaré momentum addition | `topics/why-i/kappa-poincare/` | — |
| The dial map: the framework does not make gravity a forced step | `topics/why-i/dial-map/` | — |
| String theory: what it rests on, G as (coupling)² × length^(D−2), T-duality as the k = 1 inverse | `topics/string-theory/` | — |
| Names and notation review before the 1.0 freeze (recommendations only) | `topics/names-review/` | — |
| Interim theory name: Hypercube Composition Theory (final name still open, ROADMAP item 8) | `topics/names-review/THEORY-NAME.md` | — |
| Electromagnetism through the circle operation: light sets the dial d = −1/c² | `topics/electromagnetism/` | [page](https://claude.ai/artifact/A71iTBsgUHZSkcbNUrEacY) |
| Size weighting as a dial, and what the circle operation leaves over (the 54.5%) | `topics/size-dial-residual/` | [page](https://claude.ai/artifact/PKMMfy7vXAZRK9cAucwSTg) |

## Fourteen branches of math

| Topic | Folder | Page |
|---|---|---|
| Category theory | `topics/category-theory/` | no page yet |
| Combinatorics | `topics/combinatorics/` | no page yet |
| Complex analysis | `topics/complex-analysis/` | no page yet |
| Differential geometry | `topics/differential-geometry/` | no page yet |
| Dynamical systems | `topics/dynamical-systems/` | no page yet |
| Functional analysis | `topics/functional-analysis/` | no page yet |
| Game theory | `topics/game-theory/` | no page yet |
| Graph theory | `topics/graph-theory/` | no page yet |
| Groups and Lie theory | `topics/groups-lie/` | no page yet |
| Clifford's work, piece by piece, through the circle operation | `topics/clifford/` | [page](https://claude.ai/artifact/NU8rLX3jBBeaojcq8g3hMW) |
| Information theory, full branch | `topics/information-theory-full/` | no page yet |
| Knot theory | `topics/knot-theory/` | no page yet |
| Logic and sets | `topics/logic-sets/` | no page yet |
| Measure and probability | `topics/measure-probability/` | no page yet |
| Numerical analysis | `topics/numerical-analysis/` | no page yet |

## Number theory

| Topic | Folder | Page |
|---|---|---|
| Fermat | `topics/fermat/` | [page](https://claude.ai/artifact/Q3GveFG3TjJSeqrYbXEouL) |
| Fermat's Last Theorem in pieces | `topics/fermat/` | [page](https://claude.ai/artifact/9wdazdAULhiQ32uXPpgKLa) |
| Tracking the Lean FLT project (Buzzard; futarchy-fi fork) | `topics/flt-tracking/` | no page yet |
| Elliptic curves lead | `topics/elliptic/` | [page](https://claude.ai/artifact/JhsmEQaHNRhZ9t5Ph5RZ6y) |
| The birth of a prime (Claude Doc) | `topics/primes/` | [page](https://claude.ai/artifact/SisvRHEDreffPYfN8mD1ch) |
| q across knots, primes, counting and physics | `topics/q-across/` | [page](https://claude.ai/artifact/9AtexUZrzvsPfzSA7g91dL) |
| Fermat counterexample map (Wiles, Ribet, Mazur and genus steps are named assumptions in Lean, not proofs) | `topics/fermat-counterexample-map/` | [page](https://claude.ai/artifact/G1nEFedgkjtv7X5ehD3ber) |
| Unsolved problems through the circle-operation lens | (page only) | [page](https://claude.ai/artifact/FkH6Zb3JS7URmP7S46HWoz) |

## Computation

| Topic | Folder | Page |
|---|---|---|
| P vs NP: what can't work, and where the circle operation touches it | `topics/p-vs-np/` | [page](https://claude.ai/artifact/FQr7rMChaa6ofsJV2WfDjs) |

## Physics

| Topic | Folder | Page |
|---|---|---|
| Common PDEs | (page only) | [page](https://claude.ai/artifact/6LbCxRLB9aTtBAxzQcCiXt) |
| Physical laws | (page only) | [page](https://claude.ai/artifact/DLhhZpKrY4FYYtAoaQxrQo) |
| Navier–Stokes | `topics/navier-stokes/` | [page](https://claude.ai/artifact/SkJga4ih1vuxS5ozxeZefF) |
| Navier–Stokes at each rung | `topics/ns-rungs/` | [page](https://claude.ai/artifact/35PaQq5t56ZxNhHjKdq8pJ) |
| Navier–Stokes across dimensions | `topics/ns-continuum/` | [page](https://claude.ai/artifact/GBDAr2D6pBnyedGtWAG8Va) |
| Solitons on the rungs | `topics/solitons-rungs/` | [page](https://claude.ai/artifact/5yXzCgEwC2jY5y6LkFCUXJ) |
| Relativity | `topics/relativity/` | [page](https://claude.ai/artifact/S2RVDVVQWR4UAQVtmDQQMq) |
| Twisted boosts and Gravity Probe B | `topics/twisted-boosts/` | no page yet |
| General relativity through the circle operation | `topics/general-relativity/` | [page](https://claude.ai/artifact/SHrK88RpTBLLfTaeAxgi37) |
| Curved space from the circle operation, and four bridges to Einstein's equations (Jacobson, teleparallel, Lovelock, Deser) | `topics/gr-bridges/` | [page](https://claude.ai/artifact/LhAbnuPMJnNoq9cRYBz3xc) |
| The dial's potential: a scored search of 71 guesses (Higgs-like dial, rolling dark-energy dial) | `topics/dial-potential/` | [page](https://claude.ai/artifact/WpQzqeezDjmRmUcbgTA7RT) |
| Lee–Yang zeros (cluster test, fill-then-spill) | `topics/lee-yang/` | [page](https://claude.ai/artifact/UQQzmkzekxe2hL4m8F8NQL) |
| Point processes on the dial (fermions, Poisson, bosons as k = −1, 0, 1) | `topics/dial-point-processes/` | page source `index.html` in the folder |
| Diamond hierarchical lattice on the landscape (Potts model, Julia sets) | `topics/diamond-lattice/` | page source `index.html` in the folder |
| Thermodynamics | `topics/thermodynamics/` | [page](https://claude.ai/artifact/SKKLwKiEg4nyYGfrw32sp5) |
| Entropy and quantum noise | `topics/entropy-noise/` | [page](https://claude.ai/artifact/1To24w7SazEuaE5Pjsj4YM) |
| Quantum noise | `topics/quantum-noise/` | [page](https://claude.ai/artifact/X6u4iNGch19w41ioqJ5VoW) |
| Dark sector | `topics/dark-sector/` | [page](https://claude.ai/artifact/5QncmGVWkfPj81QXcJchDo) |
| SSH chain checks | `topics/ssh-chain/` | no page yet |
| Data that doesn't fit | `topics/data-that-doesnt-fit/` | no page yet |

## Chemistry

| Topic | Folder | Page |
|---|---|---|
| Chemistry, all chapters (01 to 13) | `topics/chemistry/` | no page yet |
| 01 Quantum chemistry | `topics/chemistry/01-quantum/` | [page](https://claude.ai/artifact/T7nq4VDeDLgZcdPFunvDXd) |
| 02 Conservation of mass and charge | `topics/chemistry/02-conservation/` | [page](https://claude.ai/artifact/KrZP1fbkBG53nNKvcz2kYt) |
| 03 Periodic law | `topics/chemistry/03-periodic/` | [page](https://claude.ai/artifact/1noQ3igMrH77b9bZ7ocL8y) |
| 04 Bonding and valence | `topics/chemistry/04-bonding/` | [page](https://claude.ai/artifact/N7zionN5F7ABGi86xAW1nE) |
| 05 Stoichiometry and proportions | `topics/chemistry/05-stoichiometry/` | [page](https://claude.ai/artifact/F6KLdvnVUQ2VUS4rMbscMj) |
| 06 Gas laws | `topics/chemistry/06-gas-laws/` | [page](https://claude.ai/artifact/VRTU1XGxbaxTzGT3yiB5TK) |
| 07 Chemical thermodynamics (Gibbs) | `topics/chemistry/07-gibbs/` | [page](https://claude.ai/artifact/Q14CAvSfgRNYMtd5Q36cwn) |
| 08 Equilibrium and mass action | `topics/chemistry/08-equilibrium/` | [page](https://claude.ai/artifact/6MZCz3kfVLy4fKRLwXTdNn) |
| 09 Rates and Arrhenius | `topics/chemistry/09-rates/` | [page](https://claude.ai/artifact/4bv4XrP5PCNSJEbJxqLSf4) |
| 10 Electrochemistry and Nernst | `topics/chemistry/10-electrochemistry/` | [page](https://claude.ai/artifact/YJyPaYuGYCdwEPPwFJohp5) |
| 11 Acids and bases | `topics/chemistry/11-acids-bases/` | [page](https://claude.ai/artifact/HD9Qo9wtMqhYCqyTuNPoW4) |
| 12 Organic chemistry | `topics/chemistry/12-organic/` | [page](https://claude.ai/artifact/1EZy9KmDV2Zne7xvHoCMRL) |
| 13 Biochemistry | `topics/chemistry/13-biochemistry/` | [page](https://claude.ai/artifact/SUDjpM1zrqnDHdBXeMJvXq) |

## Biology

| Topic | Folder | Page |
|---|---|---|
| Biology, all chapters | `topics/biology/` | no page yet |
| 01 Molecular biology | `topics/biology/01-molecular/` | [page](https://claude.ai/artifact/D4N4o6M3e7vRJuU5NkDpLh) |
| 02 Metabolism | `topics/biology/02-metabolism/` | [page](https://claude.ai/artifact/8sJ7P1mBE6PeqqNEDg7JhQ) |
| 03 Cell biology | `topics/biology/03-cell/` | [page](https://claude.ai/artifact/MD9HyXQFA5M5ai5okJdWEX) |
| 04 Genetics | `topics/biology/04-genetics/` | [page](https://claude.ai/artifact/AzdqcvaM3ScGG28kVTS4No) |
| 05 Developmental biology | `topics/biology/05-development/` | [page](https://claude.ai/artifact/898LXAfMTPjWoS3wg5crNY) |
| 06 Physiology | `topics/biology/06-physiology/` | [page](https://claude.ai/artifact/4tvsRPgwNeGVhqQWPTZkCw) |
| 07 Neuroscience | `topics/biology/07-neuroscience/` | [page](https://claude.ai/artifact/DU5rWsu5k7VxUHVPFGmTwt) |
| 08 Evolution | `topics/biology/08-evolution/` | [page](https://claude.ai/artifact/7Bp87asyNsREZjqMeKPWhT) |
| 09 Population | `topics/biology/09-population/` | [page](https://claude.ai/artifact/BsidfaottuXfA8RPPGDVgT) |
| 10 Ecology | `topics/biology/10-ecology/` | [page](https://claude.ai/artifact/Xw4bKkWgC9u7G2EsK7vnwW) |
| 11 Behaviour | `topics/biology/11-behaviour/` | [page](https://claude.ai/artifact/6fsetfzEB6Yq6FfPCzMZnA) |
| Biology Lean handoff notes (history: all 11 files are now in the library) | `topics/biology/LEAN-HANDOFF.md` | — |
| Photosynthesis | `topics/photosynthesis/` | [page](https://claude.ai/artifact/Bed2mSeyhXJHP3eS8A2exN) |
| Mirrored and reversed twists on photosynthesis | `topics/photosynthesis-twists/` | [page](https://claude.ai/artifact/KvTeqUBzFEUjsKWWku9BBK) |
| Nitrogen and the Three Sisters | (page only) | [page](https://claude.ai/artifact/CKd8YAk66brByG8FwZkG1f) |

## Applied

| Topic | Folder | Page |
|---|---|---|
| Actuarial | `topics/actuarial/` | [page](https://claude.ai/artifact/AYgNMio88Y3VE9gtx11ipR) |
| Consensus beyond proof of work | `topics/consensus/` | [page](https://claude.ai/artifact/EmATd6vD1v9KBAYedm5p3b) |
| Data center energy | `topics/datacenter-energy/` | [page](https://claude.ai/artifact/EfcBMhVha47gSFvp1JZMSW) |
| Economics and forecasting | `topics/econ-forecasting/` | [page](https://claude.ai/artifact/BYsXUqrGcM7TR1k3mWRUB7) |
| Rulecraft | `topics/rulecraft/` | [page](https://claude.ai/artifact/LHvib7Q1KFLT22ZUHkwqYS) |
