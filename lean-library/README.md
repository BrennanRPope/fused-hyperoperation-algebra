# Fused Hyperoperation Algebra: the shared Lean library

One library for every Lean proof in the project, all built on one base of
mult-add laws, so each new proof starts from the same foundation.

**Status (Nov 25, 2026):** all 123 modules compile with core Lean 4.37.0-rc1
(no mathlib). 3072 theorems, none uses `sorry`, and every one uses only
Lean's three standard axioms. `Check.lean` verifies this on every build.

## Build

```
./build.sh          # needs `lean` 4.37.0-rc1 on PATH; about 10 seconds
```

`build.sh` compiles each module in order into `build/`, then runs
`Check.lean`, which fails if any theorem uses `sorry` or a non-standard axiom
and writes the full list to `build/theorems.tsv` (a copy is `theorems.tsv`
here). `python3 tools/theorems_md.py` turns that list into `THEOREMS.md`.
A `lakefile.toml` and `lean-toolchain` are included for `lake build` too.

To use the library in a new proof, put the file under `FusedAlgebra/Topics/`,
start it with `import FusedAlgebra.Base` and `open MultAdd`, add it to
`build.sh` and `FusedAlgebra.lean`, and run `./build.sh`.

## Layout

Foundation (imported by the topics):

- `FusedAlgebra/Base.lean`: mult-add a ⊞ b = a + b + ab on any commutative
  semiring (so Nat, Int and Rat at once): the lift 1 + (a ⊞ b) = (1+a)(1+b),
  commutative, associative, identity 0, the wall −1, the weighted family
  a + b + k·ab with its lift, scaling, flip and wall, the mult-add power and
  weighted power, list folds with the corner expansion, the inverse
  −u/(1+u) and mult-add subtraction on fields, and a few power laws.
- `FusedAlgebra/Rungs.lean`: the fused ladder on Nat (shifted add, mult-add,
  fused exponent, fused tetration, fused pentation), each rung as the rung
  below iterated, the laws that hold at rung 3 and fail at rung 4.
- `FusedAlgebra/Cube.lean`: the layers of the mult-add cube (Pascal),
  absorption, factorial.

Topics (each imports the foundation):

| Module | From | Note |
|---|---|---|
| `Topics/ClosingRate` | fused-constants/lean/ClosingRate.lean | |
| `Topics/Consensus` | consensus/lean/ConsensusFused.lean | one step updated for Lean 4.37 |
| `Topics/DarkSector` | dark-sector/DarkSector.lean | |
| `Topics/DataCenter` | datacenter-energy/lean/DataCenterFused.lean | |
| `Topics/E` | e/EFromMultAdd.lean | |
| `Topics/Finance` | econ-forecasting/lean/FinanceMultAdd.lean | |
| `Topics/FractionalRungs` | fractional-rungs/lean/FractionalRungs.lean | first compile, three proofs fixed |
| `Topics/FusedConstants` | fused-constants/lean/FusedConstants.lean | |
| `Topics/Infinities` | infinities/lean/Infinities.lean | |
| `Topics/InformationTheory` | information-theory/lean/InfoFused.lean | |
| `Topics/Landscape` | landscape/lean/LandscapeMap.lean | first compile, no fixes needed |
| `Topics/RungShape` | rung-shape/lean/RungShape.lean | first compile, no fixes needed |
| `Topics/Elliptic` | elliptic/lean/EllipticLead.lean | first compile, no fixes needed |
| `Topics/EllipticS3` | elliptic/lean/EllipticS3.lean | first compile, one step fixed |
| `Topics/Thermodynamics` | thermodynamics/lean/ThermoFused.lean | first compile, no fixes needed |
| `Topics/NumberPlane` | simple-laws/lean/NumberPlane.lean | |
| `Topics/GeometryTrig` | simple-laws/lean/GeometryTrig.lean | |
| `Topics/EarlyScience` | simple-laws/lean/EarlyScience.lean | |
| `Topics/Patterns` | patterns/lean/Patterns.lean | first compiled here; three steps fixed |
| `Topics/EllipticCounting` | patterns/lean/EllipticCounting.lean | first compiled here; fixes listed below |
| `Topics/Fermat/Fermat4` | fermat/lean/Fermat4.lean | |
| `Topics/Fermat/Fermat3` | fermat/lean/Fermat3.lean | |
| `Topics/Fermat/FreyPieces` | fermat/lean/FreyPieces.lean | |
| `Topics/Fermat/FermatRoute` | fermat/lean/FermatRoute.lean | imports the three above |
| `Topics/Fermat/CounterexampleMap` | fermat-counterexample-map/lean/CounterexampleMap.lean | first compiled here; one proof reworked |
| `Foundations/MultAddRing` | axioms/lean/Axioms/MultAddRing.lean | standalone; does not import Base |
| `Foundations/Independence` | axioms/lean/Axioms/Independence.lean | standalone |
| `Foundations/FusedPeano` | axioms/lean/Axioms/FusedPeano.lean | standalone |
| `Topics/Pi` | pi/lean/PiFromMultAdd.lean | |
| `Topics/Calculus` | calculus/lean/FusedCalculus.lean | |
| `Topics/AbstractAlgebra` | abstract-algebra/lean/AbstractAlgebra.lean | |
| `Topics/Statistics` | statistics/lean/Statistics.lean | |
| `Topics/NSContinuum` | ns-continuum/lean/NSContinuum.lean | first compiled here |
| `Topics/JetFields` | jet-fields/lean/JetFields.lean | |
| `Topics/OsgoodInverse` | osgood-inverse/lean/OsgoodInverse.lean | first compiled here |
| `Topics/LinearAlgebra` | linear-algebra/lean/LinearAlgebra.lean | first compiled here |
| `Topics/DialSweep` | linear-algebra/lean/DialSweep.lean | first compiled here |
| `Topics/EntropyNoise` | entropy-noise/lean/EntropyNoise.lean | 6 twins removed |
| `Topics/QAcross/QDial` | q-across/lean/QDial.lean | first compiled here |
| `Topics/QAcross/Knots` | q-across/lean/Knots.lean | first compiled here; class name fixed |
| `Topics/QAcross/Partitions` | q-across/lean/Partitions.lean | first compiled here |
| `Topics/QAcross/QCounting` | q-across/lean/QCounting.lean | first compiled here; class name fixed |
| `Topics/SolitonsRungs` | solitons-rungs/lean/SolitonsRungs.lean | first compiled here |
| `Topics/PvsNP` | p-vs-np (written here) | 20 theorems; first try |
| `Topics/Relativity` | relativity/lean/Relativity.lean | 4 twins removed |
| `Topics/LeeYang` | lee-yang/lean/LeeYang.lean | |
| `Topics/LeeYangCluster` | lee-yang/lean/LeeYangCluster.lean | |
| `Topics/LeeYangSpill` | lee-yang/lean/LeeYangSpill.lean | |
| `Topics/TwistedBoosts` | twisted-boosts/lean/TwistedBoosts.lean | first compiled here; one step fixed |
| `Topics/GravityProbeB` | twisted-boosts/gpb/lean/GravityProbeB.lean | first compiled here |
| `Topics/BinaryPulsar` | twisted-boosts/pulsar/lean/BinaryPulsar.lean | first compiled here |
| `Topics/Chemistry/ChemQuantum` | chemistry/01-quantum/lean/ChemQuantum.lean | 4 twins removed |
| `Topics/Chemistry/ChemConservation` | chemistry/02-conservation/lean/ChemConservation.lean | 1 twin removed |
| `Topics/Chemistry/ChemPeriodic` | chemistry/03-periodic/lean/ChemPeriodic.lean | |
| `Topics/Chemistry/ChemBonding` | chemistry/04-bonding/lean/ChemBonding.lean | |
| `Topics/Chemistry/ChemStoichiometry` | chemistry/05-stoichiometry/lean/ChemStoichiometry.lean | 2 twins removed |
| `Topics/Chemistry/ChemGas` | chemistry/06-gas-laws/lean/ChemGas.lean | 1 twin removed |
| `Topics/Chemistry/ChemGibbs` | chemistry/07-gibbs/lean/ChemGibbs.lean | |
| `Topics/Chemistry/ChemEquilibrium` | chemistry/08-equilibrium/lean/ChemEquilibrium.lean | |
| `Topics/Chemistry/ChemRates` | chemistry/09-rates/lean/ChemRates.lean | 3 twins removed |
| `Topics/Chemistry/ChemElectro` | chemistry/10-electrochemistry/lean/ChemElectro.lean | 3 twins removed |
| `Topics/Chemistry/ChemAcids` | chemistry/11-acids-bases/lean/ChemAcids.lean | 2 twins removed |
| `Topics/Chemistry/ChemOrganic` | chemistry/12-organic/lean/ChemOrganic.lean | 2 twins removed |
| `Topics/Chemistry/ChemBiochem` | chemistry/13-biochemistry/lean/ChemBiochem.lean | 3 twins removed |
| `Topics/ShapeDirections` | shape-directions/lean/ShapeDirections.lean | first compiled here |
| `Topics/LogicSets` | logic-sets/lean/LogicSets.lean |  |
| `Topics/CategoryTheory` | category-theory/lean/CategoryTheory.lean |  |
| `Topics/Combinatorics` | combinatorics/lean/Combinatorics.lean |  |
| `Topics/GroupsLie` | groups-lie/lean/GroupsLie.lean | 2 twins removed |
| `Topics/Clifford` | clifford/lean/Clifford.lean | first compile, Nov 25; imports GroupsLie |
| `Topics/GraphTheory` | graph-theory/lean/GraphTheory.lean | 4 twins removed |
| `Topics/MeasureProbability` | measure-probability/lean/MeasureProbability.lean |  |
| `Topics/ComplexAnalysis` | complex-analysis/lean/ComplexAnalysis.lean |  |
| `Topics/FunctionalAnalysis` | functional-analysis/lean/FunctionalAnalysis.lean |  |
| `Topics/DiffGeometry` | differential-geometry/lean/DiffGeometry.lean | 1 twin removed |
| `Topics/KnotTheory` | knot-theory/lean/KnotTheory.lean | 1 twin removed |
| `Topics/Dynamics` | dynamical-systems/lean/Dynamics.lean | 1 twin removed |
| `Topics/InformationFull` | information-theory-full/lean/InformationFull.lean | 2 twins removed |
| `Topics/Numerical` | numerical-analysis/lean/Numerical.lean |  |
| `Topics/GameTheory` | game-theory/lean/GameTheory.lean | 1 twin removed |
| `Topics/DiamondLattice` | diamond-lattice/lean/DiamondLattice.lean | |
| `Topics/DialPointProcesses` | dial-point-processes/lean/DialPointProcesses.lean | |
| `Topics/SeveralWaves` | shape-directions/lean/SeveralWaves.lean | first compiled here |
| `Topics/RuleBreakers` | shape-directions/lean/RuleBreakers.lean | first compiled here; recursion limit raised |
| `Topics/Biology/MolecularBiology` | biology/01-molecular/lean/MolecularBiology.lean | no twins |
| `Topics/Biology/Metabolism` | biology/02-metabolism/lean/Metabolism.lean | 2 twins removed |
| `Topics/Biology/CellBiology` | biology/03-cell/lean/CellBiology.lean | 2 twins removed |
| `Topics/Biology/Genetics` | biology/04-genetics/lean/Genetics.lean | 3 twins removed |
| `Topics/Biology/Development` | biology/05-development/lean/Development.lean | 3 twins removed, 1 kept |
| `Topics/Biology/Physiology` | biology/06-physiology/lean/Physiology.lean | 5 twins removed |
| `Topics/Biology/NeuroScience` | biology/07-neuroscience/lean/NeuroScience.lean | 3 twins removed |
| `Topics/Biology/Evolution` | biology/08-evolution/lean/Evolution.lean | 4 twins removed |
| `Topics/Biology/Population` | biology/09-population/lean/Population.lean | 7 twins removed |
| `Topics/Biology/Ecology` | biology/10-ecology/lean/Ecology.lean | 3 twins removed |
| `Topics/Biology/Behaviour` | biology/11-behaviour/lean/Behaviour.lean | 4 twins removed |
| `Topics/FractionalPowerset` | fractional-powerset/lean/FractionalPowerset.lean | first compiled here, unchanged |
| `Topics/PowersetDials` | fractional-powerset-dials/lean/PowersetDials.lean | 48 theorems: per-element dials, merging, the half, complement, known sets |
| `Topics/DialSurvey` | dial-survey/lean/DialSurvey.lean | first compile, Nov 25: k on per-element dials, odds, necklaces |
| `Topics/SizeDialResidual` | size-dial-residual/lean/SizeDialResidual.lean | 18 theorems, first compile Nov 25: the size weight as a twist, the weighted merge gap, sharp weights commute, distinct-times-doubled and alternating doubling, two sign dials |
| `Topics/DialWall` | dial-survey/wall/DialWall.lean | first compile, Nov 25: necklace exponents behind the wall at Re s = 0, k = ±1 special, the k = −1/2 cancellation |
| `Topics/DialPotential` | dial-potential/lean/DialPotential.lean | 24 theorems, first compile Nov 25: dial flip and rescaling, even and power-law potentials, the mirror well, exponential fixed point, Coleman-Weinberg and Veltman checks, conformal γ = 5/7 |
| `Topics/FractalBoxes` | fractal-boxes/lean/FractalBoxes.lean | 15 theorems, first compile Nov 25: signed corner layers of a half box (Newton's C(½,k)), two halves join to a stick, Cantor codes one to one in any base, no linear-size rule is a half of 2^d |
| `Topics/FractalDials` | fractal-dials/lean/FractalDials.lean | 14 theorems, first compile Nov 25: five half-dimension sets certified exactly, Moran's equation as the mult-add identity, the uneven set is non-lattice (no Minkowski wobble), weights are the dial t, joining keeps shapes self-similar, gasket touches where dust never does, the power set fits no shape class |
| `Topics/RepetitionDial` | repetition-dial/lean/RepetitionDial.lean | 22 theorems, first compile Nov 25: 1/k cells per element, at most m copies at k = 1/m, coefficients polynomial in k, k read off average copies, negative counts between the k = 1/m values |
| `Topics/ZoomDial` | zoom-dial/lean/ZoomDial.lean | 10 theorems, first compile Nov 25: zoom weight is a twist, zooms compose, words of depth n weigh x(s)ⁿ, the word series, dimension at the wall x(d) = 1 (Moran), half boxes hit the wall |
| `Topics/OrderDial` | order-dial/lean/OrderDial.lean | 33 theorems, first compile Nov 25: mult-add as composition of maps sharing the wall, stretch-and-shift composition, the order gap a·d − b·c, commuting iff shared centre, compounding corner, mutation algebras, quantum plane counts, anticommuting forces no repeats, no functor half power set |
| `Topics/Conventions` | names-review/lean/Conventions.lean | 4 theorems, first compile Nov 25: the weight-to-probability sign flip (weights compose by k, probabilities by −k) and related notation facts; recommendations only, nothing renamed |
| `Topics/GeneralRelativity` | general-relativity/lean/GeneralRelativity.lean | first compiled here; one false statement (infall_light) corrected before merge |
| `Topics/GRBridges` | gr-bridges/lean/GRBridges.lean | 22 theorems: Clifford boosts, sandwich mult-add for matrix weights, Schur factor, Raychaudhuri area terms, river constraint, Brans-Dicke; first compile, one witness fixed |
| `Topics/Electromagnetism` | electromagnetism/lean/Electromagnetism.lean | 41 theorems, first compile Nov 25: the wave equation sets the dial q = −1/c², light-cone parts of E and B scale by D and 1/D, D² − 1 composes by mult-add, invariants and null fields, unit phases as a ✳ ā = 0, plaquettes and discrete Stokes, Born–Infeld cap |
| `Topics/NavierStokes` | ns-rungs/lean/FusedRungs.lean (Burgers, shell model) | |
| `Topics/QuantumNoise` | quantum-noise/QuantumNoise.lean | |
| `Topics/Rulecraft` | rulecraft/lean/RulecraftFused.lean | first compile, one step fixed |
| `Topics/Toric/ToricCharts` | tesseract/ToricCharts.lean | |
| `Topics/Toric/ProjectivePlane` | tesseract/ProjectivePlane.lean | |
| `Topics/Topology/CubeProofs` | topology/lean/CubeProofs.lean | |
| `Topics/Topology/KummerWilson` | topology/lean/KummerWilson.lean | |
| `Topics/Topology/HopfInvariantOne` | topology/lean/HopfInvariantOne.lean | keeps its own mult-add on truncated series |

Infinities, InformationTheory and FractionalRungs were re-copied from their
finished files on Nov 25 (only InfoFused had changed: two proofs, now matched).
Landscape (27 theorems) and RungShape (17, updated version) were added later that day; both
compiled on the first try. Elliptic (22) followed and also compiled
unchanged; none of it needs mathlib. EllipticS3 (17, the 3-sphere section)
compiled after one call to Rat.le_total was given its type. Thermodynamics
(68) compiled unchanged, and so did the three simple-laws files (67).
The four Fermat files (172 theorems by the check's count, including 24
`@[simp]` and `@[ext]` lemmas on the Eisenstein integers that a plain
search for `theorem` misses) and three foundation files from the axioms
folder (123) compiled unchanged; only their `#print axioms` lines were
dropped, since `Check.lean` covers that. The axioms folder's `Rederive.lean`
and its own check stay in that folder: they import the whole library.
Pi (45 theorems), Calculus (29) and AbstractAlgebra (62 by the check's count)
compiled unchanged too.
Statistics (98 theorems) compiled unchanged as well.
NSContinuum (18) had never been compiled; it built on the first try with no
changes beyond its header note.
JetFields (46) compiled unchanged. OsgoodInverse (53) and LinearAlgebra (51) had
never been compiled either; both built on the first try with no changes beyond
their header notes. NSContinuum, JetFields and OsgoodInverse come after the
Topology modules in build.sh because NSContinuum imports CubeProofs; a clean
build from scratch (`rm -rf build && ./build.sh`) passes.
EntropyNoise came with 47 theorems. Six stated exactly what a library theorem
already states, so they were removed and a comment names the library theorem in
each place: tsallis_two_is_union, weights_compose and linear_entropy_parallel
(InfoFused.gini_simpson), tsallis_weight (InfoFused.tsallis), carnot_weight
(ThermoFused.carnot_series) and decimation (QuantumNoise.bit_flip). The other
41 compiled unchanged.
DialSweep (28) had never been compiled; it built on the first try with no changes
beyond its header note, and none of its theorems is a twin of a library theorem.
The four q-across files (QDial 16, Knots 19, Partitions 13, QCounting 13) had never
been compiled. QDial and Partitions built unchanged. Knots and QCounting each needed
one fix: `[CommRing R]` became `[Lean.Grind.CommRing R]`, core Lean's ring class.
No statement changed and none is an exact twin; QDial's parity_xor and failure_or
give the same truth tables as StatFused.xor_table and or_table but over Int, and
their docstrings say so. Partitions' theorems are finite checks (to q^30 or p(50)),
as its header says.
SolitonsRungs (26) had never been compiled; it built on the first try with no
changes and has no twins. Relativity came with 34 theorems; four stated exactly
what a library theorem already states and were removed, with a comment naming the
library theorem: lift_num (ThermoFused.lens_den), lift_den (ThermoFused.lens_num),
light_absorbs (ThermoFused.vel_wall_pos) and multAdd_wall (MultAdd.neg_one_absorbs).
The other 30 compiled unchanged. Relativity builds right after Thermodynamics.
LeeYang (39) compiled unchanged after Relativity; none of its theorems is an
exact twin of a library theorem.
TwistedBoosts (32 theorems, 36 by the check's count with its two @[ext] lemma
pairs) had never been compiled. One proof step changed, no statement: in
sphere_no_wall, two `Rat.le_total 0 x` calls became `(Rat.le_total : (0 : Rat) ≤ x ∨ x ≤ 0)`.
LeeYangCluster (10) compiled unchanged after LeeYang. It has no exact twins; its
ring-general `wAdd_units` is kept beside the Int-only Landscape.units_num.
GravityProbeB (16) had never been compiled; it compiled unchanged after
TwistedBoosts, every `grind` step closing, with no exact twins.
LeeYangSpill (10) compiled unchanged after LeeYangCluster; its Int helper
`sq_nonneg'` stays local like the others under Small leftovers.
The first three chemistry subjects go under Topics/Chemistry. All compiled as
sent. ChemPeriodic (12) is unchanged. Five theorems were exact twins of Base
and were removed with a comment naming the Base theorem: in ChemQuantum,
cisd_is_weight_zero (wAdd_weight_zero), exact_is_weight_one (wAdd_one),
n_fragments (mpow_lift) and fragments_assoc (mpow_add); in ChemConservation,
fugacity_mul (pw_add). That leaves ChemQuantum 8 and ChemConservation 10.
Patterns (32) had never been compiled and builds right after EarlyScience. Three
proof steps changed, no statement: step_is_law splits the pair before `grind`, and
the two parity proofs take their split from `omega` because `Int.emod_add_ediv` and
`Int.emod_two_eq_zero_or_one` are not in core Lean 4.37. No exact twins.
ShapeDirections (36) had never been compiled; it compiled unchanged with no exact
twins. Its Int star_def and weight_is_units sit beside the general multAdd_def and
LeeYangCluster.wAdd_units.
BinaryPulsar (10) compiled unchanged after GravityProbeB, with no exact twins.
The 14-branch pass added 14 modules (513 theorems as sent), built after the Fermat
modules so that their imports (Statistics, InformationTheory, the QAcross files) come
first. All compiled as sent. Twins were found by comparing every new theorem's
elaborated statement with every library theorem's, so a match is the same statement
up to variable names. Twelve exact twins were removed, each with a comment naming the
theorem that stays: GroupsLie scale_disc and shift_disc (ShapeDirections), GraphTheory
parallel_rel, or_line, union_multAdd and series_fails, DiffGeometry halfplane_dilation,
KnotTheory qnum_at_one, Dynamics chain_multAdd, InformationFull joint_mass and
shannon_engine, and GameTheory hedge_expert (a twin of kelly_many in the same file).
Numerical keeps multAdd_bound, the same statement as ComplexAnalysis.corner_bound,
because a later proof there uses it and the file does not import ComplexAnalysis; its
Rat helpers sq_nonneg_rat and pow_nonneg_rat stay local as under Small leftovers.
Chemistry subjects 4 to 6 compiled as sent. ChemBonding (23) is unchanged. ChemStoichiometry
lost yield_lift (StatFused.de_morgan) and any_heavy (loss_n in the same file), leaving 15;
ChemGas lost mayer_list (ThermoFused.fermion_lift), leaving 12.
Chemistry subjects 7 to 9 compiled as sent. ChemGibbs (11) and ChemEquilibrium (11) are
unchanged. ChemRates lost first_order_compose (ThermoFused.ratio_series), qlog_zero
(ThermoFused.partition_square) and qlog_two (InformationFull.qlog_neg_one), leaving 8.
DiamondLattice (48) and DialPointProcesses (40) compiled unchanged after GameTheory, with no
exact twins; DialPointProcesses keeps its Rat helper sq_nonneg' local.
SeveralWaves (23) had never been compiled; it compiled unchanged, the four `rule4` checks
included, with no exact twins.
RuleBreakers (14) had never been compiled and builds after SeveralWaves. Its seven symbolic
`grind` proofs on the expanded rule4 hit Lean's default recursion limit, so the file sets
`maxRecDepth 100000`; no proof step or statement changed otherwise. No exact twins.
Later, paw_obeys was widened by its author from one weight to four (rule4 x y 1 z 1 w = 0);
it compiles with the same setting.
Chemistry subjects 10 to 13 compiled with no proof changes; ten theorems were exact twins.
ChemElectro lost ratio_series and bv_one_sided (InfoFused.campbell) and stage_losses
(ChemStoichiometry.loss_two), leaving 9. ChemAcids lost chained_K (InfoFused.campbell) and
indep_groups (ChemStoichiometry.loss_two), leaving 7. ChemOrganic lost hammett_compose
(InfoFused.campbell) and route_losses (ChemStoichiometry.loss_two), leaving 9. ChemBiochem
lost exclusive_loewe (InfoFused.shannon_weight), proofreading (InfoFused.campbell) and
bliss_fraction (ChemStoichiometry.loss_two), leaving 8.
Biology, eleven files with 164 theorems, compiled with no proof changes and builds after
RuleBreakers; MolecularBiology also imports LeeYang. 36 theorems were exact twins and are
removed with a comment naming the survivor, leaving 128. Every file restated the ✳
definition as star_def, so the copy in MolecularBiology stays and the other ten go. The rest
mostly restate Metabolism.limitation_pow (the "or" power), ChemQuantum.fragments_multAdd,
CellBiology.bliss_or and MolecularBiology.mwc_R. Development keeps pow_mul', a twin of
Dynamics.pow_mul_R, because its later proofs use it and it does not import Dynamics.
FractionalPowerset (11) was written without a toolchain and compiled unchanged, with no
imports and no exact twins; it builds after the biology files.
GeneralRelativity (28) imports only Base and builds last. Written and compiled in the same
session: the dial's boosts close on −q·rotation, PPN fractions for dial-only theories, the river
facts (deficits are mult-add inverses, horizon at the wall, β = γ = 1 remainders, Friedmann as
Newton, Gravity Probe B's factor), and the b = −a/2 step behind G = R − ½Rg.
CounterexampleMap (24 theorems) was written without a toolchain and builds after
FermatRoute. One proof changed: germain_residues used decide on a nested bounded ∀, which
core Lean cannot decide, so it now checks the same table as a Bool over List.range 11 and
reads it back; the statement is unchanged. No exact twins. The theorem count lists 28
names for it, because the four fields of its DeepArrows structure (wiles, ribet, mazur,
dim_eq_genus) are Props. They are named hypotheses, not proofs: flt_of_deep_arrows proves
Fermat's Last Theorem only from them.
EllipticCounting (19) had never been compiled and builds after Patterns. Proof fixes:
lucas_shift_one's step case gives grind the two recurrence steps, lucas_eds_unit uses
MultAdd.one_pow', and disc needed a `;` after its `let b8` line, because the next line's
leading minus was read as part of b8. One statement was false as written and was corrected
to match its own docstring: order_five_point's last conjunct had w3 * (w3 * w3 * w3), which
is W3⁴ and evaluates to −2; it now has w3 * w3 * w3 (W3³), and the conjunct is 0 as claimed.

## What was merged

The 19 source files held 449 theorems that compile. About 100 of them were
the same few laws proved again in each file (the mult-add lift, commutative,
associative and identity laws, the weighted family, powers, the inverse, the
cube layers, `pw`, a few power helpers). Each of those now exists once, in
Base, Rungs or Cube. Where a topic gave a law a real-world name (Fisher's
rule, the MOND inverse, the union rule), the name stays as a one-line reading
of the base law.

## Not in the library

- `tesseract/TesseractCoords.lean` (14 theorems) imports mathlib, which this
  environment cannot download. It has never been compiled. Porting it to core
  Lean means replacing `Finset` sums and products.
- `proof-engine/lean/FusedPictures.sketch.lean` is a design sketch of the
  geometric proof engine. By design it imports the library, not the other way.
- Pages whose results exist only as Python checks: actuarial,
  navier-stokes (simulations), primes, toric (point counts),
  photosynthesis-twists, fractional-rungs `symbolic_proofs.py` (42 sympy
  identities, 10 of them also proved here), and the loose scripts
  `photosynthesis_calc.py`, `ssh_chain_checks.py`, `anomaly_checks.py`.

## Small leftovers

A few number-type helpers stay in their topic because they are about one
type only: `sq_nonneg'` (Rat in Finance, Int in Rulecraft), the Rat and Int
binomial sums in E and Topology, `FusedConstants.mult_add` (the rung-2 case of
its general hyperoperation), and the topology `layer_pos`.
