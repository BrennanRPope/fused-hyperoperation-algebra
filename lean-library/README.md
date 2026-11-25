# Fused Hyperoperation Algebra: the shared Lean library

One library for every Lean proof in the project, all built on one base of
mult-add laws, so each new proof starts from the same foundation.

**Status (Nov 25, 2026):** all 52 modules compile with core Lean 4.37.0-rc1
(no mathlib). 1541 theorems, none uses `sorry`, and every one uses only
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
| `Topics/Fermat/Fermat4` | fermat/lean/Fermat4.lean | |
| `Topics/Fermat/Fermat3` | fermat/lean/Fermat3.lean | |
| `Topics/Fermat/FreyPieces` | fermat/lean/FreyPieces.lean | |
| `Topics/Fermat/FermatRoute` | fermat/lean/FermatRoute.lean | imports the three above |
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
| `Topics/Relativity` | relativity/lean/Relativity.lean | 4 twins removed |
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
