# Mathlib port check: does the port prove the same things as the originals?

> **Names:** some words or symbols on this page (mult-add, fused, rung) also mean something else in standard math, physics or computing, or inside this project. [NAMING-COLLISIONS.md](../../NAMING-COLLISIONS.md) says which meaning is in play.

Checked Nov 25, 2026, against branch `mathlib-port` of BrennanRPope/fused-hyperoperation-algebra
(folder `lean-library-mathlib/`, head `95bf1f4`), compared with the core-Lean originals it was
ported from (`lean-library/` as of main `7a51396`). Those 62 original files are byte-identical to
the current library in the project folder, so nothing ported has changed since.

**How this was checked.** Every one of the 1,750 original theorems was paired with its port
counterpart using the porter's name tables, and both statements were read side by side, including
the `variable` lines, typeclass assumptions and every definition they use. Nothing was compiled:
there is no Lean on the machine that did this check, and the Mac Studio was not contacted. So every
verdict below is by careful reading, and the few items marked "check in Lean" are where a build
would settle a doubt. Nothing on the branch was changed.

## Answer

The port says what the originals said. Of 1,750 original theorems:

| Verdict | Count | Meaning |
|---|---:|---|
| SAME | 564 | Same statement up to renaming or notation |
| EQUIV | 114 | Equivalent restatement (e.g. `u⁻¹` for a hypothesis `x * u = 1`, a conjunction split into lemmas) |
| GENERAL | 715 | Port is more general or stronger, and the original follows by specialising (ℚ → any field, hypothesis removed) |
| MATHLIB | 291 | Replaced by a mathlib lemma that states the same thing |
| DROPPED-OK | 62 | Dropped as redundant, and something really covers it |
| WEAKER | 2 | Port proves less |
| MISSING | 2 | No port statement covers it |
| DIFFERENT / VACUOUS | 0 | |

So 678 match (SAME 564 + EQUIV 114), 1,068 more are covered by a generalisation, a mathlib lemma
or a redundant drop, and 4 are real changes. All four are small helper lemmas, and nothing
downstream in the port depends on what they lost.

## Real changes, worst first

1. **`FusedPeano.add_right_cancel'` is MISSING.** Original: on any fused-Peano model `N`,
   `add a c = add b c → a = b`. The table maps it to mathlib's `add_right_cancel`, but that needs an
   `IsRightCancelAdd N` instance, and the port only transports a `CommSemiring` to `N`, which does not
   give one (the porter's own note says no such instance is declared). It is true and two lines from
   `ringEquivNat.injective`, but the port no longer states it. This is the only finding in the
   foundation modules. Fix: add the lemma, or an `IsCancelAdd N` instance.
2. **`Fermat3.measure_lt` is MISSING.** Original: for integers with `x, y ≠ 0` and `N1³ < |z³|`,
   `N1 < |xyz|`. The table says `fermatLastTheoremThree` covers it, which is wrong: its hypotheses can
   hold, so FLT3 does not make it trivial. It was bookkeeping for your hand-written Euler descent,
   which the port drops in favour of mathlib's FLT3. Fix: restate it, or accept the drop with the descent.
3. **`FreyPieces.pow_dvd_bound` is WEAKER.** Original: for any `q ≥ 2` with `q ∤ m`,
   `q^a ∣ q^b · m → a ≤ b`. The mathlib replacement (`Nat.Prime.pow_dvd_iff_le_factorization`) needs
   `q` prime. The original was only ever used with primes, and that use is ported, so nothing breaks.
   Fix: restate it for `2 ≤ q`.
4. **`OsgoodInverse.const_block` is WEAKER.** Original: a constant `c` summed over a block of `2^j`
   terms is `2^j · c` for any starting index `s`. Port (`blk_const`): only for the block starting at
   `2^j`, though over any field instead of ℚ. The porter flagged this narrowing; only that start was
   ever used. Fix: one line from `Finset.sum_const` and `Nat.card_Ico`.

## Things to check in Lean (names only, statements are fine)

These rows cite a mathlib name that the reviewers could not confirm exists in the pinned mathlib.
None is used by the port, so the build passing does not test them. One `#check` each settles it:
`Nat.add_one_mul_choose_eq` (Cube.absorb), `Nat.multiplicity_mul` (KummerWilson Val.mul, low
confidence), `Int.ModEq.sum` (bsum_meq), `Nat.two_pow_succ` (OsgoodInverse pow2_succ),
`IsIdempotentElem.iff_eq_zero_or_one` (Statistics ind_cases), `neg_one_pow_eq_ite`
(pw_neg_one_parity), `fwdDiff_iter_comp_add` (Calculus dif_shift).

Two build questions also can't be answered by reading: Partitions `triple_product_30` and
`jacobi_four_squares_30` now make `decide +kernel` evaluate `IsSquare` and a sum over
`Nat.divisors`; the branch reportedly builds, which would mean they do.

## Porter's table rows that are wrong or loose (no statement lost)

- `Fermat3.primitive_of_sol`: the named mathlib lemma reduces to gcd(a, b, c) = 1, not gcd(x, y) = 1.
  The original still follows from FLT3 because its hypotheses can't hold.
- Helpers inlined with no lemma in the port, though the row names a "replacement": Fermat3
  `pow_mod_two`, `one_add_div_mul`, `Eis.small_of_sq_le_one`; KummerWilson `blockCount_eq`;
  Calculus `dif_one`; Infinities `evalL_addL` and `layers_eval` (the port never proves your recursive
  Pascal rows equal `Nat.choose`; true, but unstated).
- `le_of_dvd_two_pow_odd` maps to a mathlib lemma that gives only half of it; the other half is
  proved inline in `level_three_pow`.
- LinearAlgebra and JetFields rows marked "mathlib" actually point to the port's own `Basic.lean` lemmas.
- TwistedBoosts: 25 theorems moved from pairs of rationals to ℂ. The pair arithmetic is exactly
  complex arithmetic, but the port has no lemma that carries the ℂ statements back to rational
  pairs, so the original's form follows only by a cast argument.

## Where your own proofs were replaced

The statements survive, but in these places the port proves them by citing mathlib instead of your
argument. Your originals are untouched on main.

- **Fermat:** Euler's FLT3 descent (about 190 lines, `step` and bookkeeping) and the FLT4 descent
  are gone; the headline theorems come from mathlib's `fermatLastTheoremThree` and
  `not_fermat_42`. `key_cube` (Euler's cube lemma) and `case1` are still your proofs, but your
  hand-written Bézout is replaced by mathlib's Euclidean-domain machinery.
- **Wilson, Lucas, Fermat's little theorem, Kummer, the freshman's dream:** statements identical,
  proofs now mathlib's. The pieces of your cube route to Wilson survive as separate lemmas, but
  nothing in the port puts them together; the docstring describes the route without proving it.
- **FusedPeano:** you built + and × on a model from the fused operations and derived the Peano laws;
  the port copies ℕ's semiring across the bijection instead, so `add_assoc'`, `mul_succ'` and the
  like hold by construction. The operations are pinned to yours by `succ_add_eq_sadd` and `madd_eq`.
- **Cube:** all seven theorems are mathlib lemmas now, including your 30-line `layer_step`.
- **Elsewhere:** `sylvester2`/`sylvester3` (your 3×3 expansion), the binomial theorem in Infinities,
  E and NSContinuum, Calculus `newton`, `qnum_lift` (now via the geometric sum, not your mult-add
  lift), `addPow_eq`, and a few dozen one-liners. Each module's file in `mathlib-port-check/` marks
  these rows "route:".

## Modules not ported yet

The port was made from the 62-module library of main `7a51396`. The library has since grown to 103
modules and 2,648 theorems, so **41 modules (898 theorems) have no mathlib version**:
BinaryPulsar; Biology (11: Behaviour, CellBiology, Development, Ecology, Evolution, Genetics,
Metabolism, MolecularBiology, NeuroScience, Physiology, Population); CategoryTheory; Chemistry
(10: Acids, Biochem, Bonding, Electro, Equilibrium, Gas, Gibbs, Organic, Rates, Stoichiometry);
Combinatorics; ComplexAnalysis; DialPointProcesses; DiamondLattice; DiffGeometry; Dynamics;
EllipticCounting; FunctionalAnalysis; GameTheory; GraphTheory; GroupsLie; InformationFull;
KnotTheory; LogicSets; MeasureProbability; Numerical; RuleBreakers; SeveralWaves.

## Per-module results

Detailed rows (one per original theorem, with verdict, confidence and reason) and notes are in
`mathlib-port-check/<Module>.tsv` and `.md` next to this file.

| Module | Theorems | SAME | EQUIV | GENERAL | MATHLIB | DROPPED-OK | WEAKER | MISSING |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Base | 46 | 20 | · | 15 | 9 | 2 | · | · |
| Cube | 7 | · | · | · | 7 | · | · | · |
| Foundations.FusedPeano | 65 | 32 | 7 | · | 19 | 6 | · | 1 |
| Foundations.Independence | 7 | · | 4 | · | 3 | · | · | · |
| Foundations.MultAddRing | 51 | 28 | 2 | · | 18 | 3 | · | · |
| Rungs | 36 | 31 | 3 | · | 2 | · | · | · |
| Topics.AbstractAlgebra | 62 | 30 | 10 | 5 | 13 | 4 | · | · |
| Topics.Calculus | 29 | 18 | 1 | 3 | 3 | 4 | · | · |
| Topics.Chemistry.ChemConservation | 10 | 2 | · | 8 | · | · | · | · |
| Topics.Chemistry.ChemPeriodic | 12 | 8 | · | 4 | · | · | · | · |
| Topics.Chemistry.ChemQuantum | 8 | 5 | · | 3 | · | · | · | · |
| Topics.ClosingRate | 7 | 5 | · | 2 | · | · | · | · |
| Topics.Consensus | 27 | 9 | · | 16 | · | 2 | · | · |
| Topics.DarkSector | 17 | 4 | · | 10 | · | 3 | · | · |
| Topics.DataCenter | 16 | · | · | 15 | 1 | · | · | · |
| Topics.DialSweep | 28 | 14 | 1 | 11 | 2 | · | · | · |
| Topics.E | 37 | 5 | · | 16 | 16 | · | · | · |
| Topics.EarlyScience | 17 | 5 | · | 12 | · | · | · | · |
| Topics.Elliptic | 22 | 1 | 5 | 16 | · | · | · | · |
| Topics.EllipticS3 | 17 | 1 | 7 | 7 | 1 | 1 | · | · |
| Topics.EntropyNoise | 41 | 2 | 4 | 28 | 7 | · | · | · |
| Topics.Fermat.Fermat3 | 108 | 17 | 16 | 4 | 64 | 6 | · | 1 |
| Topics.Fermat.Fermat4 | 9 | 5 | · | · | 3 | 1 | · | · |
| Topics.Fermat.FermatRoute | 9 | · | 3 | 1 | 4 | 1 | · | · |
| Topics.Fermat.FreyPieces | 46 | 4 | 16 | 11 | 12 | 2 | 1 | · |
| Topics.Finance | 16 | · | 2 | 10 | 4 | · | · | · |
| Topics.FractionalRungs | 8 | 5 | · | 2 | · | 1 | · | · |
| Topics.FusedConstants | 20 | 12 | 1 | 4 | 2 | 1 | · | · |
| Topics.GeometryTrig | 24 | 3 | · | 21 | · | · | · | · |
| Topics.GravityProbeB | 16 | 4 | · | 12 | · | · | · | · |
| Topics.Infinities | 12 | 3 | 2 | 5 | · | 2 | · | · |
| Topics.InformationTheory | 52 | 1 | · | 48 | 2 | 1 | · | · |
| Topics.JetFields | 46 | 33 | 2 | 2 | 7 | 2 | · | · |
| Topics.Landscape | 27 | · | · | 27 | · | · | · | · |
| Topics.LeeYang | 39 | 32 | 1 | 2 | 3 | 1 | · | · |
| Topics.LeeYangCluster | 10 | 7 | · | 2 | · | 1 | · | · |
| Topics.LeeYangSpill | 10 | 5 | 2 | 2 | 1 | · | · | · |
| Topics.LinearAlgebra | 51 | 17 | 1 | 23 | 9 | 1 | · | · |
| Topics.NSContinuum | 18 | 10 | 1 | 6 | · | 1 | · | · |
| Topics.NavierStokes | 6 | 1 | 1 | 4 | · | · | · | · |
| Topics.NumberPlane | 26 | · | · | 25 | 1 | · | · | · |
| Topics.OsgoodInverse | 53 | 20 | · | 16 | 14 | 2 | 1 | · |
| Topics.Patterns | 32 | 22 | 4 | 3 | · | 3 | · | · |
| Topics.Pi | 45 | 21 | 1 | 10 | 8 | 5 | · | · |
| Topics.QAcross.Knots | 19 | 18 | · | 1 | · | · | · | · |
| Topics.QAcross.Partitions | 13 | 13 | · | · | · | · | · | · |
| Topics.QAcross.QCounting | 13 | 9 | 2 | 1 | 1 | · | · | · |
| Topics.QAcross.QDial | 16 | 6 | · | 10 | · | · | · | · |
| Topics.QuantumNoise | 21 | 4 | · | 17 | · | · | · | · |
| Topics.Relativity | 30 | · | · | 30 | · | · | · | · |
| Topics.Rulecraft | 10 | 7 | · | 2 | 1 | · | · | · |
| Topics.RungShape | 17 | 13 | 2 | · | 2 | · | · | · |
| Topics.ShapeDirections | 36 | 9 | 2 | 24 | · | 1 | · | · |
| Topics.SolitonsRungs | 26 | 2 | 1 | 22 | 1 | · | · | · |
| Topics.Statistics | 98 | · | 2 | 88 | 7 | 1 | · | · |
| Topics.Thermodynamics | 68 | · | · | 65 | 2 | 1 | · | · |
| Topics.Topology.CubeProofs | 23 | 4 | 1 | 2 | 15 | 1 | · | · |
| Topics.Topology.HopfInvariantOne | 17 | 8 | 1 | · | 7 | 1 | · | · |
| Topics.Topology.KummerWilson | 30 | 9 | 2 | 7 | 12 | · | · | · |
| Topics.Toric.ProjectivePlane | 16 | 9 | 2 | 1 | 4 | · | · | · |
| Topics.Toric.ToricCharts | 16 | 11 | 2 | 3 | · | · | · | · |
| Topics.TwistedBoosts | 36 | · | · | 31 | 4 | 1 | · | · |
| **Total** | 1750 | 564 | 114 | 715 | 291 | 62 | 2 | 2 |
