# Tracking the Lean FLT project

Checked Nov 25, 2026. Read-only survey: nothing was posted, filed or proposed on either project.

**Short answer.** Buzzard's FLT project still rests on one big open step (the "B4" boss: the Frey
curve's p-torsion is reducible), and almost everything our Fermat files prove is already in
mathlib or in their repo. The live leads are small and concrete: borrow their minimality criterion
to close the one honest gap in our mathlib port, use their reduction theory to place our circle-operation
node law exactly inside FLT (split multiplicative reduction and the Tate curve), and, if Brennan
wants a public foothold, a few elementary open lemmas upstream. A second, much larger
AI-agent-run fork (futarchy-fi/FLT) is where most of the new work is happening, and it already
covers Frey semistability.

## 1. Where it stands

### Upstream: ImperialCollegeLondon/FLT (led by Kevin Buzzard)

| | |
|---|---|
| Head | `0198058`, Nov 18, 2026 |
| Lean / mathlib | Lean `v4.35.0`, mathlib `56a0058` |
| Size | 260 Lean files under `FLT/`; 65 `sorry` tokens in 24 files |
| What `flt` rests on | `#print axioms` in `FermatsLastTheorem.lean`: `[knownin1980s, sorryAx, …]` |
| Funding / goal | EPSRC grant to Sept 2029: FLT reduced to results known by 1989, i.e. `sorryAx` gone, `knownin1980s` and the named axioms allowed |
| Activity | 29–48 commits a month in May–Aug 2026; since September almost only mathlib bumps (4 in Oct, 3 in Nov) |

The proof spine is `FLT/Proof.lean`, a chain of "boss" statements:

- **B1** FLT ⇐ **B2** FLT for primes p ≥ 5 (proved; uses mathlib's FLT for 3 and 4).
- **B2** ⇐ **B3** "there is no Frey package" (proved, elementary: `FLT/FreyCurve/FreyPackage.lean`).
- **B3** ⇐ **B4** "the p-torsion of the Frey curve is reducible" (proved, given Mazur/Serre
  irreducibility, which is `knownin1980s` in `FLT/FreyCurve/Mazur.lean`).
- **B4** is `sorry`. This is the whole of the deep mathematics. Their route (Richard Taylor's,
  Khare–Wintenberger style) avoids Ribet: show the Frey p-torsion is "hardly ramified"
  (`FreyCurve.torsion_isHardlyRamified`, sorry), then show every hardly ramified mod-p
  representation is reducible, via lifting, compatible families, and the 3-adic and mod-3 cases
  (`FLT/GaloisRepresentation/HardlyRamified/{Lift,Family,ModThree,Threeadic}.lean`, each sorry).

Assumed results, each a single `axiom` with a literature reference (`FLT/Assumptions/`):
Mazur's torsion bound, an Odlyzko discriminant bound, plus the catch-all `knownin1980s`. Planned
but not yet statable: Galois representations attached to quaternionic automorphic forms,
automorphic induction GL1→GL2, cyclic base change, Jacquet–Langlands, Poitou–Tate, Moret-Bailly,
and a class-field-theory existence result (`FLT/Assumptions/README.md`).

What is done in Lean, beyond the spine: Frey package and curve with Δ, c₄, j and
`j_valuation_of_bad_prime`; the Tate curve construction and Tate parameter (uniformisation itself
still sorry); quadratic twists to split multiplicative reduction; finite subgroups of PGL₂; odd
irreducible ⇒ absolutely irreducible; abstract patching; Haar characters, adeles, quaternion
algebras and Hecke operators (the blueprint's mini-projects). Elliptic-curve torsion of size n²
is sorry and reserved for David Angdinata (division polynomials).

### The fork: futarchy-fi/FLT

| | |
|---|---|
| Head | `70f8b79e`, Oct 22, 2026 (PR #2192) |
| Lean / mathlib | Lean `v4.35.0-rc1`, mathlib `c32e1ec` |
| Size | 20,049 Lean files under `FLT/`; 46 `sorry` tokens in 25 files |
| What `flt` rests on | `[Mazur_statement, sorryAx, …]`: no `knownin1980s` |

An agent fleet ("crew", "beads", "proof packets") working leaf by leaf, with its own maps in
`cartography/` (status per chapter: Mazur, Odlyzko, Poitou–Tate, CFT, R = T, …). It vendors
Brasca et al.'s `FltRegular` (FLT for regular primes) to handle p < 17, so it only needs a weak
Mazur statement for p ≥ 17. Its Frey work (`FLT/FreyCurve/Serre/`, ~22,000 lines) includes
`FreyPackage.good_or_multiplicative` in `Semistable.lean`: semistability of the Frey curve over
every DVR, sorry-free. A sorry count understates what is left, because some sorries are whole
inputs (e.g. the `FLT.Assembly.*` hypotheses feeding B4). Nothing on its pages says it is being
upstreamed.

## 2. Where our work connects, honestly

Our files: `lean-library/FusedAlgebra/Topics/Fermat/{Fermat3,Fermat4,FreyPieces,FermatRoute}.lean`
(core Lean, no mathlib) and their port on branch `mathlib-port` (`lean-library-mathlib/`, head
`99b18c1`, Lean `v4.37.0-rc1`).

| Our piece | Already covered by |
|---|---|
| `Fermat3.flt3`, `case1` (cubes mod 9) | mathlib `Mathlib/NumberTheory/FLT/Three.lean` (`fermatLastTheoremThree`, `fermatLastTheoremThree_case_1`) |
| `Fermat4.fermat4`, `no_x4_y4_z2` | mathlib `Mathlib/NumberTheory/FLT/Four.lean` (`fermatLastTheoremFour`, `not_fermat_42`) |
| `FreyPieces` exponent reduction, `flt_of_dvd` | mathlib `FermatLastTheorem.of_odd_primes`, `FermatLastTheoremFor.mono` |
| normalization (a ≡ 3 mod 4, b even), `FermatRoute.NoFreyCurve` | upstream `FreyPackage` and `of_not_FermatLastTheoremFor_p_ge_5`. Our `NoFreyCurve` is their **B3** in other words |
| Frey minimal model, `frey_c4`, `frey_disc` | upstream `freyCurveInt`, `FreyCurve.c₄`, `FreyCurve.Δ` (same model, same formulas) |
| `frey_disc_valuation` (odd primes' exponent in Δ is a multiple of p) | upstream `FreyCurve.j_valuation_of_bad_prime` (the same fact for j, which is what is used) |
| `frey_semistable` (no prime divides both Δ and c₄) | fork `FreyPackage.good_or_multiplicative`, stated as real good/multiplicative reduction |
| `FermatRoute` "Ribet + X₀(2)" framing | not their route: upstream replaces Ribet by "hardly ramified ⇒ reducible" |

So there is no open item in FLT that our existing files fill. "Plug our proof into theirs" does
not apply: their spine already contains our reductions, and the open part (B4) is far from
anything we have. What is ours alone is the presentation: self-contained core-Lean proofs
(ℤ[ω] built from scratch) and the circle-operation restatements.

Can we help with an open item? Judged one by one:

- **Elementary upstream sorries** in `FLT/EllipticCurve/Torsion.lean`: `group_theory_lemma`
  (a finite abelian group whose d-torsion has size d^r for every d | n has n-torsion ≅ (ℤ/n)^r),
  the four `DistribMulAction` laws for the Galois action on points ("these should all be easy"),
  and the `Module.Finite` instance. Realistic for us with mathlib on the Mac. Caveats: the fork has
  closed most of this file already, so it may reach upstream another way; contributing means
  claiming a GitHub issue and opening a PR from Brennan's account (public, his call); and I could
  not read the issue tracker (GitHub blocks automated reading), so whether these have open,
  unclaimed issues is unchecked. `n_torsion_card` is explicitly reserved; leave it alone.
- **The assumptions** (`FLT/Assumptions/`): the README says "Anybody (humans, AI) is welcome to
  work on formalizing the proofs of these assumptions", while calling it low priority. These are
  Mazur's theorem and an Odlyzko bound: multi-month research-level formalization, not a near-term
  lead for us.
- **B4 and its parts** (hardly ramified lifting, compatible families, R = T): research-level,
  tied to Buzzard's own plan. Not something we should take on.

Contribution culture: upstream accepts AI-assisted code (there is an `FLT/Slop/` folder and
"auto-task" PRs), but work flows through claimed issues and maintainer review.

## 3. Things to reuse in our own proofs

**Status, Nov 27, 2026:** items 1 to 3 are done and build with mathlib; see [`lean/`](lean/). The
Frey model is proved minimal with good or multiplicative reduction over every DVR, and mathlib's
split condition is shown to match the circle-operation node law (away from 2: split exactly when
B − A is a square mod the prime). The torsion `group_theory_lemma` and the Galois action laws were
carried over from the fork. Nothing was sent upstream.

**Update, Nov 27, 2026 (later):** item 4 is started. Brennan asked for the Tate curve as a Mathlib
candidate; it is in [`lean-library-mathlib/MathlibCandidates/`](../../lean-library-mathlib/MathlibCandidates/)
(DRAFT, built by CI): the formal curve over `R⟦q⟧`, `c₄ = E₄`, `c₆ = −E₆`, `Δ = q · unit`, split
multiplicative reduction over `k⸨q⸩`, and the link to Mathlib's modular forms `E₄`, `E₆`, `Δ`.
The uniformisation `K^×/q^ℤ ≅ E_q(K)`, where the node law would appear as the `q = 0` end, is
not done.

Paths are in the upstream repo unless marked.

1. **Minimality criterion**: `FLT/Mathlib/AlgebraicGeometry/EllipticCurve/Reduction.lean`,
   `WeierstrassCurve.isMinimal_of_valuation_c₄_eq_one` (an integral model whose c₄ is a unit is
   minimal). Our port check found that "minimal model" in `FreyPieces` is a name only. This
   lemma plus our `frey_c4_not_dvd` makes it a theorem. **Worth following**: small, closes a real
   gap, Apache 2.0.
2. **Reduction types**: mathlib `Mathlib/AlgebraicGeometry/EllipticCurve/Reduction.lean`
   (`IsMinimal`, `HasGoodReduction`, `HasMultiplicativeReduction`,
   `HasSplitMultiplicativeReduction`). Restate `frey_semistable` in these terms; the fork's
   `FLT/FreyCurve/Serre/Semistable.lean` (fork) is a worked template.
3. **Split versus nonsplit node**: upstream `FLT/Mathlib/AlgebraicGeometry/EllipticCurve/Reduction.lean`
   (`nodePoly`, `nodePoly_map_splits_iff_isSquare`, `…_of_two_eq_zero`) and
   `FLT/KnownIn1980s/EllipticCurves/QuadraticTwists/SplitMultiplicativeReduction.lean`. Our
   `node_chord_circ` assumes the two tangent slopes m₁, m₂ are rational, which is exactly
   split multiplicative reduction, where the smooth points form the multiplicative group and the
   chord law is the circle operation ∘₁ (called mult-add before 1.0). Their criterion says when that holds at a given prime; nonsplit is the
   twisted (norm-one) version. **Worth following** as a geometric bolt-in: it puts the circle operation at a
   precise place in the FLT machinery. It is known math, so it is placement, not discovery.
4. **Tate curve**: `FLT/KnownIn1980s/EllipticCurves/TateCurve*.lean`, `TateParameter.lean`.
   E(K) ≅ K^×/q^ℤ: the Tate curve is a q-deformation of the multiplicative group, and our node
   law is its q = 0 end. Same verdict as 3: a geometric anchor, not new mathematics. The analytic
   uniformisation there is still sorry.
5. **Normal form and reduction to B3**: `FLT/FreyCurve/FreyPackage.lean` (`FreyPackage`,
   `gcdab_eq_gcdac`, `of_not_FermatLastTheoremFor_p_ge_5`) and `FLT/Basic/Lemmas.lean`
   (`FermatLastTheorem.of_p_ge_5`). In the port, stating `NoFreyCurve ↔ IsEmpty FreyPackage`
   would tie our route to their B3 directly.
6. **For the counterexample-map thread**: `FltRegular/FltRegular.lean` (vendored in the fork;
   `flt_regular`: FLT for regular primes, Kummer's method). Kummer's method is a formalized
   "method that works except at irregular primes" (37, 59, 67, …): a clean example of where an
   approach provably stops.

Toolchain note: our port is on Lean `v4.37.0-rc1`; upstream FLT is on `v4.35.0`. A direct
`lake` dependency on FLT won't build until the versions match, so for now reuse means copying the
needed lemmas (with attribution) rather than importing.

## 4. How to track it

What to watch, best signal first:

1. The `#print axioms` block at the bottom of `FermatsLastTheorem.lean` (both repos). It only
   changes when what FLT rests on changes; `sorryAx` disappearing is the project's own finish line.
2. `FLT/Proof.lean` (upstream): the boss chain. A change means the route or a boss moved.
3. `FLT/Assumptions/` (upstream): a new axiom file means a new assumed theorem.
4. `lean-toolchain` (upstream): when it reaches our Lean version, we can depend on FLT directly.
5. The fork's `cartography/CAMPAIGN-PROGRESS.md` and `cartography/INDEX.md`, and whether any of its
   work starts appearing upstream.
6. Blueprint chapter `ch_freyreduction` (https://imperialcollegelondon.github.io/FLT/blueprint/ch_freyreduction.html)
   for the B4 sub-steps.

`snapshot.sh` in this folder clones both repos and prints one JSON line each (head, date,
toolchain, mathlib rev, axioms line, file and sorry counts). Appending its output to
`snapshots.jsonl` and diffing against the previous line is the whole tracking job. Weekly is
plenty at the current pace; a routine that runs it and reports only when the axioms line,
`Proof.lean`, the assumptions folder or the toolchain changes would be quiet most weeks.

First snapshot (Nov 25, 2026) is in `snapshots.jsonl`.

## Sources read

- Upstream repo at `0198058`: `README.md`, `CONTRIBUTING.md`, `GENERAL.md`, `FermatsLastTheorem.lean`,
  `FLT/Proof.lean`, `FLT/Assumptions/*`, `FLT/FreyCurve/*`, `FLT/EllipticCurve/Torsion.lean`,
  `FLT/GaloisRepresentation/HardlyRamified/Frey.lean`, `FLT/KnownIn1980s/EllipticCurves/*`,
  `FLT/Mathlib/AlgebraicGeometry/EllipticCurve/Reduction.lean`, `blueprint/src/chapter/ch03*.tex`, git log 2026.
- mathlib at `56a0058`: `Mathlib/NumberTheory/FLT/{Basic,Three,Four}.lean`,
  `Mathlib/AlgebraicGeometry/EllipticCurve/Reduction.lean`.
- Fork at `70f8b79e`: `FermatsLastTheorem.lean`, `FLT/Proof.lean`, `FLT/FreyCurve/Serre/Semistable.lean`,
  `cartography/{INDEX,CAMPAIGN-PROGRESS}.md`, its PR #89 page.
- Our files: `lean-library/FusedAlgebra/Topics/Fermat/*`, `mathlib-port-check/*Fermat*`, branch `mathlib-port`.
- Not read: the upstream GitHub issue tracker and project dashboard (automated access blocked).
