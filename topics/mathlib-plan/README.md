# Mathlib port: what it takes to merge (FHA 1.0 checklist, item 3) — Nov 25, 2026

> **Names:** some words or symbols on this page (⊞, mult-add, fused, Tate q and E₄, E₆) also mean something else in standard math, physics or computing, or inside this project. [NAMING-COLLISIONS.md](../../NAMING-COLLISIONS.md) says which meaning is in play.

> **Update (Nov 27): a Mathlib-candidates library.** `lean-library-mathlib/MathlibCandidates/` holds files written in Mathlib's own style to offer to Mathlib later, starting with the Tate curve (asked for by Brennan). CI builds it as a default target and `Check.lean` checks its axioms. Still no pull request to Mathlib; the Tate curve overlaps Buzzard's FLT project and should be merged with it first.
>
> **Update (Nov 25, ~22:30 UTC): merged into main.** Brennan asked for everything on main, so `mathlib-port` was merged and main's CI now builds `lean-library-mathlib/` (the port plus the proof-complexity study) on every push. Still no pull request to Mathlib.
>
> **Earlier update (Nov 25, ~22:00 UTC): the build check is in place and green.** The Mathlib build workflow (`.github/workflows/mathlib-build.yml`) now runs on the `mathlib-port` branch too (commit 86b5ae3 on that branch). Its first push run built the port with mathlib's cache on GitHub and checked it: 1,970 theorems, standard axioms only, no `sorry` ([run 3](https://github.com/BrennanRPope/fused-hyperoperation-algebra/actions/runs/45708349981)). The branch is still not merged, and no pull request to Mathlib is planned until the theory is more complete (Brennan). Step 2 below is done; the "can't be built in this cloud" point still holds for this workspace, but GitHub builds it.

## Answer
- **Merging the branch is easy; trusting it needs one build.** Branch `mathlib-port` (head 99b18c1) only adds the folder `lean-library-mathlib/` plus one INDEX.md line, so it merges with at most a one-line INDEX.md conflict. The branch notes say `lake build`, mathlib's linter and the theorem check all passed (62 modules, 1,970 theorems, standard axioms only). The last two commits (restoring the four helper lemmas the review found, then updating the tables) came after that, and I found no record of a build since.
- **It can't be built in this cloud.** The proxy here refuses mathlib's prebuilt-file cache, and building mathlib from source is hours of compute. Two ways to build: Brennan's Mac Studio (has mathlib), or a GitHub Actions check on the public repo that runs `lake exe cache get && lake build` (GitHub's runners can reach the cache; free for public repos). The check would make the mathlib edition self-verifying without the Mac.
- **About half the library is unported** (fig1-port-coverage.png): 62 of 119 modules ported. The 57 modules added since the port started (Clifford, groups and Lie, electromagnetism, GR, the dial pages, biology, repetition and zoom dials, …) hold about 1,280 theorems. None of the 62 ported originals has changed since, so the existing port is not stale.

## Order: names before porting
Renaming after porting means renaming twice (core and mathlib). So checklist item 4 (names and notation) should come before porting the newer modules. Two naming facts the port already surfaced:
- Mathlib already has the circle operation (called mult-add before 1.0): it is the circle operation of ring theory, the monoid on `PreQuasiregular R`; "has a way back" is mathlib's `IsQuasiregular`.
- Notation: the core library's comments write the circle operation as ⊞, the port writes ∘₁ because mathlib uses ⊞ for biproducts. One symbol should be frozen for 1.0.

## Is the port a 1.0 blocker?
Recommendation: no. The port restates the same definitions in mathlib's style and adds a folder. It changes no definition, notation or theorem name of the core library, so under the versioning rule it can land in 1.x. What 1.0 needs is the core library (core Lean, no sorry, builds) and frozen names.

## Plan
1. Item 4: freeze names and notation.
2. Add a GitHub build check for `lean-library-mathlib/` (needs Brennan's yes; it adds a workflow file).
3. Once the check is green, merge the branch as the "mathlib edition, 62 of 119 modules", with the coverage stated in its README.
4. Port the 57 newer modules in batches after the names freeze (similar size to the first port, which was 1,750 theorems), each batch checked by the same build.
