# Contributing to Fused Hyperoperation Algebra

Thank you for helping. The project's first design principle is to build on published work, and the rules below come from it. The other four are listed after it.

## Design principle: build on published work first

Follow the existing literature wherever possible: its names, notation, conventions and results. Don't invent new ones.
Anything that goes beyond the literature is welcome, but it has to be labeled clearly as new.

## The other design principles (Brennan, Nov 26, 2026)

- **New facts go to the goals.** Try each not-in-print result against quantum gravity and the theory of everything
  first, then the Millennium problems, and record the outcome in `topics/new-facts-levers/`.
- **Look at the dial map first.** For any problem, find the neighbouring work on `topics/dial-map/` before starting.
- **Neighbours include Mathlib.** Search Mathlib's Lean proofs, not just papers, for lemmas to depend on or borrow;
  record useful links in `topics/mathlib-map/`.
- **Watch for new dials** in the affine-composition vein, and add any you find to the dial map.

## What every change needs

1. **Start from published work.** Before working something out from scratch, look for where it is already done, in
   papers, textbooks or mathlib. If it exists, build on it and use its language.
2. **Cite sources.** Give author, year and where it appeared (journal, book, arXiv number, or mathlib name) for every
   result the change rests on.
3. **A "Known vs new" section** in each write-up, next to its Assumptions section:
   - *Known*: what is already in the literature, with citations.
   - *New*: what this change adds, labeled plainly, with how sure we are ("proved in Lean", "checked numerically",
     "not found stated in the literature", and so on). If nothing is new, say so. A clean re-derivation of known math
     is still useful when it is labeled that way.
4. **An Assumptions section** that says what is assumed: the setting (numbers, rings, fields), sign conventions, and
   anything taken on trust.
5. **Follow existing names.** Use the standard name and symbol for a thing when one exists. For example, the circle
   operation is written ∘ₖ on the page, always with its subscript, because a bare ∘ is function composition. If a
   name you use means something else elsewhere, add a row to [NAMING-COLLISIONS.md](NAMING-COLLISIONS.md).
6. **Lean proofs: no `sorry`, standard axioms only** (propext, Classical.choice, Quot.sound). `lean-library/build.sh`
   checks this for the core library through `Check.lean`, and the CI workflow runs it on every push and pull request,
   along with every topic's `checks.py`. The Mathlib build workflow checks the mathlib edition on GitHub.
7. **A labeled picture for each result.** Label the axes, the curves and the key points, so the result can be seen and
   not only read. Keep the script that makes the figure next to it.

## Where things go

- One folder per topic under `topics/`, with the write-up, scripts, saved script output and figures together.
- New Lean files go under `lean-library/FusedAlgebra/` and in the library list.
- Update INDEX.md when you add a topic or page. Add any new naming collision to NAMING-COLLISIONS.md.
- Don't change VERSION in a pull request. Releases are Brennan's call.

## Forking

See [FORKING.md](FORKING.md) for what forking or reusing this work involves: licenses, credit, naming, citation,
and how to stay compatible.

## Before you push

CI ([`.github/workflows/ci.yml`](.github/workflows/ci.yml)) has to pass. After changing the Lean library, copy
`lean-library/build/theorems.tsv` over `lean-library/theorems.tsv`, then run `python3 tools/theorems_md.py` and
`python3 tools/assumptions.py` and commit what they write. `python3 tools/run_checks.py` runs every topic's checks
locally, and `python3 tools/assumption_sections.py` lists write-ups still missing an Assumptions section.

## Pull requests

The pull request template asks for the items above as a checklist. A pull request that is missing one can still be
opened. Say what is missing and why.
