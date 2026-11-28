# Road to 1.0

Last updated: 2026-11-28, 08:35 UTC. Current version: 0.3.0, released Nov 28 (see VERSION and CHANGELOG.md).

## What 1.0 means

1.0 is the point where Brennan decides the framework is settled. From then on the definitions, the notation and the
theorem names are frozen. Changing any of them means 2.0. New results, new pages, new Lean theorems and new dials that
leave the existing ones untouched can all come later as 1.x releases.

So the finish line is not "every question answered". It is: **the base is settled, and every known gap is either
closed or written down as a named extension that would not change the base.**

The long-range aim, a framework that can explain all of known physics and all of known math, sits beyond 1.0. The
"After 1.0" list below tracks it, one named gap at a time.

## Design principles (the first is a 1.0 criterion)

Build on existing published work first, follow the literature's names, notation and conventions, and label anything
beyond it as new (Brennan, Nov 25, 2026). 1.0 requires that every write-up has an Assumptions section and a
Known vs new section with sources, and that CONTRIBUTING.md, the pull request template and FORKING.md say so.

Four more principles guide the work (Brennan, Nov 26, 2026): try every not-in-print fact against quantum gravity and
the theory of everything first, then the Millennium problems; look at the dial map first for neighbouring work; search
Mathlib's Lean proofs as well as the literature, as functions to depend on or borrow; and watch for new dials in the
affine-composition vein. See the README's Design principles.

## The checklist

| # | Item | Status | Where |
|---|---|---|---|
| 1 | Repetition dial: does repetition fall out of the framework? | **Done.** Repetition inside one element is k itself, with 1/k cells per element and a smooth sweep through k = 0. 22 Lean theorems, 29 checks. Crowding *between* elements (Haldane's g) stays outside; whether that blocks 1.0 is Brennan's call. | `topics/repetition-dial/` |
| 2 | Zoom dial for the fractal dials | **Done.** A zoom r is a weight r^s, and the dimension D is an output, where the zoom sum meets the wall. 10 Lean theorems, 35 checks. | `topics/zoom-dial/` |
| 3 | Mathlib port | **Built in CI; not a 1.0 blocker.** A Mathlib version of the library (`lean-library-mathlib/`) is on main and built by `mathlib-build.yml`, with a Mathlib map of 16 modules and 71 theorems. Porting the rest of the library's modules can come after 1.0. | `lean-library-mathlib/PORT-NOTES.md` |
| 4 | Names and notation before the freeze | **Done; merged in v0.2.0.** Mult-add is the circle operation ∘ₖ (Lean ◯), the plane's second dial is d, the zoom dial uses s and D, and the old Lean names still work as deprecated aliases. On Nov 26 the dials got affine-geometry names (k homothety, d curvature, p the ℓp exponent, C the cubic form), with readings such as Δ, D and the holonomy ω kept as outputs only (`topics/affine-names/`). Still open: the shorthand picked Nov 25 (ℓₖ, φₖ, a +⁽ˢ⁾ b, fuse, σ_c) is on a few pages only, the credit glossary NAMES.md isn't in the repo, and figures with baked-in labels still show old names. | `topics/names-review/`, `topics/affine-names/` |
| 5 | This roadmap, kept current | **In the repo as ROADMAP.md.** Updated after each item. | `ROADMAP.md` |
| 6 | **Naming collisions, visible to readers (required for 1.0)** | **Done.** NAMING-COLLISIONS.md is at the repo root with the old-to-new table, and the README links it twice. Every topic write-up that uses a colliding name links it near its title (added Nov 28). Two frozen folders, `topics/bridges-study/` and `topics/why-i/bridge-paper/`, carry no Names line on purpose, because their SHA256SUMS freeze the README byte for byte. | `NAMING-COLLISIONS.md` |
| 7 | Design principle, contributing guide, PR template, forking guide | **In the repo.** CONTRIBUTING.md, FORKING.md, NOTICE and ASSUMPTIONS.md are on main. Still to do: 27 of the 172 topic write-ups have a Known vs new section, so 145 need one. | `CONTRIBUTING.md`, `FORKING.md`, `NOTICE` |
| 8 | **The theory's name, fully decided (required for 1.0)** | **Done. Affine Composition Theory, built on Fused Hyperoperation Algebra as its own layer.** Picked Nov 26; Brennan reconsidered on Nov 27 and settled it ("Keep Affine Composition"). The algebra holds the operation, its laws, the rungs and the Lean library; the theory holds the shapes, dials, readings, order dial and physics (`topics/theory-name/LAYERS.md`). The README title already uses both names; the repo, CITATION.cff and .zenodo.json keep the algebra's name, since the repo is named after the algebra layer. Earlier interim name (Nov 25): Hypercube Composition Theory, dropped because the shapes now run from round to pointy. | `topics/theory-name/`, `topics/names-review/THEORY-NAME.md` |

## Left before 1.0

1. Add a Known vs new section to the 145 write-ups that lack one (item 7).
2. Finish the naming pass: apply the Nov 25 shorthand everywhere or drop it, add NAMES.md, and redraw figures that
   still show old labels (item 4).
3. Get CI green. The Lean build is the core; Brennan asked (Nov 28) to keep the CI build to unit-test-style scripts,
   since the full topic checks ran past the 120-minute limit.
4. Brennan decides whether crowding between elements (Haldane's g) must be inside the framework before 1.0, or can be
   an extension dial in 1.x.
5. Brennan sets the version to 1.0.0 and writes the release notes.

## After 1.0 (named extensions, none changes the base)

- Crowding between elements: Haldane–Wu exclusion statistics as a second dial, unless it moves into the list above.
- Porting the remaining modules to mathlib.
- Topic coverage toward "all of known physics and math": each topic page states what the framework explains there,
  what it only restates, and what it cannot reach yet.

## How to update this file

When an item changes, edit its row and the date at the top. Keep each status short: what was found, and what it waits
on. Do not change VERSION here; that is Brennan's call.
