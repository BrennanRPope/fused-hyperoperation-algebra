# Road to 1.0

Last updated: 2026-11-25, 18:30 UTC. Current version: 0.1.0 (see VERSION).

## What 1.0 means

1.0 is the point where Brennan decides the framework is settled. From then on the definitions, the notation and the
theorem names are frozen. Changing any of them means 2.0. New results, new pages, new Lean theorems and new dials that
leave the existing ones untouched can all come later as 1.x releases.

So the finish line is not "every question answered". It is: **the base is settled, and every known gap is either
closed or written down as a named extension that would not change the base.**

The long-range aim, a framework that can explain all of known physics and all of known math, sits beyond 1.0. The
"After 1.0" list below tracks it, one named gap at a time.

## The checklist

| # | Item | Status | Where |
|---|---|---|---|
| 1 | Repetition dial: does repetition fall out of the framework? | **Done.** Repetition inside one element is k itself, with 1/k cells per element and a smooth sweep through k = 0. 22 Lean theorems, 29 checks. Crowding *between* elements (Haldane's g) stays outside; whether that blocks 1.0 is Brennan's call. | `topics/repetition-dial/` |
| 2 | Zoom dial for the fractal dials | **Done.** A zoom r is a weight r^d, and the dimension is an output, where the zoom sum meets the wall. 10 Lean theorems, 35 checks. | `topics/zoom-dial/` |
| 3 | Mathlib port | **Sized; not a 1.0 blocker.** The branch `mathlib-port` merges cleanly but needs one build. A GitHub workflow now builds it with mathlib's cache (commit efa8013). 57 newer modules are unported, about 1,280 theorems. | `mathlib-port` branch, `lean-library-mathlib/PORT-NOTES` |
| 4 | Names and notation before the freeze | **Reviewed; waiting on Brennan.** Four fixes proposed: (a) ✳ as the only mult-add symbol; (b) k is defined on the quantity composed, and each page names that quantity (weights compose by ✳ₖ, probabilities by ✳₋ₖ); (c) retire q = k + 1; (d) call the plane's second dial d, after Coleman–McGuinness 1991. A credit glossary (NAMES.md) is drafted. | `topics/names-review/` |
| 5 | This roadmap, kept current | **Done when it is in the repo.** It is updated after each item. | `ROADMAP.md` |
| 6 | **Naming collisions, visible to readers (required for 1.0)** | **Drafted.** NAMING-COLLISIONS.md lists every name or symbol that means something else in standard math, physics, computing or mathlib: ⊞, ⊕, ⊛, q (twice), k's sign, the proposed d, mult-add, Fused, wall, twist. Each row says what ours means, what theirs means, and since which version, starting with the released v0.1.0. Two clashes are inside the project itself: ⊞ and ⊛ each mean two different things in the Lean library. Still to do: link the page from the README and from each write-up that uses one of these names. | `NAMING-COLLISIONS.md` |

## Left before 1.0

1. Publish NAMING-COLLISIONS.md and link it from the README and each affected write-up (item 6, required whatever
   is decided on renames).
2. Brennan decides on the four naming fixes (item 4). Adopting d for the plane needs the zoom dial's dimension
   written another way, since it also uses d.
3. If adopted: apply them across pages and the Lean library in one pass, keep the old names as dated rows on
   NAMING-COLLISIONS.md, add NAMES.md, and check that the build still
   passes. The rename of q to d also touches the relativity pages. The released v0.1.0 notes keep "q = −1/c²", and
   the collisions page records it.
4. Brennan decides whether crowding between elements (Haldane's g) must be inside the framework before 1.0, or can be
   an extension dial in 1.x.
5. Run the mathlib workflow on `mathlib-port` and merge it if it passes. Porting the remaining 57 modules can happen
   after 1.0.
6. Brennan sets the version to 1.0.0 and writes the release notes.

## After 1.0 (named extensions, none changes the base)

- Crowding between elements: Haldane–Wu exclusion statistics as a second dial, unless it moves into the list above.
- Porting the remaining modules to mathlib.
- Topic coverage toward "all of known physics and math": each topic page states what the framework explains there,
  what it only restates, and what it cannot reach yet.

## How to update this file

When an item changes, edit its row and the date at the top. Keep each status short: what was found, and what it waits
on. Do not change VERSION here; that is Brennan's call.
