# Fused Hyperoperation Algebra

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

## Lean

All the Lean proofs are gathered in [`lean-library/`](lean-library/README.md): one shared library built on a single base of mult-add laws. As of Nov 25, 2026 it holds 1,750 theorems in 62 modules, every one compiling with **core Lean 4 only** (no mathlib), none using `sorry`, and each using only Lean's three standard axioms. `lean-library/THEOREMS.md` lists them all; `./build.sh` rebuilds and re-checks them.

Some topic folders also keep the original Lean file they were first written in. The one file that imports mathlib, `topics/tesseract/TesseractCoords.lean`, writes mult-add as `✳` (mathlib already uses `⊞` for biproducts) and compiles cleanly against mathlib v4.37.0-rc1 (checked on Brennan's Mac Studio, Nov 25, 2026). A port of the whole library to mathlib is under way on its own branch.

## Where to look

- [INDEX.md](INDEX.md) lists every topic with its folder and its published page.
- [IN-PROGRESS.md](IN-PROGRESS.md) lists folders that were still being written when this repo was made.
- [DATA.md](DATA.md) lists the large files left out and where they came from.

Each folder under `topics/` usually holds the check scripts (`checks.py`), their saved output, the figure scripts and figures, any Lean proofs, and the page source (`index.html` or `page.html`).
