# Naming collisions

Some names and symbols in Fused Hyperoperation Algebra mean something else in standard mathematics, in physics, in
computing, or in Lean's mathlib library. This page lists every one we know of, so that a reader who meets one of ours
(or one of theirs) can tell which meaning is in play.

Version 0.1.0 is already released, and it uses all of the names below. Nothing here renames anything. Renames
proposed for 1.0 are in [the names review](topics/names-review/README.md). Once 1.0 is out, this page records which meaning each version
uses.

Last checked: 2026-11-25, against the repo at commit efa8013 and a fresh copy of mathlib4.

## Symbols

| Symbol | What it means here | What it means elsewhere | Where it appears here | Since |
|---|---|---|---|---|
| ⊞ | Mult-add, a ⊞ b = a + b + ab. It is the scoped Lean notation in `Base.lean` and appears about 1,630 times in the repo. **It has two other meanings inside this project too:** the rung-1 shifted add a + b + 1 in one `Rungs.lean` comment, and the product `MARing.madd` in `MultAddRing.lean`. | Mathlib: the biproduct X ⊞ Y of two objects in a category. Mathlib chose ⊞ because ⊕ was already taken. | `Base.lean`, most write-ups, `Rungs.lean`, `MultAddRing.lean` | v0.1.0 |
| ⊕ | Mult-add, in the README header and in `Rungs.lean`. Also the whole two-dial law (a + b + k·ab)/(1 − q·ab) on the landscape pages. | Lean core: the sum of two types, α ⊕ β. In mathematics generally: direct sum. In computing: exclusive or. | README, `Rungs.lean`, landscape and proof-engine pages | v0.1.0 |
| ⊛ | Two meanings here: the rung-3 operation in `Rungs.lean` comments, and mult-add itself (`opMA`) in `FunctionalAnalysis.lean`. | Mathlib: Day convolution of functors, and applying a stream of functions. In analysis and signal processing: convolution. | `Rungs.lean`, `FunctionalAnalysis.lean`, about 300 uses | v0.1.0 |
| ✳ | Mult-add. It is declared locally in about 80 Lean files and used in most write-ups (about 1,970 uses). | No use found in mathlib, but mathlib's character linter does not allow it, and Lean's editor has no shortcut for it. | Most of the repo | v0.1.0 |

## Letters

| Letter | What it means here | What it means elsewhere | Where it appears here | Since |
|---|---|---|---|---|
| q, as q = k + 1 | A second name for the dial k, shifted by 1. | Tsallis statistics uses q for the same family with **q_Tsallis = 1 − k**, so our q is 2 − q_Tsallis. Our q = 2 is their q = 0, for example. Also, q in q-analogs (q-binomials, quantum groups) is a base, an unrelated quantity. | Dial survey, information theory pages | v0.1.0 |
| q, the plane's second dial | The d-like dial in the law (a + b + k·ab)/(1 − q·ab). In relativity q = −1/c², as the v0.1.0 notes say. | Coleman and McGuinness (1991) call this dial d, with the same sign. **Inside the project it clashes with the other q** (the line above). | Landscape, relativity and electromagnetism pages, v0.1.0 release notes | v0.1.0 |
| k | The dial in a + b + k·ab, defined on whatever quantity is being composed. Pages that compose **weights** have k = +1 for "at most once". Pages that compose **probabilities** (or use the factor (1 − kz)^(−1/k)) have k = −1 for "at most once". Both are correct, because weights composing by k means probabilities compose by −k (proved in `Conventions.lean`). | Shirai and Takahashi's α for determinantal processes matches our probability side (α = −1 is determinantal, which is fermions). Tsallis's q matches 1 − k. In physics k often means Boltzmann's constant, and in algebra it often means a field. Ours means neither. | Every page with a dial | v0.1.0 |
| d (proposed for 1.0, not adopted) | The proposed new name for the plane's second dial, matching Coleman and McGuinness. | **It would clash with d for dimension in the zoom dial** (Moran's equation x(d) = 1). The zoom dial was added after v0.1.0. This has to be settled before d is adopted. | Proposed only; zoom dial uses d for dimension | not released |

## Words

| Word | What it means here | What it means elsewhere | Since |
|---|---|---|---|
| mult-add | a + b + ab. The standard algebra names are the circle operation or quasi-regular product (Perlis 1942, Jacobson 1945). Mathlib has it as `PreQuasiregular`, written y + x + xy. | In computing, "multiply-add" and "multiply-accumulate" mean a·b + c, a different operation. | v0.1.0 |
| Fused (in the project name and in "fused exponent", "fused tetration") | An ordinary operation conjugated by x ↦ 1 + x. For example, fused exponent = (a + 1)^(b + 1) − 1. | In computing, "fused multiply-add" (FMA) is the rounding-once a·b + c instruction. In physics, "fusion" means fusion rules and fusion categories. | v0.1.0 |
| wall | The point a = −1/k where 1 + k·a = 0, so there is no way back (inverse). | "Walls" and "wall-crossing" in geometry and string theory (stability conditions) are unrelated. | v0.1.0 |
| twist | Weight is twist: a ✳ₖ b rescales to mult-add (`wAdd_scale`). Twisted boosts: twist = q × area. | Drinfeld twists, twisted group algebras and cocycle twists in algebra. The twisted-boosts use may be close to a cocycle twist; this is not checked. | v0.1.0 |

## Not collisions, but private names for known things

Words such as rung (hyperoperation level), dial (parameter) and way back (inverse) name standard ideas in our own
words. They are mapped to the standard names and first sources in the credit glossary (NAMES.md), not here.

## How to use this page

- A write-up that uses one of the names above should link here the first time it uses it.
- If 1.0 adopts a rename, this page keeps both rows, marked with the version each one applies to, so 0.1.0 readers can
  still follow.
- Add a row whenever a new collision is found. Note the date and where it was found.
