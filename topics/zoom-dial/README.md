# The zoom dial (FHA 1.0 checklist, item 2) — Nov 25, 2026

> **Names:** some words or symbols on this page (mult-add, wall, twist) also mean something else in standard math, physics or computing, or inside this project. [NAMING-COLLISIONS.md](../../NAMING-COLLISIONS.md) says which meaning is in play.

Fills the one gap the fractal-dials page found ("Hidden dials of a half box", math-theory/fractal-dials/): the circle operation ∘₁ (called mult-add before 1.0) had dimension, weights, branching growth and a cousin of wobble, but no zoom (size) dial, so dimension could only be put in by hand.

## Result in plain words
- **The zoom dial plugs in as one new input, and everything after it is a dial the framework already has** (Figure 4). Give each piece a zoom r. Then:
  - a piece shrunk by r carries the weight r^s, and at the true dimension that weight is a twist, k = r^D (weight = twist, Base's circK_scale; the size dial's doubling weight 2^(−β) is zoom 1/2 at dimension β);
  - pieces side by side add their weights: x(s) = Σ rᵢˢ;
  - zooming in again and again is "any number of times", the k = −1 end of the repetition dial: all words of depth n weigh x(s)ⁿ, and all depths together weigh 1/(1 − x);
  - **dimension comes out as an output**: it is where that k = −1 element hits its wall, x(D) = 1. That is Moran's equation (1946).
- **The wobble falls out too.** The complex numbers s with x(s) = 1 are the complex dimensions. Evenly spaced zooms give an evenly spaced column (spacings 4.53, 2.86, 2.27 for sets A, B, E, matching the fractal page's wobble). The uneven set D gives a scattered set and no wobble (Figure 3).
- **What stays outside:** placement and topology (geometry that does not change counts), and the multifractal weights p when they differ from the natural weights r^D (that difference is the separate weights dial, already in the framework).
- **Honest caveat:** the circle operation itself still has no notion of size. The zoom r is a genuinely new input per piece. What this shows is that once a piece has a size, nothing else new is needed.

## Verdict for 1.0
Closed: the fractal dials table moves "zoom" from missing to has-it (one new input, the rest existing dials). The fractal-dials page should say so when it publishes.

## Known vs new (honest)
- All known: Moran 1946 (the equation), Hutchinson 1981 (self-similar sets), Lapidus and van Frankenhuijsen, *Fractal Geometry, Complex Dimensions and Zeta Functions* (2013), whose geometric zeta function for self-similar strings is this word series and whose poles are the complex dimensions.
- New only as packaging: the chain zoom → homothety → k = 0 sum → k = −1 word series → wall, which ties the fractal dials to the repetition result (item 1) and the size dial. Not a discovery.

## Figures
- fig1-zoom-is-weight.png: half boxes A and D drawn to depth 3, each piece labelled with its weight r^D; every depth adds to 1.
- fig2-dimension-at-wall.png: x(s) for six sets; four different half boxes cross the wall at the same D = 1/2; Cantor at 0.631, gasket at 1.585.
- fig3-complex-dimensions.png: the other wall points; evenly spaced columns for A, B, E; scattered for D.
- fig4-how-it-plugs-in.png: the chain, with the existing dial each step uses.

## Files (for the repo push)
- README.md, checks.py, checks_output.txt (35/35; numpy, mpmath), figures.py, fig1–fig4 PNGs.
- lean/ZoomDial.lean: 10 theorems, core Lean 4.37.0-rc1, no sorry, standard axioms only, imports Base only. Library: FusedAlgebra/Topics/ZoomDial.lean, add `FusedAlgebra.Topics.ZoomDial` to build.sh after FractalDials (or after RepetitionDial) and an import in FusedAlgebra.lean. Key theorems: words_weigh, word_series, wall_every_depth, wall_partial, half_box_walls, uneven_depth_three, zoom_is_twist, zoom_compose, halving_weight.
