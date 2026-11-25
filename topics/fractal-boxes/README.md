# Fractal boxes: does a box of fractional dimension give a fractional power set? (Nov 25, 2026)

Task: Brennan's yes on the voice call to "Should I look into whether fractal boxes give a consistent fractional power set?" Background: the circle-operation box with sides 1+a, 1+b, ... has 2^n corner pieces, one per subset ("the tesseract is a picture of a power set"). Option tested: a box of fractional dimension d, with 2^d corner pieces. Picture-led write-up at his request ("if I have more pictures to look at, I can help").

Page: https://claude.ai/artifact/NbFUcWVphAGGNvhug8XvwS (source index.html here; build_page.py inlines the SVG plates).

## Result in plain words
A box of half a dimension does have √2 corner pieces in an exact sense, and two of them join into a stick. But that is the joining half we already had (dial √2 − 1). It can't make the iterating half: the shape isn't unique, a shape's own power set needs a half-exponential, and at infinity every fractal has the same number of points (P(ℕ)), so it can't reach the choice-free open lead. Verdict: no new lead; it closes the fractal route with clear reasons.

## Findings
1. Corner layers of a box of dimension d with sides 1+a are Newton's C(d,k): (1+a)^d = Σ C(d,k) a^k = 1 + a^(∘₁d). Half box: 1, 1/2, −1/8, 1/16, −5/128, 7/256, ... (negative from layer 2, alternating); total √2.
2. Joining multiplies (convolves layers); two half boxes give 1, 1, 0, 0, ... exactly (a stick). Three third-boxes give a stick. A half box with sides 2 = one stick of side √2: the per-element dial √2 − 1 again. 2^d is forced by continuity.
3. As a shape: C_1/4 (keep two outer quarters) has dimension 1/2; its square is the four-corner (Garnett) set, dimension 1 but dust, Favard length 0 (known). Many other dimension-1/2 fractals (three pieces at 1/9; uneven 0.09/0.49); the count at scale 1/2^k divided by 2^(k/2) wobbles (period one level) and differs by shape. Not canonical as a shape; only the number 2^(1/2) is.
4. Iterating: a shape's power set is its hyperspace (closed pieces, Hausdorff distance). At a scale where the shape shows N separated points, the hyperspace shows 2^N − 1 (subsets of an ε-net are an ε-net; subsets of a separated set stay separated). A half step needs h(h(N)) ≍ 2^N: a half-exponential. These exist (Kneser 1950, base e), are not unique, and no exp-log function is one. Products, unions and Lipschitz maps change log-counts at most linearly, and twice-linear never reaches 2^N (Lean no_linear_half_pow2). Growth-rate counting does remove the n = 0 finite wall.
5. Infinite: every Cantor set of any dimension is P(ℕ) as a set (left/right code one to one, Lean code_inj). Closed sets are countable or contain a Cantor set (Cantor–Bendixson, ZF), so no fractal has size strictly between ℕ and P(ℕ). The choice-free open lead (fractional-powerset/) cannot be reached by geometric sets.
6. Other candidates: magnitude (Leinster) has magnitude dimension = Minkowski dimension for compact subsets of R^n (Meckes 2015), so it adds nothing. Groupoid with 1/2 element (point with Z/2 symmetry): maps into {0,1} number 2, not √2. Species: P = E·E under labelled product, so the half there is E (count 1 per set), not √2^n: a different joining gives a different half.

## Consistency tests
Joining yes; iterating no; canonical number yes, canonical shape no; fits the circle operation ∘₁ (called mult-add before 1.0) yes (fractional circle power); merging/complement same as the dials page; reaches the choice-free lead no.

## Files
- checks.py / checks_output.txt: 27 of 27 pass (pure Python: fractions, math, itertools).
- figures.py, figures/p1..p6 *.svg: the six plates (white ground). build_page.py assembles index.html.
- lean/FractalBoxes.lean: 15 theorems, core Lean 4, no imports, no sorry. Compiled in the shared library as Topics/FractalBoxes (Lean 4.37.0-rc1; library now 115 modules, 2,915 theorems, Check passes). build.sh, FusedAlgebra.lean, README.md and theorems.tsv in lean-library updated.

## For the repo thread (push to BrennanRPope/fused-hyperoperation-algebra)
- whole folder math-theory/fractal-boxes/ (README.md, checks.py, checks_output.txt, figures.py, figures/*.svg, build_page.py, index.html, lean/FractalBoxes.lean)
- lean-library: FusedAlgebra/Topics/FractalBoxes.lean, build.sh, FusedAlgebra.lean, README.md, theorems.tsv
- PAGES-INDEX.md (new line "fractal boxes")
- README catalog: no line added. Everything here is known math; it closes a route but isn't a new result.

## Sources
Kneser, J. reine angew. Math. 187 (1950); Wikipedia "Half-exponential function"; Nazarov–Peres–Volberg / Łaba–Zhai, arXiv 0902.0964 (four-corner Favard length); McClure, Real Analysis Exchange 21 (1995/96) (hyperspace entropy dimensions); Meckes, Potential Analysis 42 (2015) (magnitude dimension).
