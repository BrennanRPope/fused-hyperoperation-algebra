# Size weighting as a dial, and what mult-add leaves over (Nov 25, 2026)

Task (Brennan, voice call and timeline): is the 54.5% (plain ½ on the naturals by size, not 50%) a hidden dial we missed, and, across the known sets, how much does mult-add with its dials explain and what is left over? Asked during the work for clear pictures.

Page: published as an artifact (see PAGES-INDEX.md, "size dial and residual"); source index.html here (built by build_page.py from page_template.html and the three SVGs).

## Result in plain words
- The 54.5% is exact: Li2(-1/2)/Li2(-1) = -12 Li2(-1/2)/pi^2 = 0.545206.
- Counting by size adds one number, the dimension beta of the size measure (how fast the count of elements grows with size). Plain ½ keeps Li_{beta+1}(-1/2)/Li_{beta+1}(-1): 58.50% (sharp count, beta=0), 54.52% (naturals by size, beta=1), 52.42% (plane), 51.29% (3-D), exactly 50% for primes, Dirichlet weights, log density (beta=infinity).
- It is not a new kind of dial: a size weight w is the twist k = w on that element (wAdd_scale). For small dials the spread acts as one twist exactly, k = 2^-beta (q = 1 + 2^-beta). Over all dials one global twist comes within 0.06 points (k = 0.441 at beta=1) but not exactly. Weighting does not commute with merging (gap t s w (w-1)) and breaks complement on infinite sets.
- What mult-add (subsets, each element at most once) leaves over is repetition, and it is always the same set again at double size: share 2^-beta of the log-total (1/2 for partitions, 1/4 plane, 1/8 in 3-D = the Stefan-Boltzmann 7/8). For the number-theory sets (beta=infinity) it touches no growth exponent, only the constant: 6/pi^2 (naturals), prod(1-2/(p(p+1))) = 0.4717 (rationals by height), 1/zeta_K(2) (number fields). The dial ladder t, t^2, t^4, ... (Fermi-Dirac atoms) closes it exactly.
- Second kind of leftover: the class group. Z[sqrt-5] numbers are exactly 1/2 of the ideals (measured 0.4999); sign dials repair it (even_half).
- Algebraic numbers under multiplication: no atoms (every one is a square), so mult-add cannot start. As a set by height they behave like the primes.

## Verdicts
- Size weighting: real, exact, known (Fermi-Dirac functions). Keep as a geometric reading of q: q = 1 + 2^-beta.
- Missing piece: none hidden. Leftovers are named (repetition; class group) and both close with things mult-add can express.
- Lead: a continuous "how much repetition" dial (Haldane exclusion statistics) with mult-add at one end. Known physics; useful frame.

## Files (for the repo push)
- README.md, index.html, page_template.html, build_page.py
- checks.py, checks_output.txt (mpmath, numpy)
- figures.py, fig1-dimension-dial.svg, fig2-explained-residual.svg, fig3-doubling-ladder.svg
- lean/SizeDialResidual.lean: 18 theorems, core Lean 4.37.0-rc1, no sorry. Compiled against a copy of the library (add `FusedAlgebra.Topics.SizeDialResidual` after DialSurvey in build.sh and an import line in FusedAlgebra.lean; library then 2819 theorems, 112 modules, Check passes).
