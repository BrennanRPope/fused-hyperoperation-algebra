# Graph theory through the circle operation (Nov 25, 2026)

> **Names:** some words or symbols on this page (mult-add, wall) also mean something else in standard math, physics or computing, or inside this project. [NAMING-COLLISIONS.md](../../NAMING-COLLISIONS.md) says which meaning is in play.

Files
- index.html: the page (built by build.py from template.html; self-contained).
- checks.py -> checks_output.txt: 54 checks, all pass (brute-force multivariate Tutte / Potts sums, networkx, numpy, sympy).
- lean/GraphTheory.lean: 38 theorems, compiles clean (imports Topics/Statistics); lean/README.md.

Findings
- Two-terminal reliability of series-parallel networks is exactly the dial pair (x, *_{-1}); failure side swaps them (De Morgan); in odds-against u = 1/p - 1 series is the plain circle operation (called mult-add before 1.0). Wheatstone bridge needs pivoting (deletion-contraction). Known (Moore-Shannon 1956).
- Potts / random-cluster with Q colours: weights v = e^J - 1 = E(J) the circle operation in parallel; in t = v/(Q+v) series multiplies and parallel is the landscape law F_{Q-2,1-Q}, on the "or" line k + q = -1 (tangent to the border at Q = 0). Q = 1 reliability, Q = 2 velocity addition (tanh), walls t = 1 and v = -1 (proper colouring), and it is AMH union with theta = 1 - Q. Rules known (Sokal 2005); landscape placement and AMH reading ours.
- Diamond hierarchical lattice = the dial pair iterated: percolation fixed point (sqrt5-1)/2, slope 4 - 4p*, exponent 1.635; Ising t^3+t^2+t = 1; Q = 3 exactly 1/2. Known (Berker-Ostlund 1979).
- Strong product = Kronecker the circle operation A(x)I + I(x)B + A(x)B; Cartesian = k = 0, tensor = corner; degrees and eigenvalues fuse; 2m+n multiplies; OR product = k = -1 with J. Known (Sabidussi; Hammack-Imrich-Klavzar).
- Independence polynomial: union circle operations I - 1, join adds; the Lee-Yang dial gas at k = 0 is I(G). Known (Gutman-Harary 1983).
- Ihara zeta: Bass formula verified; no structure beyond the Euler product. Dropped.

Verdict: lead worth following: the diamond lattice on the landscape (zero sets / Julia sets of t -> F_{Q-2,1-Q}(t^2,t^2) as Q slides; test the Beraha values from the Lee-Yang page).
Links used: statistics FwYdFUqfFh4DEAjHMqNEsv, relativity S2RVDVVQWR4UAQVtmDQQMq, Lee-Yang UQQzmkzekxe2hL4m8F8NQL.
