# Combinatorics through the circle operation (Nov 25, 2026)

Files
- index.html: the page (built by build.py from template.html + figures.json; self-contained).
- checks.py -> checks_output.txt: 59 checks, all pass (sympy, exact fractions, networkx brute force).
- figures.py -> figures.json: flow polynomials P_n(t) (n <= 10) and the exact sign grid of P_n(t), n <= 16, t in steps of 1/40.
- lean/Combinatorics.lean: 35 theorems, compiles clean; lean/README.md.

Findings
- The cube prod(1 + x_i) - 1 is the subset generating function; the dial paints each chosen element in k colours, (1+k)^n = q^n; k = -1 is inclusion-exclusion (derangements), k = -2 parity. Known (Euler, de Montmort).
- n! - 1 = 1 * 2 * ... * (n-1) under the circle operation (called mult-add before 1.0); the cube with axis weights 1..n-1 has Stirling cycle numbers as layers; at k = -1 the axis of weight 1 sits on the wall, so even and odd permutations balance. Known (Stirling); the wall reading is ours.
- Exponential formula = rung flow E at exactly t = 1 (E(connected graphs) = all graphs; log(1+x) = E^-1). The dial version (e^{kC}-1)/k counts blocks: Stirling S(n,j). Known (Riddell-Uhlenbeck, Joyal).
- Integer rungs of E on x count nested groupings: t = 1 sets, t = 2 Bell numbers, t = 3 partitions of partitions (A000258), t = -1 signed cyclic orders (Moebius function of the partition lattice). Known (Jabotinsky, Comtet).
- Fractional rungs are not counts: P_3 < 0 on (0, 1/3), the half rung has a 0 x^4 and negative x^6 coefficient, t = 3/2 turns negative at n = 16 (checked exactly). Consistent with Baker 1958 (divergence).

Verdict: no lead worth following inside combinatorics (standard generating-function toolkit seen through the lift). One firm fact for the rung work: whole rungs mean depth of nesting; fractional rungs are only interpolation.
Links used: statistics FwYdFUqfFh4DEAjHMqNEsv, q 9AtexUZrzvsPfzSA7g91dL.
