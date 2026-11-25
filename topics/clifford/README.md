# Clifford's work through mult-add (Nov 25, 2026)
Page: https://claude.ai/artifact/NU8rLX3jBBeaojcq8g3hMW (slug add9d4ef-0573-4774-a13c-fa98c4cef8e5, version 1795625560-f8f2; source index.html, built by build_page.py from page_template.html + figures/)

Task (Brennan, Nov 25): "let's do a deep dive on Clifford's work and make sure we have an analog for each piece."
Pictures asked for (Brennan: "if I have more pictures to look at, I can help"): one labeled SVG pair per piece.

## Answer
Of 16 pieces of Clifford's work: 12 exact analogs, 2 partial (Clifford–Klein forms, space theory of matter),
1 transfer only (chain of circles), 1 none (theorem on special divisors). Plus Brennan's own question (k as a
multivector): mult-add's dial k is the null (light-cone) element k(1+e)/2 of the split-complex numbers; any
multivector dial K with plain-number arguments gives the landscape law with p = tr K, q = N K; curvature
q - p^2/4 = -(non-scalar part of K)^2. A dictionary to known math, no new prediction.

## Table (piece -> analog -> status)
| # | Clifford piece | Mult-add analog | Status |
|---|---|---|---|
| 01 | Geometric product, Cl(p,q) | landscape law = geometric product of two lifts in the even Clifford algebra of the rate form x^2+pxy+qy^2 (theta = uv) | exact (new framing of groups-lie theta-lift) |
| 02 | Three 2D algebras C, split, dual | circle / line / border by sign of p^2-4q; mult-add = split, theta=(1+e)/2; wall -1 = zero divisor (light cone) | exact |
| 03 | k as multivector (Brennan's question) | p = tr K, q = N K; mult-add's k null; Cl(3,0) general K -> complex (p,q) in span(1,I) | exact, new |
| 04 | Grades / multivectors | a1e1 ✳ ... ✳ anen keeps each subset product on its own blade; ordinary mult-add = scalar shadow | exact, new |
| 05 | Reversion, involutions | conjugation theta -> p - theta; X Xbar = norm; chart of conjugate = landscape inverse; rev(X✳_K Y) = rev Y ✳_{rev K} rev X | exact, new |
| 06 | Rotors, spinors | half-velocities = rotor charts (tan / tanh of half angle); chart forgets R vs -R, lift keeps it | exact (twisted boosts) |
| 07 | Hamilton's quaternions | q=+1 dial boost = quaternion 1+uj; corner = (1-u·v) + i(u×v) | exact (twisted boosts) |
| 08 | Grassmann's wedge | x✳y - y✳x = 2 x∧y; x✳y = x+y+x·y+x∧y | exact, new |
| 09 | Biquaternions q+ωr (1873), ω^2 = +1/0/-1 | R[ω] = theta-lift at (0,-ω^2), set by curvature sign; mult-add theta = elliptic projector (1+ω)/2 | exact, new |
| 10 | Dual numbers | dual args give first jet; border = re-charted dual numbers | exact (known; jets page) |
| 11 | Clifford parallels / Hopf | corner phase arg(X+iY) = half the velocity-triangle area = Hopf holonomy | exact, new link |
| 12 | Clifford torus | circle |t| = 1/sqrt|q|: light wall at q=-1, torus base at q=+1 | exact, small link |
| 13 | Clifford–Klein forms, elliptic space | chart = S^3/±1; tan(pi j/n) finite forms; no 2-D ✳ torus yet (needs two commuting dials + lattice) | partial |
| 14 | Space theory of matter | GR page: time dial gives half light bending, space curvature the other half; waves where river stops | partial (cited) |
| 15 | Bott periodicity | transfers via shift; wall stays a single point only in R, C, H (Cl(0,0),(0,1),(0,2)) | exact (transfer + new wall fact) |
| 16 | Chain of circles (1871) | holds in every chart by Möbius invariance; no ✳ statement | transfer only |
| 17 | Special divisors (1878) | none; would need ✳ on divisors | none |

## Leads (own merits)
- Light-cone dial: worth writing up (unifies landscape, relativity, twist rule); no new law (2D classified).
- Bolting in: landscape laws are rotor compositions in a chart, so GA / dual-quaternion libraries (robotics, graphics, physics) compute them as is.
- Hopf phase = half twist: same geometric phase as Thomas precession / Berry phase; data checks already on twisted-boosts.
- Complex landscape from general K in Cl(3,0): open, low priority (no physical meaning found).
- 2-D Clifford–Klein form from ✳: open construction, adds a picture not a prediction.

## Files (for the repo thread to push)
- README.md (this), checks.py -> checks_output.txt, results.json: 60/60 pass (sympy, numpy). Uses gr-bridges/ga.py (imported by path).
- lean/Clifford.lean: 53 theorems (incl. ext lemmas), core Lean 4.37.0-rc1, no sorry, imports Base and Topics.GroupsLie.
  MERGED into the shared library as FusedAlgebra/Topics/Clifford.lean (after GroupsLie in build.sh, import in FusedAlgebra.lean);
  full clean build: 2,876 theorems, standard axioms only. theorems.tsv and README in lean-library updated.
- figures.py -> figures/*.svg (17 pairs NN-a = Clifford, NN-b = mult-add, plus hero.svg), figures/index.json.
- build_page.py + page_template.html -> index.html (the published page), table.json.
Cited, not repeated: groups-lie, twisted-boosts, general-relativity, gr-bridges (sandwich law S6, task u0e0947bdcc22bc5ce8edb50a6fbb7654).
Source on ω^2: Clifford 1873, "Preliminary sketch of biquaternions", Proc. London Math. Soc. (ω^2 = +1 elliptic, 0 Euclidean, -1 hyperbolic).
