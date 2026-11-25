# Hidden dials of a half box

**Result.** Many fractals share dimension ½. Six known dials tell them apart: zoom, placement,
wobble, weights, topology and branching growth. Mult-add already has the two that matter for
counting (dimension = fractional mult-add power; weights = per-element dial t = p/(1−p)), has the
growth dial separately (rung flow at t = ½), and has a cousin of wobble (complex dimensions act like
a phase on the zoom, resembling twist or a Clifford rotor). It is missing one clean dial: a zoom
(size) dial, which would make dimension an output, d = log(pieces)/log(1/zoom). Placement and
topology are geometry that does not change counts (topology is silent below d = 1: every dust is
the Cantor set). Binary branching + equal zooms + fair weights pin exactly one half box, two pieces
at 1/4; joined with itself it is the binary tree, and gluing dyadic ends gives the stick. That pins
the joining half (the dial √2 − 1), not an iterating half. Shape dials keep counts c·bⁿ and joining
keeps that form; the power set's 2^(2ⁿ) fits none, so an iterating half needs the growth dial.

**On its own merits:** a small, clean addition (a zoom dial) if other parts of the theory want sizes;
no new lead on the iterating half. All dials here are known math. No line for the README
"helpful to mathematicians" list.

## Mapping

| Fractal dial | In mult-add | Status |
|---|---|---|
| dimension d | fractional power (1+a)^d | has it |
| weights p : 1−p (D_q, f(α)) | per-element t = p/(1−p) | has it |
| branching growth | rung flow at t = ½ | has it, separately |
| wobble / complex dimensions | phase on zoom; like twist, Clifford rotor | partial (analogy) |
| zoom ratio r | none: mult-add has no size | missing: one candidate |
| placement (lacunarity, Minkowski content) | geometry | outside |
| topology (touching, gluing) | silent below d = 1 | outside |

## Numbers (checks.py, 25 of 25 pass)

Five dimension-½ sets: A two at 1/4; B three at 1/9 evenly; C three at 1/9, middle off-centre;
D uneven 9/100 and 49/100; E four at 1/16.

| | A | B | C | D | E |
|---|---|---|---|---|---|
| lattice | yes | yes | yes | no | yes |
| complex-dim spacing 2π/ln(1/r) | 4.53 | 2.86 | 2.86 | — | 2.27 |
| Minkowski ratio (avg) | 2.89 wobbles | 2.97 wobbles | 2.78 wobbles | 3.00 settles | 3.10 wobbles |
| lacunarity | 4.01 | 3.78 | 4.13 | 3.76 | 3.59 |

Multifractal α range on A: p = 0.3 → [0.26, 0.87]; p = 0.1 → [0.08, 1.66].
Four-corner set vs binary tree distance ratio: 0.50 … 2.05 through depth 4 (bi-Lipschitz).

## Lean (lean/FractalDials.lean → library Topics/FractalDials, 14 theorems; whole library builds clean)

half_dimension_certificates, mulAdd_shift, mulAdd_identity, two_pow_mod7, uneven_not_lattice,
weights_odds, pieces_join_closed, four_corner_is_tree, dyadic_glue, same_dimension_gasket_dust,
gasket_vertices, gasket_touches, powerset_branch, powerset_not_selfsimilar.

## Files

- `fractals.py` — the sets and measuring tools; `checks.py` → `checks_output.txt`, `results.json`
- `figures.py` → `figures/p1…p7*.svg` (and `png/` renders)
- `build_page.py` → `index.html` (self-contained page, plates inline)
- `lean/FractalDials.lean` (same file as lean-library/FusedAlgebra/Topics/FractalDials.lean;
  library edits: build.sh, FusedAlgebra.lean, README.md, theorems.tsv)

Not pushed to the repo (Brennan will arrange that).

## Sources

Moran 1946; Lapidus & van Frankenhuijsen, *Fractal Geometry, Complex Dimensions and Zeta Functions*
(2013); Allain & Cloitre, Phys. Rev. A 44 (1991); Halsey et al., Phys. Rev. A 33 (1986); Falconer,
*Fractal Geometry*; Brouwer 1910; Kneser, J. reine angew. Math. 187 (1950); Paulsen & Cowgill,
Adv. Comput. Math. 43 (2017).
