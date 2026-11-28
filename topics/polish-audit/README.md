# DRAFT · Polishing step on every number search (Nov 26)

> **Names:** some words or symbols on this page (rung, wall, twist) also mean something else in standard math, physics or computing, or inside this project. [NAMING-COLLISIONS.md](../../NAMING-COLLISIONS.md) says which meaning is in play.

Asked for by Brennan: "add a polishing step ... to all of our methods, as long as the polishing step doesn't take too
much effort from the computer." Each search that stopped at an approximate answer got a cheap refinement (Newton steps,
a tight re-run from the found point, or a closed-form check), with its cost recorded. Machine note: the polish scripts are numpy-only;
scipy and matplotlib were installed afterwards (pip, in the sandbox) to rerun the dial search and redraw the old figures.

## Three-body (three-body/polish/)
- All 66 orbits closed exactly by Newton on the half-lap mirror condition, Taylor integrator in 34-digit arithmetic:
  26 to 32 correct digits (median 29), 10 min wall on 4 cores. No period or turn changed (max ΔT 2.5e-8).
- Re-reading loop words at finer sampling corrected 11 words: 60 orbits now in Li–Liao's 695-family list, 6 not (was 57 / 9).
- Page republished: https://claude.ai/artifact/2r1nsFubhVsJ1dsEenDwsi (section 9, Figures 10b, 12, 13; Figures 9 and 10
  redrawn from the polished orbits and re-read words by `polish/fig10_redraw.py`).

## Other folders (polish.py + polish_output.txt in each)
| folder | what changed | page |
|---|---|---|
| dial-potential | bug in quintessence.py history() bracket hid the hump fits: sombrero Δχ² 0.1 → 4.53, cosine 3.3 → 4.80; slope 0.80 → 0.82, 1σ 0.61–0.97, Δχ² 4.85, 14 of 43 | bracket fixed in quintessence.py, grid and profile polished in search.py, search rerun (184 s): the sombrero now fails the held-out test (−14.7; BAO alone drags it to w0 ≈ −0.40) and is ruled out; the circle stays alive (+14.6). README; page rebuilt from results.json and republished |
| biology/03-cell | 95% ranges for α were read off a 0.05 grid; now −0.76..0.46, −0.99..−0.83, −1.00..−0.77; plateau shift 0.012 | README/page.py; fits.json carries the polished ranges; figure f1 redrawn; page republished |
| dark-sector | free-M/L Newton total 685,295 → 685,294; NFW free-M/L 10,774 → 10,730; others confirmed | live page republished |
| why-i/e-i-pi/asymptotic-safety/beyond-one-loop | exponents turn real at N_V ≈ 846 (was ≈ 870 off a 15% grid) | README; figure b4 redrawn with the crossing bisected (≈ 846) |
| solitons-rungs | soliton jump gap 7.9e-5 → 1.3e-14 (grid argmax replaced) | no README |
| action | Snell crossing gap 3.2e-8 → 1e-16 | README Polish line |
| biology 04, 05, 06, 07, 10; information-theory-full; quantum-to-classical; ns-rungs; photosynthesis-twists; econ-forecasting | confirmed at quoted precision (some gained 8+ digits) | Polish lines where a README exists |

Closed-form, linear or already-tight searches in the rest of the 101 scripts need no polish (lists in the helpers' notes,
summarised in the thread report). Busy folders: see SUGGESTIONS-busy-folders.md (Q# and shape-dial survey already polish).
