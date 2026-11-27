"""DRAFT. Catalogue of the walls the project has run into (gathered Nov 27 2026 from the topic
folders and task notes).  Each entry: where the proof or test stopped, which dial and which
operation, the kind of break (growth, number theory, topology, positivity), and the shape(s)
it belongs to:
  I   turn by i       positivity proved on the real k axis is missing or false at imaginary k
  II  many scales     two or more independent stretch factors (rank >= 2)
  III depth           each size works, the whole family does not (pairs fine, triples the wall)
  IV  no daylight     the law is a disguised copy of something already known (rigidity)
'proved' says who proved the stopping point, where it is proved at all.
Writes walls.json.
"""
import json

W = [
 # name, topic, where it stopped, dial, operation, break, shapes, proved
 ("New ties: 19+1=20", "new-ties",
  "Every length reduces to Conjecture 21.2 (2, 5, 19 balance only in 19+1=20). Each size is closed by Baker/LLL/2-adic search; no bound covers all lengths.",
  "s", "sum of prime powers", "number theory", ["II", "III"],
  "each size: ours (Baker-type bounds are Baker, LLL is Lenstra-Lenstra-Lovasz); all lengths: open"),
 ("New ties: S-unit terms grow", "new-ties",
  "Length n gives an S-unit equation with n terms; Evertse-Schlickewei-Schmidt counts depend on n, so no uniform finish.",
  "s", "sum", "number theory", ["II", "III"], "count bound: Evertse-Schlickewei-Schmidt 2002"),
 ("New ties: v19 at five letters", "new-ties",
  "The 19-adic valuation is unbounded at five letters, so the 2-adic trick that closed shorter words does not.",
  "s", "sum", "number theory", ["II", "III"], "ours (explicit families)"),
 ("Riemann: local vs global", "riemann",
  "The elliptic object is correct near each point but does not glue to a global positivity statement.",
  "k", "zeta product", "positivity", ["I"], "ours (numerical + by hand)"),
 ("Riemann: rewind", "riemann",
  "Running de Bruijn-Newman (backward heat flow) below t=0 loses real zeros.",
  "t", "heat flow", "positivity", ["I"], "Rodgers-Tao 2020 (Lambda >= 0)"),
 ("Riemann: even half", "riemann",
  "The flip-odd half keeps every zero on the line; the flip-even half pushes zeros off.",
  "k", "even/odd split", "positivity", ["I"], "ours (numerical, n = 8 to 14)"),
 ("Riemann: Knauf chain is many-body", "riemann",
  "The chain is ferromagnetic but many-body, so pair Lee-Yang theorems do not apply; the magnetic circle holds only above a temperature threshold (fits near beta = 2), and RH sits at beta = 3/2.",
  "p", "spin chain", "positivity", ["I", "III"],
  "chain: Knauf 1993; threshold numbers ours; why: Ruelle 2010 Thm 9 (all-temperature Lee-Yang forces pair interactions)"),
 ("Riemann: limit one-sidedness", "riemann",
  "Mirror construction puts all finite zeros on the line, but the limit's one-sidedness is equivalent in strength to RH.",
  "k", "mirror", "positivity", ["I"], "known restatement (de Branges / Lagarias)"),
 ("Riemann: cusp tail", "riemann",
  "Best lead; no proof.", "k", "cusp", "positivity", ["I"], "open"),
 ("Lee-Yang p-ball, odd p > 8", "shape-dial-survey",
  "Even p proved (Polya multiplier); p <= 8 proved one p at a time (boxes); no uniform near zone for non-even p > 8 because of the gamma factors in the numerator.",
  "p", "single-site measure", "positivity", ["I", "III"], "even p: Polya; p <= 8: ours"),
 ("Lee-Yang p-ball, every p", "shape-dial-survey",
  "Hermite-Biehler start for every p (Kumaraswamy characteristic function) unproved.",
  "p", "single-site measure", "positivity", ["I"], "open"),
 ("Reflection positivity off the real axis", "rp",
  "RP holds for real k, fails off it (negative at k = -i). The remaining assumptions are independent (charge clock, mean-energy clock, superselection, spin factor V7).",
  "k", "link window 1 + k eps E >= 0", "positivity", ["I"], "ours (counterexamples)"),
 ("Klein flow peak spin", "klein",
  "Peak spin keeps rising with the grid (151 -> 233 -> 313.8, check B9 fails) while enstrophy converges.",
  "d", "flow", "growth", ["III"], "ours (numerical)"),
 ("Frame drag", "frame-drag",
  "An exact rewrite of Kerr; nothing beyond it.", "k", "rotation", "topology", ["IV"], "ours"),
 ("Quantum needs i", "why-i",
  "No path measure at imaginary k; hbar enters only as area/hbar.",
  "k", "path sum", "positivity", ["I"], "Cameron 1960 (no Feynman measure)"),
 ("Fractal boxes: iterating", "fractal-boxes",
  "Iterating needs a half-exponential: not unique, no exp-log solution.",
  "s", "iteration", "growth", ["III"], "known (Kneser; no exp-log half-exponential)"),
 ("Fractal boxes: many ratios", "fractal-boxes",
  "Finite type for one ratio (rank 1); neighbour types grow when ratios are independent (rank 2 and 3).",
  "s", "stretch group", "number theory", ["II"], "rank 1: Ngai-Wang 2001; our counts in rank.json"),
 ("Hyperoperation ladder", "hyperoperation-ladder",
  "Commutative and associative only at rungs 1 and 2; fractional heights form a different ladder.",
  "r", "rung", "growth", ["III", "IV"], "classical (rung 3 not associative)"),
 ("Hypersphere product", "hypersphere",
  "A size-keeping product exists only in dimensions 1, 2, 4, 8; only p = 2 keeps every property.",
  "p", "product", "topology", ["IV"], "Hurwitz 1898; Adams 1960"),
 ("Conley index", "conley",
  "The index cannot read stability (parity law).", "k", "flow", "topology", ["IV"], "ours (parity law)"),
 ("Dial audit", "dial-audit",
  "Scalar laws carry only kappa = d - k^2/4.", "k, d", "scalar law", "topology", ["IV"], "ours (by hand)"),
 ("New force reach", "force-reach",
  "The hold does not single out l_D; its strength is set by particle content.",
  "d", "screening", "growth", ["IV"], "ours"),
 ("G stability", "G notes",
  "No stability theorem.", "k", "law", "positivity", ["I"], "open"),
 ("Three bodies", "pairs-vs-triples census",
  "Pairs solve in closed form, three do not.", "k", "composition", "growth", ["III"], "classical (Poincare)"),
 ("Negatives: Gram, pendulum swap, half-area, temperature mirror", "riemann",
  "Gram's law fails; the dual frame breaks; -1/2 is a value at a point, not a location; mirror is RH-strength in the limit.",
  "k", "various", "positivity", ["I"], "ours"),
]

rows = [dict(name=a, topic=b, stopped=c, dial=d, operation=e, break_type=f, shapes=g, proved=h)
        for a, b, c, d, e, f, g, h in W]
if __name__ == "__main__":
    json.dump(rows, open("walls.json", "w"), indent=1)
    from collections import Counter
    print(len(rows), "walls")
    print("by break:", Counter(r["break_type"] for r in rows))
    print("by shape:", Counter(s for r in rows for s in r["shapes"]))
