"""Ratio-dial survey (Nov 26, 2026): what turning the similarity ratio r does to each of the 65 law groups.
All other dials held fixed. Verdicts: slides, plausible, partly, disguise, self, no (see ../../survey-spec rules in README).
Run: python3 laws.py  -> laws.json and the counts."""
import sys, os, json
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, "..", "dial-catalog"))
from catalog import GROUPS

NO_PLANE = "plane-dial law in (k, d); a piece's size never enters"
NO_FREE = "does not involve a piece with a size; no r to turn"
V = {  # group index (0-based) -> (verdict, candidate, reason)
 0: ("disguise", "C3", "k dial in disguise: a weight w = r^s on each element gives w(t o_w u) = (wt) o_1 (wu) (Base circK_scale = ZoomDial.zoom_is_twist); at the Moran dimension the uniform weight r^D = 1/N does not move at all"),
 1: ("disguise", "C3", "k dial in disguise: scaling every element by w turns k^(|S|-1) x^S into w (wk)^(|S|-1) x^S, a change of k (circK_scale)"),
 2: ("slides", "C2", "(1+a)^n slides to the fractional circle power (1+a)^D as r moves D = n log2/log(1/r); ends proved (Base powers, FractalBoxes.cube_layers, half_layers); Newton's binomial carries it exactly"),
 3: ("self", "C6", "ZoomDial is the r dial itself (Moran's wall, words_weigh, zoom_compose); the rung and height-flow laws in this group meet r only as the clock of the flow, which is the d dial in disguise for rho (newton-g P8, Q1)"),
 4: ("no", "", NO_FREE),
 5: ("no", "", NO_PLANE),
 6: ("no", "", NO_PLANE),
 7: ("no", "", "torus pinching is a plane-dial (d) move; no piece sizes"),
 8: ("no", "", "velocity addition is k = 0, d = -1/c^2; no piece sizes"),
 9: ("no", "", "odds multiply; the odds dial t, not a size"),
 10: ("no", "", "holonomy = curvature x area; area of a curved region, not a self-similar piece"),
 11: ("no", "", NO_FREE),
 12: ("no", "", NO_FREE),
 13: ("no", "", "Myrheim-Meyer counting reads a manifold's whole-number dimension; no similarity ratio enters"),
 14: ("partly", "C1", "FractalDials: two pieces at ratio r give D = log2/log(1/r): half box (1/4), Cantor (1/3), segment (1/2) on one closed form; past r = 1/2 the pieces overlap and D stops being the set's dimension (Hutchinson 1981); power-set dials do not depend on r"),
 15: ("no", "", "Bell, Catalan and partition counts are corner-weight (shape) laws; no sizes"),
 16: ("no", "", "[n]_q is a q (curvature-side) count; a size ratio does not enter"),
 17: ("no", "", "q^area weights are the order/curvature side; no self-similar pieces"),
 18: ("no", "", "toric charts are counted on the simplex; no sizes"),
 19: ("no", "", NO_FREE),
 20: ("no", "", "union rule is k = -1; no sizes"),
 21: ("no", "", "logistic/softmax live on the odds; r only enters through C3's weights"),
 22: ("plausible", "C4", "Renyi entropy becomes the Renyi dimension D_q = H_q(p)/log(1/r) at a uniform ratio (exact); with unequal ratios the multifractal formalism (Cawley-Mauldin 1992) carries it, not proved in the repo"),
 23: ("slides", "C4", "Kraft sum 1 for every full binary tree (InformationFull.kraft_one, r = 1/2, D = 1) slides to the Moran wall (ZoomDial.wall_every_depth) for any branch ratios at their D; checked on random trees along r"),
 24: ("partly", "C5", "Gentile cap (repetition-dial: at most m copies = any-number element at size 1 over one at size m+1): leftover share r^beta with r = 1/(m+1); Fermi (r = 1/2, 7/8 at beta = 3) to Bose (r -> 0); holds only at r = 1/(m+1)"),
 25: ("disguise", "C3", "temperature in disguise: r_i = exp(-E_i) are Boltzmann factors and Moran's equation says the free energy is 0 (Bowen's equation, 1979); D comes out as the inverse temperature, a reading"),
 26: ("no", "", "Lee-Yang zeros of spin models; no self-similar pieces"),
 27: ("no", "", "Born rule lives on the sphere; no sizes"),
 28: ("no", "", "CHSH is a shape (p) law"),
 29: ("no", "", "entropy cones are polymatroid shapes; no sizes"),
 30: ("no", "", "replicator flow on the simplex; no sizes"),
 31: ("no", "", NO_FREE),
 32: ("no", "", "strain triangle; self-similar blow-up ansatz is not a law proved here (ns-same-shape)"),
 33: ("no", "", NO_FREE),
 34: ("no", "", "newton-g's clock r is treated under C6 (group 4 and 64)"),
 35: ("no", "", "Fermat curves are p-circles; no sizes"),
 36: ("no", "", NO_FREE),
 37: ("no", "", NO_FREE),
 38: ("no", "", NO_FREE),
 39: ("no", "", "corner volumes 1/n! of the simplex; whole-number dimension only"),
 40: ("no", "", "ball volumes; whole-number dimension only"),
 41: ("no", "", NO_FREE),
 42: ("no", "", NO_FREE),
 43: ("no", "", NO_FREE),
 44: ("no", "", NO_FREE),
 45: ("no", "", NO_FREE),
 46: ("no", "", NO_FREE),
 47: ("no", "", NO_FREE),
 48: ("no", "", NO_FREE),
 49: ("no", "", "whole rungs of the ladder; no sizes"),
 50: ("no", "", NO_PLANE),
 51: ("no", "", NO_PLANE),
 52: ("no", "", "jet fields are a k setting; no sizes"),
 53: ("no", "", "holonomy is a reading of the order dial q; no sizes"),
 54: ("no", "", "states that a weight is a choice of units; C3 leans on it, and turning r leaves it as it is"),
 55: ("no", "", NO_FREE),
 56: ("slides", "C2", "fractional boxes (FractalBoxes: half_layers, half_join_stick, third_join_stick) slide along D as r turns; joining at a common ratio adds D (Moran); the size weight part is C3 (disguise), corner weights do not depend on r"),
 57: ("no", "", "monogamy from the odds; no sizes"),
 58: ("no", "", NO_PLANE),
 59: ("no", "", NO_PLANE),
 60: ("no", "", "prime products; the repo's complex dimensions (zoom-dial) do not touch Lapidus-Maier (see the fractal-strings section)"),
 61: ("no", "", "kappa = 1 - 2m is a d move; no sizes"),
 62: ("disguise", "C6", "for the two-loop number rho, changing the flow's clock r is the curvature dial d in disguise: rho = (ln r)/3 = 1/3 - 4d/(3k^2) at d = k^2(1 - ln r)/4 (newton-g P8, Q1); rho = 0 is clock-free (Q4)"),
 63: ("no", "", NO_FREE),
 64: ("no", "", NO_FREE),
}
assert len(GROUPS) == 65 and set(V) == set(range(65))
rows = []
for i, (law, mods, shape) in enumerate(GROUPS):
    v, c, why = V[i]
    rows.append({"n": i + 1, "law": law, "modules": mods, "verdict": v, "candidate": c, "reason": why})
json.dump(rows, open(os.path.join(here, "laws.json"), "w"), indent=1, ensure_ascii=False)
from collections import Counter
cnt = Counter(r["verdict"] for r in rows)
for k in ("slides", "plausible", "partly", "disguise", "self", "no"):
    print(f"{k:10s} {cnt[k]:3d}  " + ", ".join(str(r['n']) for r in rows if r['verdict'] == k and k != 'no'))
print("total", sum(cnt.values()))
