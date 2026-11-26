"""The 65 law groups (topics/dial-catalog/catalog.py) sorted by what turning the order dial q does to each,
every other dial held fixed. Run: python3 laws.py -> laws.json, prints the counts."""
import json, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "dial-catalog"))
from catalog import GROUPS

V = {"slides": "carries a proved law into another proved law; link checked along the dial",
     "plausible": "numbers slide, the law in between is not proved here",
     "partly": "survives on part of the dial, breaks on the rest",
     "disguise": "the move is really another dial",
     "self": "is the order dial itself",
     "no": "no slide (another dial, or does not depend on q)"}

# index (0-based in GROUPS): (verdict, candidate, one-line reason)
S = {
 0: ("slides", "C1", "Over q-commuting inputs o_k stays associative (A1); order gap k(1-q)ab (A2): order-free at q = 1 "
                    "(Lean circK_comm), exclusion at q = -1 (Lean anticommute_self_zero, pauli)"),
 1: ("slides", "C2", "Run the chain in reversed order over x_j x_i = q x_i x_j: rank-m corners get q^C(m,2) (B1), the "
                    "soliton corner form (Lean tau3_equal_corner); sorted order keeps the plain expansion (B2)"),
 2: ("slides", "C3", "(x + y)^n with yx = q xy is the q-binomial (Potter 1950, Schutzenberger 1953; C1); q = 1 binomial "
                    "(Lean gauss_at_one), q = 0 sorted words, q -> inf reversed word, roots of unity q-Lucas"),
 3: ("no", "", "Rungs live on the half-line and the k dial; no product of two moves to reorder"),
 4: ("disguise", "C12", "Jackson's q-derivative gives [n]_q x^(n-1) (L4), but its q = 1 + h is a step size: the k dial"),
 5: ("disguise", "C12", "The face weight q = beta/alpha is a reading of (k, d): (1+q)^2/q = k^2/d (L1); turning it is "
                       "turning d (curvature dial in disguise)"),
 6: ("partly", "C9", "Elements that square to numbers allow only q = +1 or -1 (I1, order-dial G1); in between, Clifford breaks"),
 7: ("no", "", "The Tate parameter q of the torus (QDial tate_glue_n) is the curve's period, not a swap factor"),
 8: ("no", "", "Velocity addition is a ball law on the plane dial; no q"),
 9: ("no", "", "Bayes in odds form is commutative and has no swap factor"),
 10: ("disguise", "C12", "Turn = d x area flattens into phase = area/hbar with hbar = 1/(j d) (quantum-action A3): "
                        "the d dial and spin, not q at fixed d"),
 11: ("no", "", "Maxwell's law does not depend on q (magnetic translations obey a Weyl rule, but that is not proved here)"),
 12: ("no", "", "Einstein tensor: no swap factor"),
 13: ("no", "", "Causal-diamond counting: no swap factor"),
 14: ("slides", "C6", "Routes through the cube (the permutohedron's vertices) weighted q^inv sum to [n]_q! (F1; "
                     "cube-weights E1): n! at q = 1, one route at q = 0, all cancel at a cube root of 1 from n = 3 (F2)"),
 15: ("no", "", "q-Catalan and q-Stirling numbers exist in print (Carlitz) but are not built in this repo"),
 16: ("slides", "C4", "[n]_q = sum over words with one y of q^inv (D1): n at q = 1, 1 at q = 0, parity at q = -1, 0 at a "
                     "root of unity (Lean QDial); same numbers as the k dial at k = q - 1 (Lean qnum_lift)"),
 17: ("slides", "C5", "Face weights q^area: inversions = path area (E2), sum = [n,m]_q (E1): counts C(n,m) at q = 1 (Lean "
                     "gauss_at_one), keeps the one minimum at q -> 0 (E3; Lean theta2_face, orbit_product_q0)"),
 18: ("disguise", "C12", "Toric point counts run through a = q - 1 points on the circle (Lean fan_count): q is the field "
                        "size, reached through the k dial"),
 19: ("plausible", "C10", "Jones polynomial: Lean proves t = 1 (components) and t = -1 (determinant); t = e^{2pi i/3} "
                         "and t = i checked (J2, J3); that t is the order dial goes through U_q(sl2), not proved here"),
 20: ("no", "", "The union rule is the k dial at k = -1"),
 21: ("no", "", "Logistic and softmax are commutative sums"),
 22: ("no", "", "Tsallis's q is k + 1 (the homothety dial), the first clashing meaning of q, not the order dial"),
 23: ("no", "", "Shannon and code lengths: commutative"),
 24: ("partly", "C9", "Repeats survive only at q = 1; any q != 1 on self-swaps forces x^2 = 0 (I2, Lean anticommute_self_zero):"
                     " Bose to Fermi is a jump; anyons in between not proved"),
 25: ("no", "", "Gibbs and temperature are the Maslov (k) dial"),
 26: ("no", "", "Lee-Yang: no swap factor"),
 27: ("no", "", "Born rule: no swap factor"),
 28: ("no", "", "CHSH uses anticommuting observables (q = -1), but nothing slides with q in the repo"),
 29: ("no", "", "Entropy cones: commutative"),
 30: ("no", "", "Replicator dynamics: no swap factor"),
 31: ("no", "", "Point processes are the k dial (alpha-determinantal)"),
 32: ("no", "", "Navier-Stokes: no swap factor"),
 33: ("no", "", "Dark sector: no q"),
 34: ("no", "", "Newton's G: no q"),
 35: ("no", "", "Fermat: no q"),
 36: ("no", "", "Frey curves: no q"),
 37: ("no", "", "Topology modules: no q"),
 38: ("no", "", "P vs NP: no q"),
 39: ("plausible", "C11", "e_q(x + y) = e_q(x) e_q(y) for yx = q xy (Schutzenberger 1953; K1 to degree 7): e at q = 1 "
                         "(Lean E), 1/(1 - z) at q = 0; the q-law is not proved in the repo"),
 40: ("no", "", "pi and the Gaussian: no q"),
 41: ("no", "", "Euler's identity marks the stop q = e^{i pi} = -1, but the law itself has no q"),
 42: ("no", "", "Conservation laws: commutative"),
 43: ("no", "", "Chemistry: no q"),
 44: ("no", "", "Biology: no q"),
 45: ("no", "", "Fisher's rule: commutative (compounding order is different centres, order-dial section 2)"),
 46: ("no", "", "The Boom hierarchy is order on/off, not a dial (order-dial section 7)"),
 47: ("no", "", "Functional analysis etc.: no q"),
 48: ("plausible", "C7", "Rulecraft's order of reforms weighted q^(out-of-order pairs): best order at q = 0 (H5a), plain "
                        "average at q = 1, worst at q -> inf; monotone in between checked (G2, G3), not proved here"),
 49: ("no", "", "Whole rungs: k dial"),
 50: ("disguise", "C12", "Mobius translations' multiplier ratio is q = beta/alpha (anchors, rulecraft-history H15): the "
                        "d dial in disguise"),
 51: ("no", "", "Number planes write u^2 = k u - q, where q is the curvature dial d (a third clashing meaning)"),
 52: ("no", "", "Jet fields: no swap factor"),
 53: ("self", "C8", "The order dial: swap gap = oriented area (Lean order_gap), Weyl phase = area (Lean swap_costs_area, "
                   "weyl), Moyal at hbar (order-dial H1-H5)"),
 54: ("no", "", "Units and roots: no q"),
 55: ("no", "", "Dial potential and checkerboard: k and d dials"),
 56: ("no", "", "Size weights: r and beta"),
 57: ("no", "", "Holographic gluing: odds"),
 58: ("no", "", "Force reach: plane dial"),
 59: ("no", "", "Tunnelling: o_k with k"),
 60: ("no", "", "Riemann rounds: no swap factor"),
 61: ("no", "", "Elliptic modulus kappa = 1 - 2m: the nome is a period, not a swap factor"),
 62: ("no", "", "Beyond one loop: plane dial"),
 63: ("no", "", "Boolean lattice: commutative"),
 64: ("no", "", "Warp bubbles: no q"),
}

rows = []
for i, (law, mods, shape) in enumerate(GROUPS):
    v, c, r = S[i]
    rows.append({"n": i + 1, "law": law, "modules": mods, "shape": shape, "verdict": v, "candidate": c, "reason": r})

if __name__ == "__main__":
    assert len(rows) == 65 == len(S)
    json.dump({"verdicts": V, "rows": rows}, open("laws.json", "w"), indent=1, ensure_ascii=False)
    from collections import Counter
    cnt = Counter(r["verdict"] for r in rows)
    for kk in V:
        print(f"{kk:10s} {cnt[kk]:3d}")
    print("total", sum(cnt.values()))
