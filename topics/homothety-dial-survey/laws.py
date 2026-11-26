"""The 65 law groups of topics/dial-catalog/catalog.py, and what turning the homothety k (d = 0, every other dial fixed)
does to each. Run: python3 laws.py -> laws.json, prints the counts.
Verdicts: slides, plausible, partly, disguise, self, no (see V)."""
import json, sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "dial-catalog"))
from catalog import GROUPS

V = {"slides": "carries a proved law into another proved law; the link is checked along the whole dial",
     "plausible": "numbers slide, the law in between is not proved",
     "partly": "survives on part of the dial, breaks on the rest",
     "disguise": "the move is really another dial",
     "self": "the group is the k dial itself",
     "no": "no slide (lives on another dial, or k is only a choice of units there)"}

# (verdict, candidate, one-line reason), in catalog order
J = [
 ("self", "", "circK_assoc, circK_comm, circK_zero hold at every k; circK_scale makes k a unit choice; the wall -1/k moves (A1-A4)"),
 ("slides", "C1", "1 + k(x1 o ... o xn) = prod(1 + k xi): inclusion-exclusion (-1), sums (0), products (1); circK_lift, kfold_lift, incl_excl3 (B1-B3)"),
 ("slides", "C1", "wpow_lift at every k; k = -1 is 1 - (1 - a)^n (wpow_neg_one), k = 0 is n a, k = 1 is (1 + a)^n - 1 (B4)"),
 ("no", "", "every k != 0 is the k = 1 ladder rescaled (fexpK_scale); k = 0 cannot be reached (lift_zero_not_injective)"),
 ("slides", "C5", "D_k F = F'/(1 + kF) obeys the product rule and the fundamental theorem at every k: ordinary (0), geometric/product (1), hazard (-1) (F1-F4)"),
 ("no", "", "lives on the curvature dial d; at d = 0 it is group 1 (turning k at fixed d != 0 crosses the parabola: the d survey's business)"),
 ("no", "", "signatures are readings of d per axis; k plays no part"),
 ("no", "", "the nome and the pinching live on d and the period lattice"),
 ("no", "", "velocity addition is (k, d) = (0, -1/c^2); turning k at that d keeps the same kind and only moves the walls (a change of coordinates)"),
 ("no", "", "odds multiply at k = 1 in odds units; any k != 0 is a unit change; log-odds adding is the ruler, not the k = 0 stop"),
 ("no", "", "turn = d x area: curvature dial and order dial"),
 ("no", "", "no k"),
 ("no", "", "no k"),
 ("no", "", "causal counting lives on the shape and order dials"),
 ("disguise", "", "corner weight k^(|S|-1) t^S = (k t)^S / k: turning k rescales the odds t (t dial in disguise)"),
 ("slides", "C9", "Touchard T_n(k): Bell numbers at k = 1, complementary Bell at k = -1 (dial_blocks_minus_one); rising product prod(1 + jk) (kfold_lift) (J4, J5)"),
 ("slides", "C9", "[n]_q = n copies of 1 under o_k with k = q - 1 (qnum_eq_wpow); stops proved in qnum_at_* (J1-J3)"),
 ("no", "", "the q of q^area is the order dial, not k"),
 ("no", "", "toric charts live on the shape (corner-weight) dial"),
 ("no", "", "knot and graph polynomials use their own q; the k used there is fixed at -1 or a unit choice"),
 ("slides", "C4", "Sugeno lambda-measure: normalisation pins lambda = k; plausibility side (-1 < k < 0), probability (0), belief side (k > 0); lam_* theorems (E1-E4)"),
 ("slides", "C1", "normaliser 1 + (prod(1 + kx) - 1)/k: softmax with a none option (k -> 0) to independent logits (k = 1) (B6)"),
 ("slides", "C3", "ln_q is the dial log with k = 1 - q_T; Tsallis composes by o_k, Renyi = L_k(Tsallis) is additive; tsallis_engine, renyi_additive (D1-D3)"),
 ("slides", "C3", "Shannon is the k = 0 stop of the Tsallis/Renyi slide (D4); Campbell's code lengths carry Renyi's order"),
 ("slides", "C2", "per-mode count = Katz/Panjer (a,b,0) class: binomial (k = -1/m), Poisson (0), negative binomial (k > 0); kRise_* and RepetitionDial (C1-C8)"),
 ("slides", "C2", "Fermi (-1), Boltzmann (0), Bose (+1) occupations 1/(e^x - k) in one formula; fermion_lift, boson_mode, boson_fermion_dual (C6)"),
 ("partly", "C8", "on a graph the edge weight k is not a unit: Lee-Yang holds for k >= 1 (attraction), fails for k < 1 (I1-I3)"),
 ("slides", "C10", "the noise channels: amplitude damping (-1), depolarizing (-4/3), bit flip (-2) are stops of k = -q/(q-1) (qary_cascade; K1-K4)"),
 ("no", "", "CHSH lives on the shape dial p"),
 ("partly", "C11", "corner sum with weight k^(rank-1): sign-definite on geometries only at k = -1 (MMI); positive for every k >= sqrt6 - 3; between, either sign (L1-L4)"),
 ("no", "", "replicator and Fisher's theorem do not depend on k; dial_game_union is C4's law"),
 ("slides", "C2", "alpha-determinantal processes are the per-mode law with k = alpha; k < 0 exists only at -1/m (Shirai-Takahashi 2003)"),
 ("no", "", "no k"), ("no", "", "no k"), ("no", "", "no k"),
 ("no", "", "Fermat lives on the roundness dial p"), ("no", "", "no k"), ("no", "", "no k"),
 ("slides", "C9", "witness fold [n]_(1+k): OR at -1 (NP), count at 0 (#P), parity at -2 (parity P); witness_fold, np_point, sharpP_point, parityP_point"),
 ("slides", "C5", "e = limit of (1 + k)^(1/k) as the grain k -> 0: compound to continuous (F6)"),
 ("no", "", "pi and the Gaussian live on the sphere (shape dial)"),
 ("no", "", "EIPi uses k only as a direction (dial_is_direction): a unit choice"),
 ("no", "", "uses k = -1 only"), ("no", "", "no k"),
 ("slides", "C10", "Haldane's map r1 + r2 - 2 r1 r2 is the q = 2 stop of the alphabet ladder; the cure model k = -1/pi (where-k-moves) sits on the same side"),
 ("slides", "C5", "Fisher's rule (1 + n) = (1 + r)(1 + i) with a period of length k: rates add exactly as k -> 0 (continuous compounding) (F5)"),
 ("partly", "C9", "bits stay bits only at k = -1 and k = -2 (idem_defect, bits_closed_iff); matrix walls at k = -1/lambda (det(1 + kA), Dk_hom) (J6, B5)"),
 ("no", "", "Neumann series is 'no wall while |k||a| < 1', a unit statement (mathlib-map)"),
 ("no", "", "same-target rules are k = -1 in strength units; any k < 0 is the same law; the history result lives on the order dial"),
 ("no", "", "every k != 0 is the k = 1 ladder rescaled (circK_scale, fexpK_scale); the partner a + b + 1/k blows up at k = 0"),
 ("no", "", "Mobius plane law lives on (k, d) with d != 0"),
 ("slides", "C5", "u^2 = k u (q = 0): split numbers for k != 0, dual numbers at k = 0 (split_point, dual_point)"),
 ("slides", "C5", "R[e]/(e^2 - k e): k = 0 jets (derivative), k = h two-point difference (ofTwoPt_twoPt, dual_units_add) (F7)"),
 ("no", "", "the order dial q"),
 ("self", "", "the group that proves k is units (ShapeDirections: F_(lam k)(a, b) = F_k(lam a, lam b)/lam)"),
 ("no", "", "dial_rescale, dial_flip: k is a unit choice there"),
 ("disguise", "", "size and corner weights: k rescales the per-element weight (t dial in disguise)"),
 ("no", "", "cross-ratios glue by o_-1 and odds by o_1; the flip x/(1 - x) links them, not a turn of k"),
 ("no", "", "conics on the d dial"),
 ("slides", "C7", "averaged law of cosines composes by o_k with k = -2K: piling-up (-2), variances add (0), mean resistance (+2) (H1-H4)"),
 ("partly", "C6", "Z_k(s) = prod_N zeta(Ns)^(M(k,N)/k) for every k, but it continues past Re s = 0 only at k = +-1 (Dahlquist 1952) (G1-G5)"),
 ("no", "", "kappa = 1 - 2m is an elliptic modulus, not k"),
 ("no", "", "plane-law flows with d != 0"),
 ("partly", "C9", "the Boolean lattice is the k = -1 point; inclusion-exclusion slides with C1, Boolean closure holds only at k = -1, -2"),
 ("no", "", "energy conditions do not involve k"),
]
assert len(J) == len(GROUPS) == 65
rows = [dict(law=g[0], modules=g[1], verdict=v, candidate=c, reason=r) for g, (v, c, r) in zip(GROUPS, J)]
counts = {key: sum(r["verdict"] == key for r in rows) for key in V}
json.dump({"dial": "k (homothety), d = 0", "verdicts": V, "rows": rows, "counts": counts}, open("laws.json", "w"), indent=1)
print(len(rows), "groups;", counts)
