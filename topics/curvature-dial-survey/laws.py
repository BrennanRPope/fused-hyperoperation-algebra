"""Every law group of topics/dial-catalog/catalog.py, and what turning the curvature dial d does to it
(k held fixed; the survey is done at k = 0 first, then k != 0).
Run: python3 laws.py  -> laws.json, prints the counts.
The plane law is a o b = (a + b + k ab)/(1 - d ab). Curvature = -Delta/4 with Delta = k^2 - 4d is a reading."""
import json, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "dial-catalog"))
from catalog import GROUPS

V = {"slides": "slides into another law (both ends proved, link checked along the whole dial)",
     "plausible": "numbers slide; the law in between is not proved",
     "partly": "survives on part of the d dial, breaks on the rest",
     "disguise": "the move is really another dial or parameter",
     "self": "is the curvature dial itself",
     "no": "no slide (lives on another dial, or does not depend on d)"}

# index in GROUPS -> (verdict, candidate, one-line reason grounded in the repo)
J = {
 0: ("slides", "C1", "d = 0 is the core law; turning d moves its wall -1/k and brings the second wall in from infinity (Landscape assoc_num/assoc_den, DialSweep)"),
 1: ("slides", "C2", "the corner sum becomes a ratio of two corner sums with Lucas weights U_|S|(k,d); at d = 0 U_n = k^(n-1) gives the cube back (checks B1-B3)"),
 2: ("slides", "C2", "1 + k a^(n) = (1 + k a)^n becomes V(a^(n)) A_n^2 = V(a)^n, the norm of the number plane (NumberPlane norm_mul; check B5)"),
 3: ("plausible", "C8", "at fixed k, d sets the rung height s = k^2/Delta (Landscape rung_equivalent); the formal residue Delta/(3k^2) = 1/(3s) matches the rung's, analytic match not proved"),
 4: ("no", "", "the Calculus module is built on o_k steps at d = 0; a plane-law calculus would go through the ruler dx/V(x), not in the repo"),
 5: ("self", "C1", "the (k, d) plane law; d is one of its two coordinates (Landscape, DialSweep, DialWall, DialSurvey)"),
 6: ("slides", "C4", "each generator squares to -d at k = 0: Clifford's omega^2 = +1/0/-1 (1873) is d = -1/0/+1 (Clifford division_cases, sig_examples, dual_nilpotent)"),
 7: ("slides", "C5", "the pinched cubic y^2 + a1 xy = x^3 + a2 x^2 sits at (k, d) = (-a1, -a2): split node, cusp, non-split node as d crosses the parabola (Elliptic node_param, walls_are_branches; elliptic checks 1f)"),
 8: ("slides", "C1", "d = -1/c^2 is Einstein, d = 0 Galileo, d > 0 tangent addition (NumberPlane split/dual/complex_point, GroupsLie theta_*)"),
 9: ("partly", "C6", "odds updating needs two anchors; it survives (in a Moebius chart) while Delta > 0 and breaks at the parabola (rulecraft-history H15, H16)"),
 10: ("slides", "C3", "tan(w/2) = d(u x v)/(1 - d u.v): no twist at d = 0, Thomas at d = -1, Foucault/Berry at d = +1 (TwistedBoosts galilean_no_twist, twist_half_angle)"),
 11: ("partly", "C7", "the boost-invariant operator is d_x^2 + d d_t^2: a wave equation only for d < 0, degenerate at d = 0 (Galilean limits), elliptic for d > 0 (electromagnetism: wave equation survives a dial boost only at d = -1/c^2)"),
 12: ("no", "", "the Einstein tensor and Lovelock's step read the same in every signature; Schur pins d constant in 3D+ (gr-bridges S5), it does not slide"),
 13: ("partly", "C7", "causal-diamond counting needs a light cone, which exists only on the two-wall side; every d < 0 gives the same counting (units), d = 0 flattens the cone"),
 14: ("no", "", "corner weights and power-set dials are the shape and odds dials"),
 15: ("slides", "C9", "the genus slide: at fixed k, d sets q = beta/alpha of 1 + kx + dx^2 = (1 + alpha x)(1 + beta x); h-vectors at q (shape-dial survey C5, Landscape chi_y_line)"),
 16: ("slides", "C9", "[n]_q with q = beta/alpha; the Lucas weight is U_n = alpha^(n-1) [n]_q (check B6)"),
 17: ("plausible", "C9", "cube face weights q^area ride the same q, but the law in between is not proved (as in the shape-dial survey)"),
 18: ("slides", "C9", "toric h-vectors [n+1]_q at q set by d/k^2: 1 at d = 0, corner count on the parabola"),
 19: ("no", "", "the diamond lattice's Potts Q moves k and d together along k + d = -1, with Delta = Q^2 >= 0 (diamond-lattice finding 2); d alone freezes Q"),
 20: ("partly", "C6", "the union rule's anchor 1 survives as a moved anchor for every d < 1/4 at k = -1 (a Moebius chart); at d = 1/4 and beyond no bounded slider"),
 21: ("no", "", "logistic vs softmax is the k dial; d only changes the chart of the logit"),
 22: ("no", "", "Renyi = L_k(Tsallis) is the k dial"),
 23: ("no", "", "Shannon and Campbell live on the k dial"),
 24: ("no", "", "Fermi / classical / Bose is k = -1, 0, +1 (DialWall); d is not used"),
 25: ("no", "", "Gibbs and the Maslov temperature are k-dial and tropical limits"),
 26: ("no", "", "Lee-Yang zeros: turning d at fixed k moves the root modulus 1/sqrt(d), off the circle except at d = 1"),
 27: ("no", "", "the Born rule lives on the shape dial (sphere vs simplex)"),
 28: ("no", "", "CHSH lives on the shape dial"),
 29: ("no", "", "entropy cones do not use d"),
 30: ("no", "", "replicator dynamics do not use d"),
 31: ("no", "", "consensus and dial point processes are k-dial laws"),
 32: ("no", "", "restricted Euler's discriminant 27R^2 + 4Q^3 is a cubic's, not the plane law's"),
 33: ("no", "", "the dark-sector notes do not use d"),
 34: ("plausible", "C8", "the second-loop number rho = 1/3 - 4d/(3k^2) slides with d (newton-g fig n7), so gravity's 3/5 is a choice d = -k^2/5, not a law"),
 35: ("no", "", "Fermat curves live on the roundness dial p"),
 36: ("no", "", "Frey curves reduce to nodes, never cusps (elliptic check 3g): a placement off the parabola, not a slide"),
 37: ("no", "", "Hopf invariant one / Bott: d picks the square of one generator, the theorem is about how many"),
 38: ("no", "", "P vs NP does not use d"),
 39: ("no", "", "corner volumes 1/n! use the simplex"),
 40: ("partly", "C10", "pi is the turning side's ruler length pi/sqrt(d); it goes to infinity at d = 0 and has no meaning on the two-wall side"),
 41: ("partly", "C10", "e^(pi u/sqrt d) = -1 holds exactly for d > 0; the dual (d = 0) and split (d < 0) exponentials never return"),
 42: ("no", "", "conservation laws do not use d"),
 43: ("no", "", "equilibrium laws live on a log axis"),
 44: ("no", "", "biology modules do not use d"),
 45: ("no", "", "Fisher's real-rate rule is the k line"),
 46: ("no", "", "logic, sets, categories, algebra do not depend on d"),
 47: ("no", "", "functional analysis, numerics, dynamics do not depend on d"),
 48: ("no", "", "data centres, Rulecraft, early science do not depend on d"),
 49: ("no", "", "the hyperoperation ladder is the height / k dial"),
 50: ("slides", "C6", "the plane law is a Moebius product; turning d carries hyperbolic -> parabolic -> elliptic one-parameter groups (GroupsLie mob_dial, theta_tan, theta_velocity, theta_dual)"),
 51: ("slides", "C4", "number planes u^2 = ku - d: split, dual, complex (NumberPlane split/dual/complex_point, norm_walls/norm_border/norm_nonneg_turning); Kepler's velocity space too (Hypersphere fock_on_sphere; gravity-from-theory N4-N6)"),
 52: ("slides", "C11", "jets R[e]/(e^2 - ke) are the d = 0 stop of R[u]/(u^2 - ku + d): two sample points -> one point with a derivative -> complex step (jet-fields 'two sample points k apart')"),
 53: ("no", "", "the order dial's oriented area lives on q; d enters only through the swap rule (group 11)"),
 54: ("no", "", "units: d -> lambda^2 d (Landscape units_den) says d's size is a unit choice at k = 0; used in every candidate, not a slide"),
 55: ("no", "", "the dial potential and checkerboard are not run along d here"),
 56: ("no", "", "size and box weights do not use d"),
 57: ("no", "", "holographic gluing does not use d"),
 58: ("slides", "C12", "carrier mass m^2 = Delta/4 = -d at k = 0: Yukawa (d < 0) -> Coulomb (d = 0) -> closed ruler / tower (d > 0) (ForceReach drift_removal, double_root_iff_border)"),
 59: ("slides", "C13", "a barrier is B_{-1}(u), a beam splitter B_{+1}(u); averaged, probabilities compose by o_k with k = -2d (Tunnelling mean_law_curvature, rung_two_barrier, rung_two_mixer)"),
 60: ("slides", "C14", "the heat flow on a pair of zeros keeps their sum (k) and moves their product (d) by 2t; crossing the parabola is a collision (RiemannRound3 heat_pair_disc, fused_disc, conj_disc)"),
 61: ("disguise", "C15", "kappa = 1 - 2m moves the modulus m (the t^4 coefficient, off the (k, d) plane); it meets the d road only at kappa = +1 and -1 (EllipticModulus)"),
 62: ("slides", "C1", "the flow x' = 1 + d x^2: tanh -> linear -> tan; wall exponents +-2 sqrt(-d) shrink to 0 at d = 0 (BeyondOneLoop plane_generator, plane_exponents_opposite)"),
 63: ("no", "", "Boolean lattice and fuzzy frame are k = -1 laws"),
 64: ("no", "", "energy conditions are not run along d here"),
}
assert len(GROUPS) == 65 and set(J) == set(range(65))
rows = []
for i, (law, mods, shape) in enumerate(GROUPS):
    v, c, why = J[i]
    rows.append(dict(n=i + 1, law=law, modules=mods, verdict=v, candidate=c, reason=why))
counts = {k: sum(r["verdict"] == k for r in rows) for k in V}
json.dump({"dial": "d (curvature dial), k fixed", "verdicts": V, "rows": rows, "counts": counts},
          open("laws.json", "w"), indent=1)
print(len(rows), "groups;", counts)
