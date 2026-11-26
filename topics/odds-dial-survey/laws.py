"""The 65 law groups (topics/dial-catalog/catalog.py) and what turning the odds t does to each.
t is one number per element: how strongly that element takes part (t = 0 drops it, t = 1 plain power set,
t -> infinity forces it). Two ways to turn it: all elements together (uniform) or one element alone.
Run: python3 laws.py  -> laws.json, prints the counts."""
import json, sys
sys.path.insert(0, "../dial-catalog")
from catalog import GROUPS

V = {"slides": "carries a proved law into another proved law; link checked along the whole dial",
     "plausible": "numbers slide; the law in between not proved",
     "partly": "survives on part of the dial, breaks on the rest",
     "disguise": "the move is really another dial",
     "self": "the group is the odds dial itself",
     "no": "no slide (t is only the law's argument, or the law lives on another dial)"}

# group index (1-based) -> (verdict, candidate, one-line reason grounded in the repo)
R = {
 1: ("no", "", "t is what the circle operation composes (DialSurvey.odds_union); the law holds at every t, so nothing moves"),
 2: ("self", "", "with x_i = t_i the corner expansion IS the dialled power-set count (PowersetDials.total_W_chain)"),
 3: ("slides", "C1", "uniform t = 2^a - 1 carries 1 + t^(n) = (1+t)^n onto the exponent ruler 2^(an) (DialSurvey.exponent_ruler)"),
 4: ("no", "", "rungs live on the rung position; only the circle-operation rung's fractional powers meet t, and that is C1"),
 5: ("no", "", "calculus laws do not depend on a per-element weight"),
 6: ("no", "", "plane law (k, d) is a two-dial law; t is its argument"),
 7: ("no", "", "Clifford signatures live on d per axis"),
 8: ("no", "", "torus laws; no per-element weight"),
 9: ("slides", "C2", "odds (1+v)/(1-v) is the Doppler lift; multiplying t is a boost (GameTheory.odds_is_lift, InformationTheory.odds_to_velocity)"),
 10: ("slides", "C2", "multiplying one element's odds by L is Bayes' update (Statistics.bayes_odds, evidence_circ, repeated_evidence)"),
 11: ("no", "", "swap rule is d x area; no per-element weight"),
 12: ("no", "", "Maxwell; no per-element weight"),
 13: ("no", "", "Einstein tensor; no per-element weight"),
 14: ("no", "", "causal-diamond counting is a shape (C) / order question"),
 15: ("self", "", "PowersetDials, FractionalPowerset, FractalDials define the per-element odds (total_W, conv_W, weights_odds)"),
 16: ("no", "", "set partitions and Catalan counts carry no per-element odds in the Lean statements"),
 17: ("slides", "C9", "geometric odds t_i = q^i t turn the binomial count into the q-binomial theorem (QAcross.QCounting.q_binomial_theorem)"),
 18: ("slides", "C9", "odds t_i = -q^i, n -> infinity, give Euler's product (QAcross.Partitions.pentagonal_30, euler_is_circ_30)"),
 19: ("no", "", "toric charts are shapes; the genus slide belongs to the shape and plane dials"),
 20: ("slides", "C10", "one bond's odds 0 -> infinity slides deletion into contraction (GraphTheory.wheatstone, series_odds)"),
 21: ("self", "", "independent OR is the circle operation on odds (MeasureProbability.odds_or), exactly how two power-set dials merge"),
 22: ("slides", "C3", "chance = t/(1+t) runs 0 -> 1/2 -> 1; cube vs simplex reading of the same odds (Statistics, compl_fixed_dial)"),
 23: ("slides", "C2", "belief propagation in this module: posterior odds multiply (InformationTheory.variable_node); Renyi itself does not move"),
 24: ("slides", "C3", "entropy of the dialled power set is n h(t/(1+t)): 0 at t = 0 and infinity, n bits at t = 1 (Thermodynamics.landauer_states 2^n)"),
 25: ("slides", "C4", "occupation n = t/(1+kt): small odds bring Fermi, Boltzmann, Bose together (RepetitionDial.inverse_gap, fluctuation)"),
 26: ("slides", "C4", "Boltzmann odds in the Fermi function (Thermodynamics.fermi_velocity_num, boson_mode, lens_fixed_point)"),
 27: ("slides", "C5", "uniform odds is the fugacity: zeros on the circle |normalised t| = 1, empty phase <-> full phase (LeeYang, lee-yang page)"),
 28: ("no", "", "Born rule is a p = 2 / p = 1 matter (shape-dial survey C3)"),
 29: ("no", "", "CHSH lives on the roundness dial; the 22.5 degree coincidence is checked in N2 and is not a law"),
 30: ("no", "", "region-entropy cones and soft min cuts move on k and a softness (Maslov) dial, not on t"),
 31: ("slides", "C2", "selection multiplies an allele's odds by (1+s) each generation; Elo odds multiply (Biology.Evolution.generations_add, GameTheory.odds_circ)"),
 32: ("slides", "C8", "tilting a point process by t^N multiplies each mode's odds mu/(1+k mu) by t: same family (DialPointProcesses; L-ensembles)"),
 33: ("no", "", "strain triangle; no per-element weight"),
 34: ("no", "", "no per-element weight"),
 35: ("no", "", "no per-element weight"),
 36: ("no", "", "Fermat lives on the roundness dial"),
 37: ("no", "", "Frey pieces; no per-element weight"),
 38: ("no", "", "topology; no per-element weight"),
 39: ("no", "", "no per-element weight"),
 40: ("slides", "C4", "n elements at odds 1/n count (1 + 1/n)^n -> e (E.stair, stair_le_s): the small-odds (Poisson) limit"),
 41: ("no", "", "pi and the Gaussian; no per-element weight in the Lean statements"),
 42: ("no", "", "complex units; t enters only as a complex fugacity, which is C5"),
 43: ("no", "", "conservation laws do not depend on t"),
 44: ("slides", "C3", "per-site binding odds [L]/K: saturation t/(1+t), cooperative sites have no real zero (Chemistry.ChemEquilibrium.binding_dial, cooperative_no_real_zero)"),
 45: ("no", "", "Bliss and Loewe are the union law and addition on effect odds: they hold at every t (same as group 21)"),
 46: ("no", "", "Fisher's real-rate rule; t is the argument"),
 47: ("no", "", "Hamacher's family acts on odds with k = gamma (LogicSets.hamS_odds): a k matter, t is its argument"),
 48: ("no", "", "no per-element weight"),
 49: ("no", "", "stage overhead odds stack by the circle operation; t is the argument"),
 50: ("no", "", "whole rungs; no per-element weight"),
 51: ("no", "", "Mobius products live on (k, d)"),
 52: ("no", "", "number planes live on (k, d)"),
 53: ("no", "", "jet fields live on (k, d)"),
 54: ("no", "", "order dial q"),
 55: ("no", "", "conventions; no law moves"),
 56: ("no", "", "insertion cubes weighted by potential steps, not by per-element odds"),
 57: ("no", "", "for solitons a uniform t in det(I + tM) is a gauge: it only slides the waves (corner-weights README)"),
 58: ("slides", "C7", "odds X of the two RT phases: disconnected (I = 0) <-> connected (I = log X), glued by X1 o_1 X2 (HoloGluing.monogamy_glued)"),
 59: ("no", "", "force reach lives on (k, d)"),
 60: ("no", "", "barrier odds R/T compose by o_K; t is the argument (Tunnelling)"),
 61: ("slides", "C6", "t_p = p^(-s): power set of primes -> squarefree numbers, zeta(s)/zeta(2s) (PowersetDials.euler_squarefree, fermi_dirac)"),
 62: ("no", "", "elliptic modulus; no per-element weight"),
 63: ("no", "", "plane-law flows"),
 64: ("no", "", "Boolean lattice and fuzzy frame: t-free statements"),
 65: ("no", "", "energy conditions; the TFD cross-ratio odds is an argument"),
}
assert len(GROUPS) == 65 and set(R) == set(range(1, 66))
rows = []
for i, (law, mods, shape) in enumerate(GROUPS, 1):
    v, c, why = R[i]
    rows.append(dict(id=i, law=law, modules=mods, verdict=v, candidate=c, reason=why))
counts = {k: sum(r["verdict"] == k for r in rows) for k in V}
json.dump({"dial": "t, the odds", "verdicts": V, "rows": rows, "counts": counts}, open("laws.json", "w"), indent=1)
print(len(rows), "groups;", counts)
