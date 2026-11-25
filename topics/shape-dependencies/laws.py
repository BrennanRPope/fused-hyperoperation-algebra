"""laws.py: the science laws in the library, their textbook dependencies, collapsed to shapes.

Nodes  : library modules that each state one body of law (physics, chemistry, biology).
Edges  : A -> B when B is derived from, or needs, A. Tag "textbook" = standard teaching and
         derivation order; "reading" = my own reading, flagged on the page; "lean" = the library's
         own imports agree. Hand-made: this table is the only hand step.
Shapes : each law's body and family = the most common among its module's theorems in the
         proof gallery (rerun on all 3,003 theorems), with the share so a weak majority shows.
Output : laws.json, laws_output.txt."""
import json, sys, random, collections, pathlib
HERE = pathlib.Path(__file__).resolve().parent
G = json.load(open(sys.argv[1]))
mods = collections.defaultdict(list)
for p in G["proofs"]:
    mods[p["module"].split("Topics.")[-1]].append(p)

N = {  # id: (label, [modules], field)
 "prod": ("Product laws (Newton, Ohm, Boyle, Kepler)", ["EarlyScience"], "physics"),
 "thermo": ("Thermodynamics", ["Thermodynamics"], "physics"),
 "em": ("Electromagnetism (Maxwell)", ["Electromagnetism"], "physics"),
 "sr": ("Special relativity", ["Relativity"], "physics"),
 "thomas": ("Twisted boosts (Thomas)", ["TwistedBoosts"], "physics"),
 "gr": ("General relativity", ["GeneralRelativity", "GRBridges"], "physics"),
 "gpb": ("Gravity Probe B", ["GravityProbeB"], "physics"),
 "pulsar": ("Binary pulsar", ["BinaryPulsar"], "physics"),
 "info": ("Information theory", ["InformationTheory"], "physics"),
 "qnoise": ("Quantum noise", ["QuantumNoise"], "physics"),
 "entropy": ("Entropy and noise", ["EntropyNoise"], "physics"),
 "leeyang": ("Lee–Yang zeros", ["LeeYang", "LeeYangCluster", "LeeYangSpill"], "physics"),
 "fluids": ("Fluids (Navier–Stokes)", ["NavierStokes", "NSContinuum"], "physics"),
 "soliton": ("Solitons (KdV)", ["SolitonsRungs"], "physics"),
 "dark": ("Dark sector (MOND, redshift)", ["DarkSector"], "physics"),
 "cq": ("Quantum chemistry", ["Chemistry.ChemQuantum"], "chemistry"),
 "ccons": ("Conservation of mass and charge", ["Chemistry.ChemConservation"], "chemistry"),
 "cper": ("Periodic law", ["Chemistry.ChemPeriodic"], "chemistry"),
 "cbond": ("Bonding and valence", ["Chemistry.ChemBonding"], "chemistry"),
 "cstoi": ("Stoichiometry", ["Chemistry.ChemStoichiometry"], "chemistry"),
 "cgas": ("Gas laws", ["Chemistry.ChemGas"], "chemistry"),
 "cgibbs": ("Gibbs energy", ["Chemistry.ChemGibbs"], "chemistry"),
 "ceq": ("Equilibrium, mass action", ["Chemistry.ChemEquilibrium"], "chemistry"),
 "crate": ("Rates, Arrhenius", ["Chemistry.ChemRates"], "chemistry"),
 "celec": ("Electrochemistry, Nernst", ["Chemistry.ChemElectro"], "chemistry"),
 "cacid": ("Acids and bases", ["Chemistry.ChemAcids"], "chemistry"),
 "corg": ("Organic chemistry", ["Chemistry.ChemOrganic"], "chemistry"),
 "cbio": ("Biochemistry", ["Chemistry.ChemBiochem"], "chemistry"),
 "bmol": ("Molecular biology", ["Biology.MolecularBiology"], "biology"),
 "bmet": ("Metabolism", ["Biology.Metabolism"], "biology"),
 "bcell": ("Cell biology", ["Biology.CellBiology"], "biology"),
 "bgen": ("Genetics", ["Biology.Genetics"], "biology"),
 "bdev": ("Development", ["Biology.Development"], "biology"),
 "bphys": ("Physiology", ["Biology.Physiology"], "biology"),
 "bneuro": ("Neuroscience", ["Biology.NeuroScience"], "biology"),
 "bevo": ("Evolution", ["Biology.Evolution"], "biology"),
 "bpop": ("Population growth", ["Biology.Population"], "biology"),
 "beco": ("Ecology", ["Biology.Ecology"], "biology"),
 "bbeh": ("Behaviour", ["Biology.Behaviour"], "biology"),
}
T, Rd, L = "textbook", "reading", "lean"
E = [  # (needed, user, tag, why)
 ("prod", "thermo", T, "gas laws and mechanics underlie kinetic theory"),
 ("prod", "em", T, "Coulomb's inverse square and Ohm's law come first"),
 ("em", "sr", T, "light speed is the same in every frame (Maxwell), Einstein 1905"),
 ("prod", "sr", T, "Newton's mechanics is the slow limit SR must reproduce"),
 ("sr", "thomas", L, "two boosts in different directions; the library imports Relativity"),
 ("sr", "gr", T, "GR makes SR's boosts local"),
 ("prod", "gr", T, "Newton's gravity is GR's weak-field limit"),
 ("gr", "gpb", T, "geodetic and frame-dragging turns"),
 ("thomas", "gpb", T, "the velocity twist is Thomas precession"),
 ("gr", "pulsar", T, "spin precession and orbital decay"),
 ("thomas", "pulsar", T, "the velocity twist"),
 ("thermo", "info", Rd, "Shannon's entropy borrows Boltzmann's form (historically; logically independent)"),
 ("info", "qnoise", T, "channels and their capacity"),
 ("qnoise", "entropy", L, "the library imports QuantumNoise"),
 ("info", "entropy", L, "the library imports InformationTheory"),
 ("thermo", "entropy", T, "thermodynamic entropy"),
 ("thermo", "leeyang", T, "partition functions of statistical mechanics"),
 ("prod", "fluids", T, "Newton's second law for a fluid parcel"),
 ("thermo", "fluids", T, "energy equation, viscosity"),
 ("fluids", "soliton", T, "KdV is shallow water"),
 ("prod", "dark", T, "MOND modifies Newton's gravity"),
 ("gr", "dark", Rd, "cosmological redshift composes by the circle operation"),
 ("em", "cq", T, "the atom's Coulomb potential"),
 ("em", "ccons", Rd, "charge conservation is Maxwell's continuity equation"),
 ("cq", "cper", T, "shells fill in order"),
 ("cq", "cbond", T, "bonds are shared orbitals"),
 ("ccons", "cstoi", T, "balancing equations"),
 ("prod", "cgas", T, "Boyle, Charles"),
 ("thermo", "cgas", T, "ideal gas from kinetic theory"),
 ("thermo", "cgibbs", T, "G = H − TS"),
 ("cgibbs", "ceq", T, "ΔG° = −RT ln K"),
 ("cstoi", "ceq", T, "mass action counts by the equation's coefficients"),
 ("thermo", "crate", T, "the Boltzmann factor in Arrhenius"),
 ("cgibbs", "celec", T, "ΔG = −nFE"),
 ("em", "celec", T, "potential and charge"),
 ("ceq", "celec", T, "Nernst has the reaction quotient"),
 ("ceq", "cacid", T, "Ka is an equilibrium constant"),
 ("cbond", "corg", T, "carbon's four bonds"),
 ("corg", "cbio", T, "biomolecules are organic"),
 ("cgibbs", "cbio", T, "coupled reactions"),
 ("crate", "cbio", T, "enzyme kinetics"),
 ("cbio", "bmol", T, "DNA, RNA, proteins"),
 ("cbio", "bmet", T, "pathways"),
 ("cgibbs", "bmet", T, "ATP's free energy"),
 ("bmol", "bcell", T, "the cell's machinery"),
 ("bmol", "bgen", T, "genes are DNA"),
 ("bgen", "bdev", T, "gene regulation"),
 ("bcell", "bdev", T, "cell division and signalling"),
 ("bcell", "bphys", T, "tissues"),
 ("bmet", "bphys", T, "energy budgets"),
 ("bphys", "bneuro", T, "neurons are cells"),
 ("celec", "bneuro", T, "Nernst and Goldman potentials"),
 ("bgen", "bevo", T, "heredity"),
 ("bpop", "bevo", T, "Fisher–Haldane–Wright population genetics"),
 ("bpop", "beco", T, "Lotka–Volterra builds on logistic growth"),
 ("bevo", "bbeh", T, "behavioural ecology, kin selection"),
 ("bneuro", "bbeh", T, "behaviour runs on nervous systems"),
]
lines = []
def say(*a):
    s = " ".join(map(str, a)); print(s); lines.append(s)
shape = {}
for k, (lab, ms, f) in N.items():
    ps = [p for m in ms for p in mods.get(m, [])]
    b = collections.Counter(p["body"] for p in ps); c = collections.Counter(p["cluster"] for p in ps)
    bb, bn = b.most_common(1)[0]; cc, cn = c.most_common(1)[0]
    shape[k] = dict(label=lab, field=f, n=len(ps), body=bb, body_share=bn / len(ps), family=cc,
                    family_share=cn / len(ps), bodies=b.most_common())
    say(f"{lab:<44} n={len(ps):>3} body {bb:<6} {100*bn/len(ps):>3.0f}%  family {cc:<14} {100*cn/len(ps):>3.0f}%")
BO = {"point": 0, "rod": 1, "rect": 2, "box": 3, "hyper": 3, "cube": 4}  # box and hyperbox one layer
def flow(lab):
    up = dn = eq = 0
    for a, b, *_ in E:
        x, y = BO[lab[a]], BO[lab[b]]
        if x < y: up += 1
        elif x > y: dn += 1
        else: eq += 1
    return up, dn, eq
lab = {k: v["body"] for k, v in shape.items()}
up, dn, eq = flow(lab)
random.seed(3); ks = list(lab); vals = [lab[k] for k in ks]; null = []
for _ in range(20000):
    random.shuffle(vals); u2, d2, _e = flow(dict(zip(ks, vals))); null.append(u2 - dn if False else u2 - d2)
obs = up - dn; mu = sum(null) / len(null); sd = (sum((v - mu) ** 2 for v in null) / len(null)) ** .5
p = sum(1 for v in null if v >= obs) / len(null)
say(f"\nBody flow along textbook edges ({len(E)}): up {up}, down {dn}, same layer {eq}; net {obs};"
    f" null (bodies shuffled across laws) mean {mu:.2f}, sd {sd:.2f}, z {(obs-mu)/sd:.2f}, p {p:.3f}")
# collapse: body graph and family graph
BG = collections.Counter((lab[a], lab[b]) for a, b, *_ in E if lab[a] != lab[b])
say("Body-level edges (needed body -> user body: count):", dict(BG))
FG = collections.Counter((shape[a]["family"], shape[b]["family"]) for a, b, *_ in E if shape[a]["family"] != shape[b]["family"])
two_way = sorted({tuple(sorted(k)) for k in FG if (k[1], k[0]) in FG})
say("Family-level two-way pairs (cycles of length 2):", two_way)
# body two-way
say("Body-level two-way pairs:", sorted({tuple(sorted(k)) for k in BG if (k[1], k[0]) in BG}))
# exceptions
say("\nEdges that run against the body order (a law built on a law of bigger body):")
ex = []
for a, b, t, why in E:
    if BO[lab[a]] > BO[lab[b]]:
        say(f"  {shape[a]['label']} ({lab[a]}) -> {shape[b]['label']} ({lab[b]})  [{t}]")
        ex.append(dict(a=a, b=b))
# chains: longest textbook path
import functools
succ = collections.defaultdict(list)
for a, b, *_ in E: succ[a].append(b)
@functools.lru_cache(None)
def longest(v):
    best = [v]
    for w in succ[v]:
        c = [v] + list(longest(w))
        if len(c) > len(best): best = c
    return tuple(best)
roots = [k for k in N if not any(b == k for _, b, *_ in E)]
say("\nRoots (need nothing in this table):", [shape[r]["label"] for r in roots])
chain = max((longest(r) for r in roots), key=len)
say("Longest chain:", " -> ".join(f"{shape[k]['label']} [{lab[k]}]" for k in chain))
json.dump(dict(nodes=shape, edges=[dict(a=a, b=b, tag=t, why=w) for a, b, t, w in E],
               flow=dict(up=up, down=dn, same=eq, null_mean=mu, null_sd=sd, p=p),
               body_edges={f"{a}>{b}": w for (a, b), w in BG.items()}, against=ex,
               chain=list(chain), roots=roots), open(HERE / "laws.json", "w"), indent=1, ensure_ascii=False)
open(HERE / "laws_output.txt", "w").write("\n".join(lines) + "\n")
