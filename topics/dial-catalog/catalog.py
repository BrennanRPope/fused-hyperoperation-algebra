"""The law groups every dial survey sorts (Nov 26, 2026).
Groups 1-49 are the shape-dial survey's list (topics/shape-dial-survey/laws.py), copied, not changed.
Groups 50-65 cover the Lean modules that list left out (checked against lean-library/theorems.tsv).
Each entry: (law group, Lean modules, the shape it lives on)."""
GROUPS = [
 ('Core law o_k: order free, identity 0, the wall', 'Base, Foundations', 'cube'),
 ('Corner expansion: 1 + (x1 o ... o xn) = sum over corners', 'Base, Cube', 'cube'),
 ('Powers: 1 + a^(n) = (1 + a)^n', 'Base', 'cube'),
 ('Rungs, height flow, Abel ruler', 'Rungs, FractionalRungs, RungShape, ZoomDial', 'half-line'),
 ('Calculus: fundamental theorem with o-steps', 'Calculus, ClosingRate', 'half-line'),
 ('Plane law (k, d): circle, parabola, hyperbola', 'Landscape, DialSweep, DialWall, DialSurvey', 'conics'),
 ('Clifford signatures', 'Clifford', 'conics per axis'),
 ('Elliptic curves: pinching the torus', 'Elliptic, EllipticS3, EllipticCounting', 'torus'),
 ('Einstein velocity addition', 'Relativity, Electromagnetism', 'ball'),
 ("Bayes' rule in odds form", 'Statistics', 'simplex'),
 ('Swap rule: turn = d x area (Thomas, Foucault, pulsars)', 'TwistedBoosts, GRBridges, GravityProbeB, BinaryPulsar', 'sphere / saddle'),
 ('Maxwell as one equation, F as a wedge', 'Electromagnetism', 'light cone'),
 ("Einstein tensor, Lovelock's step", 'GeneralRelativity, GRBridges', 'curved space'),
 ('Curvature by counting (2D causal diamonds)', 'curved-space page', 'cube order (triangle light cone)'),
 ('Corner-weight shapes, power-set dials', 'PowersetDials, FractionalPowerset, FractalDials', 'generalized permutohedra'),
 ('Bell numbers, Catalan, set partitions', 'Combinatorics', 'associahedron / permutohedron'),
 ('q-counting [n]_q', 'QAcross.QCounting, QAcross.QDial', 'simplex'),
 ('Cube face weights q^area, solitons', 'SolitonsRungs, QAcross.Partitions', 'cube faces / Grassmannian'),
 ('Toric charts, projective plane', 'Toric.ToricCharts, Toric.ProjectivePlane', 'simplex'),
 ('Knots, graphs, diamond lattice, patterns', 'KnotTheory, GraphTheory, DiamondLattice, Patterns', 'various'),
 ('Union rule o_-1, independence', 'MeasureProbability, Foundations.Independence', 'cube'),
 ('Logistic (independent yes/no) and softmax (one of many)', 'Statistics', 'cube / simplex'),
 ('Renyi = L_k(Tsallis)', 'InformationTheory', 'simplex'),
 ("Shannon, Campbell's code lengths", 'InformationFull', 'simplex'),
 ('Repetition dial: Fermi, classical, Bose', 'RepetitionDial', 'cube / orthant'),
 ('Thermodynamics: Gibbs, temperature as Maslov dial', 'Thermodynamics, EntropyNoise', 'free-energy bodies'),
 ('Lee-Yang circle theorem', 'LeeYang, LeeYangCluster, LeeYangSpill', 'circle'),
 ('Born rule (quantum noise)', 'QuantumNoise', 'sphere'),
 ('Quantum correlations: CHSH <= 2 sqrt 2', 'beyond-the-cube page', 'elliptope'),
 ('Entropy of regions as a shape, holographic rules', 'exploration, polar-recast pages', 'polymatroid'),
 ("Game theory and evolution: replicator, Fisher's theorem", 'GameTheory, Biology.Evolution', 'simplex'),
 ('Consensus, dial point processes', 'Consensus, DialPointProcesses', 'various'),
 ('Navier-Stokes rungs, strain', 'NavierStokes, NSContinuum', 'strain triangle'),
 ('Dark sector', 'DarkSector', '-'),
 ("Newton's G lead", 'newton-g page', '-'),
 ('Fermat n = 3, n = 4, Pythagorean triples', 'Fermat.Fermat3, Fermat.Fermat4', 'p-circle with p = n'),
 ('Frey pieces, counterexample map, FLT route', 'Fermat.FreyPieces, Fermat.CounterexampleMap, Fermat.FermatRoute', 'torus of angles'),
 ('Hopf invariant one, Kummer-Wilson', 'Topology', '-'),
 ("P vs NP: what can't work", 'PvsNP', '-'),
 ('e: corner volumes 1/n!', 'E, FusedConstants', 'simplex'),
 ('pi and the Gaussian', 'Pi, GeometryTrig', 'sphere'),
 ("i, Euler's identity, complex analysis", 'EIPi, WhyI, ComplexAnalysis', 'sphere'),
 ('Conservation of mass and charge, stoichiometry', 'Chemistry.ChemConservation, ChemStoichiometry', 'simplex'),
 ('Equilibrium, Nernst, Arrhenius, pH, gas laws', 'Chemistry (other modules)', 'log axis'),
 ('Biology: genetics, physiology, ecology, cells', 'Biology (other modules)', 'various'),
 ("Finance: Fisher's rule for real rates", 'Finance', 'log axis'),
 ('Logic, sets, categories, algebra, linear algebra', 'LogicSets, CategoryTheory, AbstractAlgebra, LinearAlgebra', '-'),
 ('Functional analysis, Osgood inverse, numerics, dynamics', 'FunctionalAnalysis, OsgoodInverse, Numerical, Dynamics, Infinities', '-'),
 ('Data centres, Rulecraft, rule breakers, early science', 'DataCenter, Rulecraft, RuleBreakers, EarlyScience', '-'),
 ('Hyperoperation ladder: whole rungs fused', 'HyperLadder', 'half-line'),
 ('Lie groups: the plane law as a Mobius product', 'GroupsLie', 'conics'),
 ('Number planes a + b u with u^2 = k u - q; quaternions, precession', 'NumberPlane, DiffGeometry, Hypersphere, Precession', 'conics / sphere'),
 ('Jet fields R[e]/(e^2 - k e); jet river: Schwarzschild, Kerr, Kerr-Newman', 'JetFields, JetRiver, SpinningRiver', 'curved space'),
 ('Order dial: the swap gap is an oriented area; Weyl moves', 'OrderDial, CutIndependence, QuantumAction', 'plane of moves'),
 ('Units and roots: the weight is units, every law is addition plus a corner', 'ShapeDirections, ShapeDependencies, Conventions, PathNames', '-'),
 ('Dial potential, several waves, Dirac checkerboard', 'DialPotential, SeveralWaves, Checkerboard, CheckerboardPotential', 'cube of insertions'),
 ('Size weight, fractional boxes, corner weights', 'SizeDialResidual, FractalBoxes, CornerWeights', 'box'),
 ('Holographic gluing: monogamy from the odds', 'HoloGluing, HoloGluingField', 'cube'),
 ('Force reach, the hidden-circle force, flux-tube dial, Yang-Mills gap', 'ForceReach, NewForce, FluxTubeDial, YangMillsGap', 'conics'),
 ('Tunnelling: barriers compose by o_k', 'Tunnelling', 'conics'),
 ('Riemann rounds: prime cube fold, fused pairs', 'Riemann, RiemannElliptic, RiemannRound3, RiemannRound4', 'cube / conics'),
 ('Elliptic modulus: kappa = 1 - 2m', 'EllipticModulus', 'torus'),
 ('Beyond one loop: plane-law flows and their walls', 'BeyondOneLoop', 'conics'),
 ('Boolean lattice, fresh look, fuzzy frame', 'BooleanLattice, FreshLook, FuzzyFrame', 'cube'),
 ('Warp bubbles and wormholes: energy conditions', 'WarpWormhole', 'strain triangle'),
]

def check_coverage(tsv="../../lean-library/theorems.tsv"):
    """Return the Lean module families in theorems.tsv that no group names."""
    import os
    here = os.path.dirname(os.path.abspath(__file__))
    covered = set()
    for _, mods, _ in GROUPS:
        for m in mods.replace(" ", "").split(","):
            covered.add(m.split(".")[0])
    missing = set()
    for line in open(os.path.join(here, tsv)):
        parts = line.split("\t")[0].split(".")
        key = parts[2] if len(parts) > 2 and parts[1] == "Topics" else parts[1]
        if key not in covered:
            missing.add(key)
    return sorted(missing)

if __name__ == "__main__":
    print(len(GROUPS), "groups; Lean families not named:", check_coverage())
