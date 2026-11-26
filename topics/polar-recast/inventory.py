"""Group the compiled Lean theorems (lean-library/theorems.tsv) into families, fundamental first."""
import collections
FAM = [
 ("1 Core law and its algebra", ["Base","Foundations.","Cube","AbstractAlgebra","LinearAlgebra","CategoryTheory","LogicSets","Conventions","PathNames"]),
 ("2 Rungs, height and calculus", ["Rungs","FractionalRungs","RungShape","Topics.E","OsgoodInverse","Dynamics","Calculus","Numerical","FunctionalAnalysis","JetFields","ZoomDial","Infinities","ClosingRate"]),
 ("3 Plane dial and composition shapes", ["Landscape","ShapeDirections","NumberPlane","GeometryTrig","Elliptic","EllipticS3","EIPi","WhyI","Pi","Clifford","ComplexAnalysis","GroupsLie","DiffGeometry","DialWall","DialSweep","DialSurvey","DialPotential","OrderDial","Relativity","FusedConstants"]),
 ("4 Swap rule: precession and curvature", ["TwistedBoosts","GravityProbeB","BinaryPulsar","GeneralRelativity","GRBridges","Electromagnetism"]),
 ("5 Corners, shapes and counting", ["PowersetDials","FractionalPowerset","Combinatorics","DiamondLattice","GraphTheory","KnotTheory","QAcross.","Toric.","Topology.CubeProofs","Patterns","Checkerboard","CheckerboardPotential","FractalBoxes","FractalDials","SizeDialResidual","SolitonsRungs","SeveralWaves"]),
 ("6 Probability, information, statistics, heat", ["MeasureProbability","Statistics","InformationTheory","InformationFull","EntropyNoise","QuantumNoise","DialPointProcesses","LeeYang","LeeYangCluster","LeeYangSpill","Consensus","GameTheory","RepetitionDial","Thermodynamics"]),
 ("7 Fields and flows", ["NavierStokes","NSContinuum","DarkSector"]),
 ("8 Number theory and complexity", ["Fermat.","KummerWilson","HopfInvariantOne","EllipticCounting","PvsNP"]),
 ("9 Applied sciences", ["Chemistry.","Biology.","EarlyScience","Finance","DataCenter","Rulecraft","RuleBreakers"]),
]
def fam(mod):
    m = mod.replace("FusedAlgebra.","")
    m2 = m.replace("Topics.","")
    best=None
    for name, keys in FAM:
        for k in keys:
            if (k.endswith(".") and (m.startswith(k) or m2.startswith(k) or ("."+k) in "."+m)) or m2==k or m==k or (k=="Topics.E" and m=="Topics.E") or (k=="Topology.CubeProofs" and m2==k):
                if best is None or len(k)>len(best[1]): best=(name,k)
    if m2.startswith("Topology.") and best is None: best=("8 Number theory and complexity","")
    return best[0] if best else "unsorted"
cnt=collections.Counter(); mods=collections.defaultdict(set)
for line in open("../lean-library/theorems.tsv"):
    mod=line.split("\t")[0]; f=fam(mod); cnt[f]+=1; mods[f].add(mod.replace("FusedAlgebra.","").replace("Topics.",""))
tot=sum(cnt.values())
with open("inventory.md","w") as out:
    out.write(f"# Lean theorem families (compiled snapshot, {tot} theorems)\n\n| Family | Theorems | Modules |\n|---|---|---|\n")
    for name,_ in FAM+[("unsorted",[])]:
        if cnt[name]: out.write(f"| {name} | {cnt[name]} | {', '.join(sorted(mods[name]))} |\n")
for name,_ in FAM+[("unsorted",[])]:
    if cnt[name]: print(name, cnt[name], len(mods[name]))
print("total", tot)
