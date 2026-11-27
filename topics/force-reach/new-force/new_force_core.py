"""Shared inputs for the new-force checks, computed here so that no check script has to read another's output.
Each check script writes only its own results_<part>.json (atomically); figures.py merges them into results.json.
This lets CI run the scripts in any order or in parallel."""
import glob, json, math, os, tempfile
import mpmath as mp

lP = 1.616255e-35                       # Planck length, m
MP = 1.220890e19                        # Planck mass, GeV (non-reduced)
Mbar = MP/math.sqrt(8*math.pi)          # reduced Planck mass, GeV
ainv = 137.035999084
R_lP = 2*math.sqrt(ainv)                # part G: the electric circle's radius, Planck lengths
LAB_RANGE_M = 38.6e-6                   # Eot-Wash 2020, gravitational-strength Yukawa range limit

def range_lP(NF=2, Rl=R_lP):
    """Part H: the size carrier's range in Planck lengths for NF periodic 5D Dirac fermions (graviton 5 dof)."""
    Cv = 3*float(mp.zeta(5))*(4*NF - 5)/(64*math.pi**6)
    m = math.sqrt(20*Cv)/(Rl**2/math.sqrt(8*math.pi))
    return 1/m, Cv

def radion_GeV(NF=2):
    return MP/range_lP(NF)[0]

def copy1_GeV():
    return MP/R_lP

def lab_over_range(NF=2):
    return LAB_RANGE_M/(range_lP(NF)[0]*lP)

def inflation_bounds(As=2.1e-9, rmax=0.036):
    """Part I: largest inflationary Hubble rate and energy scale allowed by r < rmax (GeV)."""
    H = math.pi*Mbar*math.sqrt(As*rmax/2)
    return H, (3*H**2*Mbar**2)**0.25

def write_json_atomic(path, data):
    d = os.path.dirname(os.path.abspath(path))
    fd, tmp = tempfile.mkstemp(dir=d, prefix=".tmp_", suffix=".json")
    with os.fdopen(fd, "w") as f:
        json.dump(data, f, indent=1, default=str)
    os.replace(tmp, path)

def merged_results(write=True):
    """Merge every results_<part>.json into one dict (and into results.json when write is True)."""
    out = {}
    for p in sorted(glob.glob("results_*.json")):
        out.update(json.load(open(p)))
    if write:
        write_json_atomic("results.json", out)
    return out
