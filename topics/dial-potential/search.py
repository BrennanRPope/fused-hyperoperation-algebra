"""Search for the dial's potential and what it couples to: 20 guesses a round, keep the closest, refine.

Run:  python3 -I search.py      (writes results.json and search_output.txt; about two minutes)

Question (Brennan, Nov 25 2026): if the circle-operation dial k (a + b + k·ab) behaves like the Higgs field, with a
hump at k = 0 and dips at ±k0, what is its potential and what does it couple to?  Make 20 guesses, score
them against measured numbers, keep the closest, iterate.

SCORING, in bits (log2 of a likelihood ratio), so every guess is judged on its own evidence:
  * A prediction p for a measured value x with total log-uncertainty s (experiment plus the guess's own
    theory error) earns   bits = log2( W / (s·sqrt(2π)) ) - z^2 / (2 ln 2),   z = ln(p/x)/s,
    against a prior that spreads the value evenly in log over a width W (7 e-folds for O(1) ratios, i.e.
    0.01 to 10; 92 e-folds for hierarchy numbers).  A close, sharp prediction scores high; a miss
    goes very negative.  |z| > 3 is marked ruled out.
  * A target used to fix a free parameter scores nothing (that is the parameter's cost).
  * Look-elsewhere: picking a constant out of a list of G candidates costs log2(G) bits, and the whole
    search costs log2(number of guesses tried).  "Net bits" subtracts both.
  * Cosmology guesses are scored by chi^2 against Pantheon+ and DESI DR2 BAO relative to a cosmological
    constant, bits = Δchi^2/(2 ln 2), with free parameters fitted on one data set and scored on the other.
  * Chance calibration: for the formula families, the same pick-the-closest procedure is run against
    random targets to see how often chance alone does as well.
"""
import json, math, os, sys, itertools, time
import numpy as np
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import sm_rge as R
import quintessence as Q

LN2 = math.log(2)
OUT = []


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s); OUT.append(s)


# ------------------------------------------------------------------ measured numbers (PDG 2024 unless noted)
v = 246.21965                     # (sqrt2 G_F)^-1/2, GeV
MH, sMH = 125.20, 0.11
MW, sMW = 80.3692, 0.0133
MZ, sMZ = 91.1880, 0.0020
MT, sMT = 172.57, 0.29
AS, sAS = 0.1180, 0.0009
ALPHA_INV = 137.035999177
S2W_MS = 0.23129                  # sin^2 theta_W, MS-bar at M_Z
ALPHA_MZ_INV = 127.930
MPL_RED = 2.43536e18
RHO_L14 = 2.24e-12                # dark energy density^(1/4), GeV (Planck 2018: Omega_L 0.685, H0 67.4)

TARGETS = {
    #  name          value               rel sigma (exp)        prior width (e-folds)  note
    "mH/v":      (MH / v,             sMH / MH,              math.log(1e3), "Higgs mass over the vev"),
    "mW/v":      (MW / v,             sMW / MW,              math.log(1e3), "W mass over the vev"),
    "mZ/v":      (MZ / v,             sMZ / MZ,              math.log(1e3), "Z mass over the vev"),
    "mt/v":      (MT / v,             sMT / MT,              math.log(1e3), "top mass over the vev"),
    "sin2W":     (S2W_MS,             0.00004 / S2W_MS,      math.log(1e3), "weak mixing angle (MS-bar, M_Z)"),
    "alpha_s":   (AS,                 sAS / AS,              math.log(1e3), "strong coupling at M_Z"),
    "alpha":     (1 / ALPHA_INV,      1.5e-10,               math.log(1e3), "fine-structure constant"),
    "v/MPl":     (v / MPL_RED,        0.0,                   math.log(1e40), "weak scale over reduced Planck mass"),
    "rhoL/MPl":  (RHO_L14 / MPL_RED,  0.01,                  math.log(1e40), "dark energy scale over reduced Planck mass"),
}
# linear-scale targets
GAMMA, sGAMMA = 1 + 2.1e-5, 2.3e-5          # PPN gamma, Cassini (Bertotti et al. 2003)
DALPHA, sDALPHA = 1.3e-6, 1.4e-6            # Δalpha/alpha at z = 1.15, ESPRESSO (Murphy et al. 2022)


def bits_log(name, pred, th_rel):
    """Bits for a log-scale prediction with theory error th_rel (relative)."""
    x, se, W, _ = TARGETS[name]
    if pred is None or not np.isfinite(pred) or pred <= 0:
        return -math.inf, math.inf
    s = math.sqrt(se ** 2 + th_rel ** 2)
    z = math.log(pred / x) / s
    return math.log2(W / (s * math.sqrt(2 * math.pi))) - z * z / (2 * LN2), z


def bits_lin(pred, x, s, W):
    z = (pred - x) / s
    return math.log2(W / (s * math.sqrt(2 * math.pi))) - z * z / (2 * LN2), z


TREE = 0.02          # theory error assigned to a tree-level relation between pole masses (loop corrections ~1-3%)

# ------------------------------------------------------------------ guess machinery
class Guess:
    def __init__(self, gid, family, name, shape, coupling, geometry, grade, n_free=0, choice_bits=0.0, meta=None):
        self.gid, self.family, self.name = gid, family, name
        self.shape, self.coupling, self.geometry, self.grade = shape, coupling, geometry, grade
        self.n_free, self.choice_bits, self.meta = n_free, choice_bits, meta or {}
        self.preds = {}          # target -> (pred, theory_rel)  or ('lin', pred, x, s, W)
        self.fitted = []
        self.extra = {}

    def score(self):
        tot, worst, rows = 0.0, 0.0, []
        for t, p in self.preds.items():
            if p[0] == "lin":
                b, z = bits_lin(*p[1:])
                val = p[1]
            else:
                b, z = bits_log(t, p[0], p[1])
                val = p[0]
            rows.append({"target": t, "pred": val, "z": z, "bits": b})
            tot += b
            worst = max(worst, abs(z))
        if "cosmo_bits" in self.extra:
            tot += self.extra["cosmo_bits"]
            worst = max(worst, self.extra.get("cosmo_z", 0))
        self.bits, self.rows = tot, rows
        self.ruled_out = worst > 3 or tot == -math.inf
        self.max_z = worst
        return tot

    def as_dict(self):
        return {"id": self.gid, "family": self.family, "name": self.name, "shape": self.shape,
                "coupling": self.coupling, "geometry": self.geometry, "grade": self.grade,
                "free_params": self.n_free, "fitted_on": self.fitted, "choice_bits": self.choice_bits,
                "bits": self.bits, "net_bits": getattr(self, "net", None), "ruled_out": self.ruled_out,
                "max_abs_z": self.max_z, "rows": self.rows, "extra": self.extra, "meta": self.meta}


# ------------------------------------------------------------------ grammars (for honest look-elsewhere counts)
def grammar_constants(maxpq=8, pis=(-2, -1, 0, 1, 2)):
    """Numbers p/q · pi^a, small p, q.  Each is one 'geometric-looking' normalisation."""
    seen = {}
    for p in range(1, maxpq + 1):
        for q in range(1, maxpq + 1):
            if math.gcd(p, q) != 1:
                continue
            for a in pis:
                val = p / q * math.pi ** a
                lab = f"{p}/{q}" + ("" if a == 0 else f"·π^{a}")
                seen.setdefault(round(val, 12), lab)
    return sorted((v_, l) for v_, l in seen.items())


GRAM = grammar_constants()
G_BITS = math.log2(len(GRAM))

SUM_RULES = []      # c_W, c_Z, c_t in {0,1,2} for m_H^2 + cW mW^2 + cZ mZ^2 + ct mt^2 = cv v^2  (or with sign)
for cW, cZ, ct, sgn in itertools.product((0, 1, 2), (0, 1, 2), (0, 1, 2, 3, 4), (1, -1)):
    for cv in (Fraction(1, 4), Fraction(1, 2), Fraction(1), Fraction(2), Fraction(0)):
        SUM_RULES.append((cW, cZ, ct, sgn, cv))
SR_BITS = math.log2(len(SUM_RULES))


def sumrule_mh(rule):
    cW, cZ, ct, sgn, cv = rule
    # m_H^2 = cv v^2 - sgn (cW mW^2 + cZ mZ^2 + ct mt^2)
    m2 = float(cv) * v * v - sgn * (cW * MW ** 2 + cZ * MZ ** 2 + ct * MT ** 2)
    return math.sqrt(m2) if m2 > 0 else float("nan")


def sumrule_label(rule):
    cW, cZ, ct, sgn, cv = rule
    terms = [f"{c if c > 1 else ''}m{n}²" for c, n in ((cW, "W"), (cZ, "Z"), (ct, "t")) if c]
    side = (" + " if sgn > 0 else " − ") + " + ".join(terms) if terms else ""
    if sgn < 0 and terms:
        side = " − (" + " + ".join(terms) + ")"
    return f"mH²{side} = {cv}·v²"


# ------------------------------------------------------------------ physics helpers
def flat_dial_prediction(scale="MPl", with_beta=False):
    """m_H (and m_t) from 'the dial is flat at the top scale'.  Returns dict with predictions and errors."""
    if scale == "MPl" and not with_beta:
        mh = R.mh_for_flat(MT, AS)
        dm_t = 0.0066 * sMT / 0.0029
        dm_a = 0.0018 * (sAS / 0.0007) / 0.0029
        th = math.sqrt(dm_t ** 2 + dm_a ** 2 + (0.0003 / 0.0029) ** 2)   # Buttazzo theory error ~3e-4 in lambda
        mh2 = R.mh_for_flat(MT, AS, mu=R.MPL_FULL, use_buttazzo=False)    # our two-loop run, cross-check
        return {"mh": mh, "mh_err": th, "mh_own_rge": mh2}
    mu = {"MPl": R.MPL_FULL, "MPl_red": R.MPL_RED, "GUT": 2e16, "string": 1e17, "seesaw": 1e14}[scale]
    if not with_beta:
        mh = R.mh_for_flat(MT, AS, mu=mu, use_buttazzo=False, offset=0.0002)   # +0.0002 = our offset vs NNLO
        # error: propagate top mass by finite difference
        mh_p = R.mh_for_flat(MT + sMT, AS, mu=mu, use_buttazzo=False, offset=0.0002)
        mh_a = R.mh_for_flat(MT, AS + sAS, mu=mu, use_buttazzo=False, offset=0.0002)
        err = math.sqrt((mh_p - mh) ** 2 + (mh_a - mh) ** 2 + 0.4 ** 2)
        return {"mh": mh, "mh_err": err}
    # lambda(mu) = 0 and beta_lambda(mu) = 0: solve for (m_H, m_t)
    from scipy.optimize import fsolve
    def eqs(p):
        l, b = R.lambda_at(p[0], p[1], AS, mu)
        return [(l + 0.0002) * 100, b * 1e4]
    sol, info, ier, msg = fsolve(eqs, [128.0, 171.0], full_output=True)
    mh, mt = sol
    solA = fsolve(eqs, [128.0, 171.0], full_output=False)
    # theory error: vary alpha_s
    def eqs_a(p):
        l, b = R.lambda_at(p[0], p[1], AS + sAS, mu)
        return [(l + 0.0002) * 100, b * 1e4]
    sa = fsolve(eqs_a, sol)
    return {"mh": mh, "mt": mt, "mh_err": math.hypot(sa[0] - mh, 0.6), "mt_err": math.hypot(sa[1] - mt, 0.4),
            "ier": int(ier)}


def gw_cw_mass2():
    """Gildener-Weinberg: m_h^2 = (6 mW^4 + 3 mZ^4 - 12 mt^4)/(8 pi^2 v^2) for a Higgs that is a pure dilaton."""
    return (6 * MW ** 4 + 3 * MZ ** 4 - 12 * MT ** 4) / (8 * math.pi ** 2 * v * v)


def unification_s2w():
    """Predict sin^2 theta_W(M_Z) from alpha_em(M_Z), alpha_s(M_Z) if the three couplings meet (SM, 1 loop)."""
    from scipy.optimize import fsolve
    b = np.array([41 / 10, -19 / 6, -7.0])
    def eqs(p):
        s2, L = p
        inv = np.array([0.6 * (1 - s2) * ALPHA_MZ_INV, s2 * ALPHA_MZ_INV, 1 / AS]) - b * L / (2 * math.pi)
        return [inv[0] - inv[1], inv[1] - inv[2]]
    s2, L = fsolve(eqs, [0.21, 30.0])
    return s2, MZ * math.exp(L)


# cosmology fits are cached because they are the slow part
_CF = {}
PHI0 = {"exp": 0.0, "hilltop": 0.1, "cosine": 1.0}   # starting field: hilltop 0.1 M_Pl off the top; cosine 1 rad... see note


def cosmo_fit(shape, par, which):
    key = (shape, None if par is None else round(par, 6), which)
    if key not in _CF:
        bounds = {"exp": (0.0, 2.5), "hilltop": (0.4, 4.0), "cosine": (0.3, 3.0)}[shape]
        if shape != "exp" and par is None:
            _CF[key] = Q.fit_grid(shape, which, np.linspace(bounds[0], bounds[1], 12), PHI0[shape])
        else:
            _CF[key] = Q.fit(shape, par, which, par_bounds=bounds, phi0=PHI0[shape])
    return _CF[key]


LCDM = {}


def lcdm(which):
    if which not in LCDM:
        LCDM[which] = Q.fit("exp", 0.0, which)
    return LCDM[which]


def cosmo_score(g, shape, par):
    """Held-out score: fit on SN, score BAO; fit on BAO, score SN; plus full-data fit (reported)."""
    out = {}
    for train, test in (("sn", "bao"), ("bao", "sn")):
        f = cosmo_fit(shape, par, train)
        l = lcdm(train)
        # recompute the held-out chi2 at the trained parameters
        E, w = Q.history(f["Om"], shape, f["par"], PHI0[shape])
        sn, bao = Q.chi2s(E)
        El, _ = Q.history(l["Om"], "exp", 0.0)
        snl, baol = Q.chi2s(El)
        test_g, test_l = (bao, baol) if test == "bao" else (sn, snl)
        out[f"train_{train}"] = {"par": f["par"], "Om": f["Om"], "heldout_chi2": test_g,
                                 "heldout_chi2_lcdm": test_l, "w0": float(w[0])}
    full = cosmo_fit(shape, par, "both")
    fl = lcdm("both")
    E, w = Q.history(full["Om"], shape, full["par"], PHI0[shape])
    ho = sum(out[k]["heldout_chi2_lcdm"] - out[k]["heldout_chi2"] for k in out)
    out["full"] = {"par": full["par"], "Om": full["Om"], "chi2": full["chi2"], "chi2_lcdm": fl["chi2"],
                   "dchi2": fl["chi2"] - full["chi2"], "w0": float(w[0]),
                   "wa_proxy": float((np.interp(1.0, Q.ZG, w) - w[0]) / 0.5)}
    out["heldout_dchi2_sum"] = ho
    g.extra.update(out)
    g.extra["cosmo_bits"] = ho / (2 * LN2)
    g.extra["cosmo_z"] = 0.0 if ho > -9 else 3.1
    return out


def dalpha_if_light_dial(shape, par, Om, n_alpha=0.5):
    """If the dial is light's dial (k = 1/c^2 and alpha ∝ k^(1/2)), and the canonical field is
    phi = sqrt(3/2) ln k (conformal coupling), then Δalpha/alpha = n_alpha · Δ ln k = n_alpha·Δphi/sqrt(3/2)."""
    E, w = Q.history(Om, shape, par)
    phi = Q.history.phi
    dphi = float(np.interp(1.15, Q.ZG, phi) - phi[0])
    return n_alpha * dphi / math.sqrt(1.5), dphi


# ------------------------------------------------------------------ round 1: twenty hand-made guesses
def round1():
    G = []
    def add(*a, **k):
        g = Guess(*a, **k); G.append(g); return g

    g = add("P1", "tree_lam", "Quartic sombrero, height free", "V = (λ/4)(k² − k0²)²", "k = Higgs magnitude, k0 = v",
            "flip symmetry k → −k (proved) makes V even; a quartic is the simplest even well", "B", n_free=1)
    g.fitted = ["mH/v"]
    g = add("P2", "tree_lam", "Lift-normalised well", "V = (k ∘₁ (−k) + k0²)² = (k0² − k²)²",
            "k = Higgs magnitude", "k ∘₁ (−k) = −k² (proved): the well is the lift of a dial against its mirror, coefficient 1",
            "B", meta={"lam": 4.0})
    g.preds["mH/v"] = (math.sqrt(2 * 4.0), TREE)
    g = add("P3", "tree_lam", "Higgs at half the vev", "V = (1/32)(k² − k0²)²  (λ = 1/8)", "k = Higgs magnitude",
            "none beyond the even well; λ = 1/8 is picked as a small fraction", "C", choice_bits=G_BITS, meta={"lam": 0.125})
    g.preds["mH/v"] = (math.sqrt(2 * 0.125), TREE)
    fd = flat_dial_prediction("MPl")
    g = add("P4", "flat", "Dial flat at the Planck scale", "λ(M_Pl) = 0: the well vanishes at the top scale",
            "k = Higgs magnitude, Standard Model running",
            "at k's symmetric point the circle operation is plain addition; 'flat at the top' = no hump left at M_Pl", "B",
            meta={"scale": "MPl", "beta": False, "own_rge_mh": fd["mh_own_rge"]})
    g.preds["mH/v"] = (fd["mh"] / v, fd["mh_err"] / fd["mh"])
    fdb = flat_dial_prediction("MPl", with_beta=True)
    g = add("P5", "flat", "Dial flat and level at the Planck scale", "λ(M_Pl) = 0 and β_λ(M_Pl) = 0",
            "k = Higgs magnitude", "flat and level: the hump and its slope both vanish at the top scale", "B",
            meta={"scale": "MPl", "beta": True})
    g.preds["mH/v"] = (fdb["mh"] / v, fdb["mh_err"] / fdb["mh"])
    g.preds["mt/v"] = (fdb["mt"] / v, fdb["mt_err"] / fdb["mt"])
    m2 = gw_cw_mass2()
    g = add("P6", "cw", "Scale-symmetric dial (Coleman–Weinberg)", "V = λ k⁴ (ln(k²/k0²) − ½), no tree mass",
            "k = Higgs magnitude", "the dial's rescaling symmetry (proved) forbids a tree-level k0; only quantum running can set it",
            "A", meta={"m2_GeV2": m2})
    g.preds["mH/v"] = (math.sqrt(m2) / v if m2 > 0 else None, TREE)
    vm = math.sqrt(4 * MT ** 2 - 2 * MW ** 2 - MZ ** 2)
    g = add("P7", "veltman", "Dial blind to the cutoff (Veltman)", "quadratic sensitivity cancels: mH² = 4mt² − 2mW² − mZ²",
            "k = Higgs magnitude", "the symmetric point k = 0 is not pulled by heavy physics", "C")
    g.preds["mH/v"] = (vm / v, TREE)
    rule = (1, 1, 1, 1, Fraction(1))
    g = add("P8", "sumrule", "Mass budget sum rule", "mH² + mW² + mZ² + mt² = v²", "k sets every mass",
            "none: a budget of corners adding to the vev", "C", choice_bits=SR_BITS, meta={"rule": list(map(str, rule)), "rule_t": rule})
    g.preds["mH/v"] = (sumrule_mh(rule) / v, TREE)
    g = add("P9", "yt", "Top at the dial's unit", "y_t = 1, so mt = v/√2", "k sets the top Yukawa",
            "unit coupling = the lift's own scale", "C", choice_bits=G_BITS, meta={"yt2": 1.0})
    g.preds["mt/v"] = (1 / math.sqrt(2), TREE)
    g = add("P10", "angle", "Weinberg angle as a dial quarter", "sin²θ_W = 1/4", "k mixes hypercharge and weak isospin",
            "a quarter turn of the dial's rotation", "C", choice_bits=G_BITS, meta={"s2": 0.25})
    g.preds["sin2W"] = (0.25, 0.005)
    s2u, MU = unification_s2w()
    g = add("P11", "unif", "One dial for all three forces", "the three gauge couplings meet (SM running)",
            "k is the single coupling above the meeting scale", "one dial, three rungs", "B", meta={"M_meet": MU})
    g.preds["sin2W"] = (s2u, 0.003)
    a_num = 1 / (4 * math.pi ** 3 + math.pi ** 2 + math.pi)
    g = add("P12", "alpha", "Fine-structure constant from π", "1/α = 4π³ + π² + π", "k = electric charge",
            "none (a known numerology)", "C", choice_bits=12.0)
    g.preds["alpha"] = (a_num, 0.0)

    # ---- gravity and cosmology
    g = add("C1", "gravity", "Dial alone as gravity", "k(x) = 1/c(x)², clocks only", "k couples to time",
            "a varying dial changes clock rates (relativity pass)", "A")
    g.preds["gamma"] = ("lin", 0.0, GAMMA, sGAMMA, 3.0)
    g = add("C2", "gravity", "Dial boosts plus the river (= general relativity)", "local boosts by the dial, river at escape speed",
            "k couples to clocks and rulers equally", "time and space deficits are circle inverses (proved)", "A")
    g.preds["gamma"] = ("lin", 1.0, GAMMA, sGAMMA, 3.0)
    g = add("C3", "lightdial", "Light's dial rolling as dark energy", "exponential V, λ = √(2/3), α ∝ √k",
            "k = 1/c² couples to light", "the same k sets c and the potential", "B", meta={"lam": math.sqrt(2 / 3)})
    g = add("C4", "bd", "Dial coupled to curvature (Brans–Dicke)", "k·R, coupling ω free", "k couples to curvature",
            "dial as a local G", "B", n_free=1)
    g.fitted = ["gamma"]
    g = add("C5", "quint", "Scale-symmetric dial today, slope free", "V = V0 e^(−λφ), φ = √(3/2) ln k",
            "k couples only to itself and gravity", "rescaling symmetry (proved) makes ln k the natural field, so V is exponential",
            "A", n_free=1, meta={"lam": None})
    g = add("C6", "quint", "Dial today, V ∝ k", "V ∝ k, λ = √(2/3) = 0.816", "k is the metric's conformal factor",
            "conformal factor's kinetic term fixes the slope at n·√(2/3), n = 1", "B", meta={"lam": math.sqrt(2 / 3)})
    g = add("C7", "quint", "Dial today, V ∝ k²", "V ∝ k², λ = 2√(2/3) = 1.633", "k is the metric's conformal factor",
            "n = 2: an energy density scales as k² when k ∝ 1/length²", "B", meta={"lam": 2 * math.sqrt(2 / 3)})
    g = add("C8", "hier", "Dark energy from the weak scale", "ρΛ^(1/4) = v²/M_Pl", "k's dip height set by the vev",
            "seesaw between the two scales", "C")
    g.preds["rhoL/MPl"] = (v * v / MPL_RED / MPL_RED, 0.05)
    return G


def finish_cosmo(G):
    for g in G:
        if g.family == "quint" or g.family == "quint_shape":
            shape = g.meta.get("shape", "exp")
            cosmo_score(g, shape, g.meta.get("par", g.meta.get("lam")))
        if g.family == "lightdial":
            f = cosmo_fit("exp", g.meta["lam"], "both")
            da, dphi = dalpha_if_light_dial("exp", g.meta["lam"], f["Om"])
            g.extra.update({"dalpha_pred": da, "dphi": dphi})
            g.preds["dalpha"] = ("lin", da, DALPHA, sDALPHA, 2.0)
            # how weakly would the dial have to couple to light to pass?  |n_alpha| bound
            g.extra["n_alpha_bound_2sigma"] = (abs(DALPHA) + 2 * sDALPHA) / abs(dphi / math.sqrt(1.5))


# ------------------------------------------------------------------ children (refinements) of a guess
def children(g, used, rnd=2, seen_vals=None):
    seen_vals = seen_vals if seen_vals is not None else set()
    out = []
    fam = g.family
    def new(**kw):
        c = Guess(**kw); out.append(c); return c
    if fam == "tree_lam" or fam == "yt" or fam == "angle":
        key = {"tree_lam": "lam", "yt": "yt2", "angle": "s2"}[fam]
        x0 = g.meta.get(key)
        if x0 is None:
            return out
        # four grammar neighbours of the parent's constant (neighbours in value, not chosen by the data)
        vals = [v_ for v_, _ in GRAM]
        i = min(range(len(vals)), key=lambda j: abs(math.log(vals[j] / x0)))
        cand = [j for j in (i - 2, i - 1, i + 1, i + 2) if 0 <= j < len(vals)]
        for j in cand:
            val, lab = GRAM[j]
            gid = f"{fam}:{lab}"
            if gid in used or (fam, round(val, 9)) in seen_vals:
                continue
            seen_vals.add((fam, round(val, 9)))
            base = {"tree_lam": "Higgs from a grammar quartic", "yt": "Top from a grammar Yukawa", "angle": "Angle from a grammar constant"}[fam]
            c = new(gid=gid, family=fam, name=f"{base}: {key} = {lab}", shape=g.shape.split("(")[0] + f"({key} = {lab})",
                    coupling=g.coupling, geometry="grammar neighbour of the parent constant", grade="C",
                    choice_bits=G_BITS, meta={key: val})
            if fam == "tree_lam":
                c.preds["mH/v"] = (math.sqrt(2 * val), TREE)
            elif fam == "yt":
                c.preds["mt/v"] = (math.sqrt(val / 2), TREE)
            else:
                c.preds["sin2W"] = (val, 0.005)
    elif fam == "sumrule":
        rule = tuple(g.meta["rule_t"]) if "rule_t" in g.meta else (1, 1, 1, 1, Fraction(1))
        if rnd > 2:
            return out
        nbrs = []
        for r in SUM_RULES:
            d = sum(a != b for a, b in zip(r, rule))
            if d == 1:
                nbrs.append(r)
        nbrs = sorted(nbrs, key=lambda r: SUM_RULES.index(r))[:4]     # first four in grammar order: no data peeking
        for r in nbrs:
            gid = f"SR:{sumrule_label(r)}"
            if gid in used:
                continue
            c = new(gid=gid, family=fam, name=f"Sum rule {sumrule_label(r)}", shape=sumrule_label(r), coupling=g.coupling,
                    geometry="neighbour in the sum-rule grammar", grade="C", choice_bits=SR_BITS,
                    meta={"rule": [str(x) for x in r], "rule_t": r})
            c.preds["mH/v"] = (sumrule_mh(r) / v, TREE)
    elif fam == "flat":
        for scale, beta in (("MPl_red", False), ("GUT", False), ("string", False), ("MPl_red", True),
                            ("GUT", True), ("seesaw", False)):
            gid = f"flat:{scale}:{'beta' if beta else 'lam'}"
            if gid in used or gid in seen_vals or (scale == g.meta["scale"] and beta == g.meta["beta"]):
                continue
            seen_vals.add(gid)
            if len(out) >= 4:
                break
            fd = flat_dial_prediction(scale, with_beta=beta)
            c = new(gid=gid, family=fam, name=f"Dial flat{' and level' if beta else ''} at the {scale} scale",
                    shape=f"λ({scale}) = 0" + (" and β_λ = 0" if beta else ""), coupling=g.coupling,
                    geometry=g.geometry, grade="B", choice_bits=math.log2(6), meta={"scale": scale, "beta": beta})
            c.preds["mH/v"] = (fd["mh"] / v, fd["mh_err"] / fd["mh"])
            if beta:
                c.preds["mt/v"] = (fd["mt"] / v, fd["mt_err"] / fd["mt"])
    elif fam in ("quint", "quint_shape"):
        lam = g.meta.get("lam")
        if lam is None:
            lam = g.extra.get("full", {}).get("par", 0.8)
        step = 0.3 / 2 ** (rnd - 2)
        if fam == "quint":
            opts = [("exp", round(max(0.0, lam - step), 4)), ("exp", round(lam + step, 4))]
            if rnd == 2:
                opts += [("hilltop", None), ("cosine", None)]
        else:
            opts = []
        for shape, par in opts:
            gid = f"quint:{shape}:{'free' if par is None else round(par, 3)}"
            pk = ("quint", shape, None if par is None else round(par, 2))
            if gid in used or pk in seen_vals:
                continue
            seen_vals.add(pk)
            nm = {"exp": f"Dial today, exponential slope {par if par is None else round(par, 3)}",
                  "hilltop": "Dial today near the top of its hump (hilltop), width free",
                  "cosine": "Dial today on a circle (q < 0, compact dial), period free"}[shape]
            c = new(gid=gid, family="quint" if shape == "exp" else "quint_shape", name=nm,
                    shape={"exp": "V0 e^(−λφ)", "hilltop": "V0 (1 − φ²/2μ²)", "cosine": "V0 (1 + cos φ/f)"}[shape],
                    coupling="k couples only to itself and gravity",
                    geometry={"exp": "exponential from rescaling symmetry",
                              "hilltop": "the 1D sombrero: field sits near the hump at k = 0 and starts to roll",
                              "cosine": "negative dial: the dial's rotations close on a circle"}[shape],
                    grade="B", n_free=1 if par is None else 0,
                    meta={"shape": shape, "par": par, "lam": par})
    return out


# ------------------------------------------------------------------ chance calibration for formula families
def chance_formula(family_vals, target, th, n_trials=4000, seed=1):
    """Fraction of random targets (log-uniform over the prior) for which the best of `family_vals`
    scores at least `real_best` bits.  family_vals: list of predicted values for the target."""
    x, se, W, _ = TARGETS[target]
    s = math.sqrt(se ** 2 + th ** 2)
    rng = np.random.default_rng(seed)
    lnv = np.log(np.array(family_vals))
    lo = math.log(x) - W / 2
    best = []
    for _ in range(n_trials):
        t = lo + W * rng.random()
        z = (lnv - t) / s
        best.append(math.log2(W / (s * math.sqrt(2 * math.pi))) - np.min(z * z) / (2 * LN2))
    return np.array(best)


# ------------------------------------------------------------------ main
def main():
    t0 = time.time()
    say("=" * 100)
    say("Dial potential search, Nov 25 2026.  Scores in bits; ruled out = some prediction more than 3σ off.")
    say(f"grammar of constants p/q·π^a: {len(GRAM)} members ({G_BITS:.1f} bits to pick one); "
        f"sum-rule grammar: {len(SUM_RULES)} rules ({SR_BITS:.1f} bits)")
    say("=" * 100)
    say("Targets:")
    for k, (x, s, W, note) in TARGETS.items():
        say(f"  {k:9s} = {x:.6g}  (±{s:.1e} rel)   {note}")
    say(f"  gamma     = 1 + (2.1 ± 2.3)e-5 (Cassini);  Δα/α(z=1.15) = (1.3 ± 1.4)e-6 (ESPRESSO)")
    say(f"  Pantheon+ supernovae ({len(Q.mb)}) and DESI DR2 BAO ({len(Q.bv)} numbers)")

    rounds, used, n_tried = [], set(), 0
    global SEEN
    SEEN = set()
    G = round1()
    for g in G:
        if g.family == "quint" and g.meta.get("lam") is not None:
            SEEN.add(("quint", "exp", round(g.meta["lam"], 2)))
        if g.family == "flat":
            SEEN.add(f"flat:{g.meta['scale']}:{'beta' if g.meta['beta'] else 'lam'}")
    finish_cosmo(G)
    for rnd in range(1, 5):
        for g in G:
            used.add(g.gid)
            g.score()
        n_tried += len(G)
        for g in G:
            g.net = g.bits - g.choice_bits - math.log2(n_tried)
        G.sort(key=lambda g: (-g.bits if math.isfinite(g.bits) else 1e9))
        say("")
        say(f"---- Round {rnd}: {len(G)} guesses (search so far: {n_tried} guesses, look-elsewhere {math.log2(n_tried):.1f} bits)")
        say(f"{'id':34s} {'bits':>8s} {'net':>8s} {'max|z|':>7s}  verdict   guess")
        for g in G:
            verdict = "ruled out" if g.ruled_out else ("no test" if not g.rows and "cosmo_bits" not in g.extra else "alive")
            b = f"{g.bits:8.1f}" if g.bits > -1e5 else "   -huge"
            nb = f"{g.net:8.1f}" if g.net > -1e5 else "   -huge"
            say(f"{g.gid[:34]:34s} {b} {nb} {g.max_z:7.1f}  {verdict:9s} {g.name}")
            for r in g.rows:
                if r['pred'] is None or not math.isfinite(r['bits']):
                    say(f"{'':36s}   {r['target']:9s} no real prediction (mass squared negative)")
                else:
                    bb = f"{r['bits']:+.1f}" if r['bits'] > -1e5 else "< -100000"
                    say(f"{'':36s}   {r['target']:9s} pred {r['pred']:.6g}  z = {r['z']:+.2f}  ({bb} bits)")
            if "full" in g.extra:
                f = g.extra["full"]
                say(f"{'':36s}   full fit: slope/width {f['par']}, Ωm {f['Om']:.3f}, Δχ² vs Λ {f['dchi2']:+.2f}, "
                    f"w0 {f['w0']:.3f}; held-out Δχ² (SN→BAO + BAO→SN) {g.extra['heldout_dchi2_sum']:+.2f}")
            if "dalpha_pred" in g.extra:
                say(f"{'':36s}   Δα/α predicted {g.extra['dalpha_pred']:.3g}; coupling to light must be below "
                    f"{g.extra['n_alpha_bound_2sigma']:.2g} of natural to pass")
        rounds.append([g.as_dict() for g in G])
        if rnd == 4:
            break
        # keep the closest: the five best-scoring guesses with something to refine
        keep = [g for g in G if math.isfinite(g.bits) and (g.rows or "cosmo_bits" in g.extra)][:5]
        fams = {k.family for k in keep}
        for g in G:                       # near misses (|z| < 4) keep their family alive
            if g.family not in fams and g.max_z < 4 and (g.rows or "cosmo_bits" in g.extra):
                keep.append(g); fams.add(g.family)
        say(f"  kept for refinement: {', '.join(k.gid for k in keep)}")
        nxt = []
        for k in keep:
            nxt += children(k, used, rnd + 1, SEEN)
        extra = [g for g in G if g not in keep and math.isfinite(g.bits)]
        for g in extra:
            if len(nxt) >= 20:
                break
            nxt += children(g, used, rnd + 1, SEEN)
        nxt = nxt[:20]
        finish_cosmo(nxt)
        G = nxt
        if not G:
            break

    # ---------------- held-out and chance checks on the survivors
    say("")
    say("=" * 100)
    say("Held-out and chance checks")
    allg = [g for r in rounds for g in r]
    res = {"rounds": rounds}
    # 1. formula families: chance calibration
    lam_vals = [math.sqrt(2 * v_) for v_, _ in GRAM]
    best_tree = max((g for g in allg if g["family"] == "tree_lam" and g["rows"]), key=lambda g: g["bits"])
    dist = chance_formula(lam_vals, "mH/v", TREE)
    p = float(np.mean(dist >= best_tree["bits"]))
    say(f"  m_H from a grammar constant: best real hit {best_tree['name']} ({best_tree['bits']:.1f} bits). "
        f"Random targets do as well {100 * p:.0f}% of the time.")
    res["chance_tree_lam"] = {"best": best_tree["id"], "bits": best_tree["bits"], "p": p}
    sr_vals = [sumrule_mh(r) / v for r in SUM_RULES if np.isfinite(sumrule_mh(r))]
    best_sr = max((g for g in allg if g["family"] == "sumrule"), key=lambda g: g["bits"])
    dist = chance_formula(sr_vals, "mH/v", TREE)
    p = float(np.mean(dist >= best_sr["bits"]))
    say(f"  m_H from a mass sum rule: best real hit {best_sr['name']} ({best_sr['bits']:.1f} bits). "
        f"Random targets do as well {100 * p:.0f}% of the time.")
    res["chance_sumrule"] = {"best": best_sr["id"], "bits": best_sr["bits"], "p": p}
    yt_vals = [math.sqrt(v_ / 2) for v_, _ in GRAM]
    best_yt = max((g for g in allg if g["family"] == "yt"), key=lambda g: g["bits"])
    dist = chance_formula(yt_vals, "mt/v", TREE)
    p = float(np.mean(dist >= best_yt["bits"]))
    say(f"  m_t from a grammar Yukawa: best real hit {best_yt['name']} ({best_yt['bits']:.1f} bits). "
        f"Random targets do as well {100 * p:.0f}% of the time.")
    res["chance_yt"] = {"best": best_yt["id"], "bits": best_yt["bits"], "p": p}

    # 2. flat dial: m_t as the held-out number
    fdb = flat_dial_prediction("MPl", with_beta=True)
    fd = flat_dial_prediction("MPl")
    say(f"  Flat dial at M_Pl (λ = 0): m_H = {fd['mh']:.2f} ± {fd['mh_err']:.2f} GeV (NNLO fit; our two-loop run "
        f"{fd['mh_own_rge']:.2f}) vs measured {MH} ± {sMH}: {(fd['mh'] - MH) / math.hypot(fd['mh_err'], sMH):+.1f}σ")
    say(f"  Flat and level at M_Pl (λ = β_λ = 0): m_H = {fdb['mh']:.2f} ± {fdb['mh_err']:.2f}, "
        f"m_t = {fdb['mt']:.2f} ± {fdb['mt_err']:.2f} GeV (held out: measured m_t {MT} ± {sMT}, "
        f"{(fdb['mt'] - MT) / math.hypot(fdb['mt_err'], sMT):+.1f}σ)")
    # where does lambda cross zero, and where is beta zero, with measured inputs?
    mus, ys = R.run(MH, MT, AS, R.MPL_FULL, n=2000)
    lam = ys[4]
    i0 = np.argmax(lam < 0) if np.any(lam < 0) else None
    bet = np.array([R.betas(0, ys[:, j])[4] for j in range(ys.shape[1])])
    jb = int(np.argmin(np.abs(bet[200:]))) + 200
    say(f"  With measured inputs: λ turns negative at μ ≈ {mus[i0]:.2e} GeV; β_λ is smallest at μ ≈ {mus[jb]:.2e} GeV "
        f"(λ there {lam[jb]:+.4f}); λ(M_Pl) = {lam[-1]:+.4f}")
    res["flat_dial"] = {"lambda0": fd, "lambda_beta0": fdb, "lambda_negative_from": float(mus[i0]),
                        "beta_min_at": float(mus[jb]), "lambda_at_beta_min": float(lam[jb]),
                        "lambda_MPl_measured_inputs": float(lam[-1])}
    np.savez(os.path.join(HERE, "running.npz"), mu=mus, y=ys)

    # 3. how well is the exponential slope pinned down?  profile chi2 over lambda (Om refitted)
    prof = []
    for lam_ in np.linspace(0.0, 1.6, 17):
        f = cosmo_fit("exp", float(lam_), "both")
        prof.append((float(lam_), f["chi2"]))
    lams = np.array([p_[0] for p_ in prof]); ch = np.array([p_[1] for p_ in prof])
    fine = np.linspace(0, 1.6, 1601); chf = np.interp(fine, lams, ch)
    i_ = int(np.argmin(chf)); inside = fine[chf <= chf[i_] + 1]
    say(f"  Exponential slope from all data: best λ = {fine[i_]:.2f}, 1σ range {inside.min():.2f} to {inside.max():.2f}; "
        f"√(2/3) = {math.sqrt(2/3):.3f}; Λ (λ = 0) is Δχ² = {ch[0] - chf[i_]:.1f} worse")
    res["slope_profile"] = {"lam": lams.tolist(), "chi2": ch.tolist(), "best": float(fine[i_]),
                            "one_sigma": [float(inside.min()), float(inside.max())]}
    # chance that a 'geometric' slope n·sqrt(2/3) or sqrt(n/m)-type value lands inside the 1σ range
    geo = sorted({round(math.sqrt(p_ / q_), 6) for p_ in range(1, 9) for q_ in range(1, 9)})
    hits = [g_ for g_ in geo if inside.min() <= g_ <= inside.max()]
    say(f"  Of {len(geo)} simple slopes √(p/q), p, q ≤ 8, {len(hits)} fall inside that 1σ range")
    res["slope_geo_hits"] = {"n": len(geo), "inside": hits}

    # 4. held-out test of the conformal reading: the factor that gives slope n·sqrt(2/3) also couples to matter
    a0sq = 1 / 6
    gam = 1 - 2 * a0sq / (1 + a0sq)
    zg = (gam - GAMMA) / sGAMMA
    say(f"  Conformal dial (slope n·√(2/3)) couples to matter with α0² = 1/6, so PPN γ = 1 − 2α0²/(1+α0²) = 5/7 = {gam:.4f}: "
        f"Cassini puts that {abs(zg):.0f}σ away unless the force is screened")
    say(f"  So the conformal slope survives only if the dial is screened in the solar system (as f(R) gravity must be), "
        f"or if the dial is dark-only, in which case the slope is free (best {fine[i_]:.2f})")
    res["conformal_gamma"] = {"gamma": gam, "z": zg}

    # 5. the 'light dial' bound in words
    res["targets"] = {k: {"value": x, "rel_sigma": s, "prior_efolds": W, "note": n_} for k, (x, s, W, n_) in TARGETS.items()}
    res["n_guesses"] = n_tried
    res["runtime_s"] = time.time() - t0
    json.dump(res, open(os.path.join(HERE, "results.json"), "w"), indent=1, default=str)
    open(os.path.join(HERE, "search_output.txt"), "w").write("\n".join(OUT) + "\n")
    say(f"done in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
