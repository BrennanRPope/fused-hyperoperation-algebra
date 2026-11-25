"""Pictures for the height-flow write-up. Run after height_entropy.py (reads height_entropy.json and
../height-running-coupling/running_coupling.json). Writes figures/f1 ... f8."""
import os, sys, json, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "fractional-rungs"))
import height_entropy as he   # re-runs the checks (about 10 s) and gives the same functions
from regular_iteration import Es, B
FIG = os.path.join(HERE, "figures"); os.makedirs(FIG, exist_ok=True)
D = json.load(open(os.path.join(HERE, "height_entropy.json")))
RC = json.load(open(os.path.join(HERE, "..", "height-running-coupling", "running_coupling.json")))

# reference palette (dataviz skill, light mode): categorical in fixed order; heights on one sequential hue
BLUE, ORANGE, AQUA, YELLOW, MAGENTA, GREEN, VIOLET, RED = "#257de0", "#f47230", "#00b894", "#f2a600", "#ee7fac", "#008a17", "#6a4ce2", "#ec4647"
INK, INK2, GRID = "#191919", "#666666", "#e6e6e6"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
                     "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 2,
                     "legend.frameon": False, "axes.titleweight": "bold", "axes.titlesize": 11})
def hcol(s, lo=-0.5, hi=3.0):
    """one hue, light -> dark with height"""
    return cm.Blues(0.35 + 0.6*(s - lo)/(hi - lo))
def save(fig, name):
    fig.savefig(os.path.join(FIG, name), dpi=150, bbox_inches="tight"); plt.close(fig); print("wrote", name)
G = he.G

# ---------- f1: the logs and where they sit among known group logarithms ----------
fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
X = np.linspace(1.0, 8, 400)
for s in [-0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 3.0]:
    y = he.ln_s(s, X)
    ax[0].plot(X, y, color=hcol(s), label=f"s = {s:g}" + {0.0: "  (X − 1)", 1.0: "  (ln X)", 2.0: "  (ln(1+ln X))"}.get(s, ""))
ax[0].set_ylim(-0.2, 7.5); ax[0].set_xlabel("X  (1/p, or a number of states)"); ax[0].set_ylabel("ln_s(X)")
ax[0].set_title("A. Supernatural logs on X ≥ 1 (all start at (1, 0), slope 1)")
ax[0].legend(fontsize=8.5, loc="upper left")
t = np.linspace(0, 3, 300); r = 0.5
ax[1].plot(t, t, color=INK2, lw=1.2, ls=":", label="G(t) = t  (Shannon / Rényi, known)")
ax[1].plot(t, G(1.0, t), color=INK2, lw=1.2, ls="--", label="G(t) = e^t − 1  (Tsallis q = 0, known)")
ax[1].plot(t, np.expm1(r*t)/r, color=ORANGE, label="Tsallis / Sharma–Mittal dial, k = 1/2  (known)")
ax[1].plot(t, np.sinh(r*t)/r, color=AQUA, label="Kaniadakis, κ = 1/2  (known)")
ax[1].plot(t, G(r, t), color=BLUE, lw=2.6, label="height: G = E^(1/2), half-rung  (ours)")
ax[1].set_ylim(0, 8); ax[1].set_xlabel("t = ln X"); ax[1].set_ylabel("G(t)")
ax[1].set_title("B. Each entropy family is one choice of G in ln_G(X) = G(ln X)")
ax[1].legend(fontsize=8.5, loc="upper left")
save(fig, "f1_logs_and_G.png")

# ---------- f2: composition for independent systems ----------
rng = np.random.default_rng(3)
s = 0.5; rr = 1 - s
pts = {"Hill form H (ours)": [], "trace form T": [], "escort form (q = 2)": []}
for _ in range(60):
    pa, pb = he.rand_p(rng.integers(2, 7)), he.rand_p(rng.integers(2, 7)); pab = np.outer(pa, pb).ravel()
    pts["Hill form H (ours)"].append((float(B(rr, he.H(s, pa), he.H(s, pb))), he.H(s, pab)))
    pts["trace form T"].append((float(B(rr, he.T(s, pa), he.T(s, pb))), he.T(s, pab)))
    pts["escort form (q = 2)"].append((float(B(rr, he.Esc(s, pa, 2), he.Esc(s, pb, 2))), he.Esc(s, pab, 2)))
fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
for (name, P), c, m in zip(pts.items(), [BLUE, ORANGE, AQUA], ["o", "s", "^"]):
    P = np.array(P); ax[0].scatter(P[:, 0], 100*(P[:, 1] - P[:, 0])/P[:, 0], s=26, color=c, marker=m, edgecolor="white", linewidth=0.8, label=name, zorder=3)
ax[0].axhline(0, color=INK2, lw=1, zorder=1)
ax[0].set_xlabel("predicted by our operation: B_(1−s)(S_A, S_B)"); ax[0].set_ylabel("actual minus predicted, % of predicted")
ax[0].set_title("A. Two independent systems at height s = 1/2"); ax[0].legend(fontsize=9, loc="lower left")
ax[0].text(0.98, 0.93, "on the zero line = combines exactly", transform=ax[0].transAxes, ha="right", color=INK2, fontsize=9)
comp = D["composition_errors"]; hs = sorted(comp["H"], key=float)
xs = np.arange(len(hs)); w = 0.2
for i, (f, c, lab) in enumerate([("H", BLUE, "Hill form H"), ("Z", VIOLET, "Tempesta Z form"), ("T", ORANGE, "trace form T"), ("Esc", AQUA, "escort form")]):
    vals = [max(comp[f][h], 1e-13) for h in hs]
    ax[1].bar(xs + (i - 1.5)*w, vals, width=w*0.9, color=c, label=lab)
ax[1].set_yscale("log"); ax[1].set_xticks(xs); ax[1].set_xticklabels([f"{float(h):g}" for h in hs])
ax[1].axhline(1e-7, color=INK2, lw=1, ls="--"); ax[1].text(len(hs) - 1.2, 2e-7, "exact (rounding only) below", fontsize=8, color=INK2, ha="right")
ax[1].set_xlabel("height s"); ax[1].set_ylabel("worst relative error of the combining law")
ax[1].set_title("B. Exact at every height for H and Z; trace and escort only at s = 0, 1")
ax[1].legend(fontsize=8.5, loc="center left")
save(fig, "f2_composition.png")

# ---------- f3: a coin, at each height; and the height flow of the entropy ----------
fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
pp = np.linspace(1e-4, 1 - 1e-4, 400)
for s in [0.25, 0.5, 1.0, 1.5, 2.0, 3.0]:
    y = np.array([he.T(s, np.array([q, 1 - q])) for q in pp])
    ax[0].plot(pp, y/y.max(), color=hcol(s), label=f"s = {s:g}" + ("  (Shannon)" if s == 1 else ""))
ax[0].set_xlabel("chance of heads p"); ax[0].set_ylabel("T_s(coin) ÷ its maximum")
ax[0].set_title("A. Trace form on a coin: same peak, flatter shoulders at low height")
ax[0].legend(fontsize=8.5, loc="lower center")
ss = np.linspace(-0.5, 12, 300)
for R, c in zip([0.4, 1.3, 2.5], [AQUA, BLUE, VIOLET]):
    Hs = np.array([float(G(1 - s, np.array([R]))[0]) for s in ss])
    ax[1].plot(ss, Hs, color=c, label=f"Shannon value R = {R:g} at s = 1")
    cst = he.Phi(R) + 1
    ss2 = ss[ss > cst + 0.6]
    ax[1].plot(ss2, 2/(ss2 - cst), color=c, lw=1.2, ls="--")
ax[1].scatter([1, 1, 1], [0.4, 1.3, 2.5], color=INK, s=18, zorder=4)
ax[1].set_ylim(0, 4); ax[1].set_xlabel("height s"); ax[1].set_ylabel("H_s = E^(1−s)(R)")
ax[1].set_title("B. Moving along height runs the entropy like a coupling")
ax[1].text(5.5, 1.1, "dashed: one-loop law 2/(s − c),\nc = Φ(R) + 1 the invariant", fontsize=9, color=INK2)
ax[1].legend(fontsize=8.5, loc="upper right")
save(fig, "f3_coin_and_flow.png")

# ---------- f4: known vs ours along height, and the state-space ladder ----------
fig, ax = plt.subplots(1, 2, figsize=(12, 4.6), gridspec_kw={"width_ratios": [1.1, 1]})
a = ax[0]; a.grid(False); a.set_yticks([]); a.spines["left"].set_visible(False)
a.plot([-1, 3.2], [0, 0], color=INK2, lw=1)
for lo, hi in [(-1, 0), (0, 1), (1, 2), (2, 3.2)]:
    a.plot([lo + 0.05, hi - 0.05], [0, 0], color=BLUE, lw=7, solid_capstyle="butt", alpha=0.55)
known = [(0, 0.55, "s = 0\nTsallis q = 0,\nSharma–Mittal β = 0"),
         (1, 0.12, "s = 1\nShannon, Rényi"),
         (2, 0.55, "s = 2\nlog(1 + Rényi),\nTempesta–Jensen 2019")]
for x, yy, lab in known:
    a.scatter([x], [0], s=90, color=ORANGE, zorder=4, edgecolor="white", linewidth=2)
    a.plot([x, x], [0.04, yy - 0.02], color=INK2, lw=0.8)
    a.text(x, yy, lab, ha="center", va="bottom", fontsize=8.5, color=INK)
a.text(1.1, -0.12, "blue: every fractional height is a new member (ours),\nand one flow links them all (ours)", ha="center", va="top", fontsize=9, color=BLUE)
a.set_xlim(-1.1, 3.3); a.set_ylim(-0.5, 1.15); a.set_xlabel("height s")
a.set_title("A. Known points (orange) and what is ours (blue)")
N = np.geomspace(1, 1000, 300)
for s, lab in [(0.0, "s = 0: W ≈ N (one walker on N sites)"), (0.5, "s = 1/2: between polynomial and exponential"),
               (1.0, "s = 1: W = e^N (N independent units)"), (1.5, "s = 3/2: between exponential and double exp."),
               (2.0, "s = 2: W ~ exp(e^N), like Boolean functions of N bits")]:
    lnW = G(s - 1, N)
    ax[1].plot(N, np.log10(lnW), color=hcol(s), label=lab)
ax[1].set_xscale("log"); ax[1].set_ylim(-0.5, 12)
ax[1].set_xlabel("system size N"); ax[1].set_ylabel("log10( ln W(N) )")
ax[1].set_title("B. Where S_s grows in step with size: one height per tower rung")
ax[1].legend(fontsize=8, loc="upper left")
save(fig, "f4_known_vs_ours_and_ladder.png")

# ---------- f5: Lesche stability ----------
fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
deltas = [1e-2, 1e-3, 1e-4, 1e-5]
for s in [0.0, 0.25, 0.5, 1.0, 2.0]:
    sup = [max(D["lesche"][f"T|{s}"][str(d)]) for d in deltas]
    ax[0].plot(deltas, sup, marker="o", ms=6, color=hcol(s), label=f"s = {s:g}" + ("  (Tsallis q = 0)" if s == 0 else ""))
ax[0].set_xscale("log"); ax[0].set_yscale("log"); ax[0].invert_xaxis()
ax[0].set_xlabel("size of the nudge δ (shrinking →)"); ax[0].set_ylabel("worst |ΔS| ÷ S_max over all W tried")
ax[0].set_title("A. Trace form: error shrinks with the nudge\n(stable for every s > 0)")
ax[0].legend(fontsize=8.5)
LH = np.array(D["lesche_hill_L"])
for s in [0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0]:
    ax[1].plot(LH, D["lesche_hill_long"][str(s)], color=hcol(s), label=f"s = {s:g}", lw=2.6 if s == 1 else 2)
ax[1].set_xscale("log"); ax[1].set_ylim(0, 1.05)
ax[1].set_xlabel("ln W  (number of states, as a log)"); ax[1].set_ylabel("|ΔS| ÷ S_max at δ = 0.01")
ax[1].set_title("B. Hill form: only s = 1 stays flat\n(every other height climbs: unstable)")
ax[1].legend(fontsize=8, ncol=2, loc="upper left")
save(fig, "f5_lesche.png")

# ---------- f6: maximum-entropy distributions ----------
fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
for s in [0.25, 0.5, 0.75, 1.0]:
    e = np.array(D["maxent"][str(s)]["eps"]); p = np.array(D["maxent"][str(s)]["p"])
    m = (e > 0.5) & (p > 0)
    ax[0].plot(e[m], -np.log(p[m]), color=hcol(s), label=f"s = {s:g}" + ("  (Boltzmann)" if s == 1 else ""))
e = np.geomspace(0.5, 400, 200)
ax[0].plot(e, 2*np.log(e) + 2, color=ORANGE, lw=1.3, ls="--", label="power law (−ln p ∝ ln ε)")
ax[0].plot(e, 2*np.sqrt(e), color=AQUA, lw=1.3, ls="--", label="stretched exponential (√ε)")
ax[0].set_xscale("log"); ax[0].set_yscale("log"); ax[0].set_xlabel("energy ε"); ax[0].set_ylabel("−ln p(ε)")
ax[0].set_title("A. Max-entropy tails of the trace form"); ax[0].legend(fontsize=8.5)
for s in [0.25, 0.5, 0.75, 1.0]:
    e = np.array(D["maxent"][str(s)]["eps"]); p = np.array(D["maxent"][str(s)]["p"])
    m = (p > 0) & (e < 60)
    ax[1].plot(e[m], p[m], color=hcol(s), label=f"s = {s:g}")
ax[1].set_yscale("log"); ax[1].set_xlabel("energy ε"); ax[1].set_ylabel("p(ε)")
ax[1].set_title("B. Same, plain: heavier than Boltzmann, lighter than power laws")
ax[1].legend(fontsize=8.5)
save(fig, "f6_maxent.png")

# ---------- f7: the beta function and two-loop running (from the running-coupling handoff) ----------
fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
x = np.linspace(-1.5, 3, 400)
vx = np.array([he.v_of(xx) if xx > 0 else np.nan for xx in x])
xs_neg = x[x < 0]
# for x < 0 use v(x) = v(E^n x) / prod E'(E^k x) via the map toward 0
def v_neg(xx, n=4000):
    z = xx; logp = 0.0
    for _ in range(n):
        logp += z; z = math.expm1(z)
    return he.v_small(z)*math.exp(-logp)
vx[x < 0] = [v_neg(xx) for xx in xs_neg]; vx[x == 0] = 0
ax[0].plot(x, vx, color=BLUE, label="β(x) = v(x), Julia's iterative log of e^x − 1")
ax[0].plot(x, x**2/2, color=ORANGE, lw=1.3, ls="--", label="one loop: x²/2")
ax[0].plot(x, x**2/2 - x**3/12, color=AQUA, lw=1.3, ls="--", label="two loops: x²/2 − x³/12")
ax[0].axhline(0, color=INK2, lw=1)
for xx in [-1.0, -0.5, 0.6, 1.5]:
    ax[0].annotate("", xy=(xx + 0.18, 0.05), xytext=(xx - 0.18, 0.05), arrowprops=dict(arrowstyle="->", color=INK, lw=1.4))
ax[0].set_ylim(-0.3, 4); ax[0].set_xlabel("coupling x"); ax[0].set_ylabel("dx/ds")
ax[0].set_title("A. Height flow as a beta function (double zero at 0)"); ax[0].legend(fontsize=8.5, loc="upper left")
run = RC["running"]
n = np.array([r_["n"] for r_ in run]); xn = np.array([r_["x"] for r_ in run])
e1 = np.abs(np.array([r_["one_loop"] for r_ in run]) - xn); e2 = np.abs(np.array([r_["two_loop_with_constant"] for r_ in run]) - xn)
ax[1].plot(n, e1, marker="o", ms=5, color=ORANGE, label="one-loop formula −2/n")
ax[1].plot(n, np.maximum(e2, 1e-20), marker="s", ms=5, color=BLUE, label="two-loop formula (with the ρ = 1/3 log)")
ax[1].set_xscale("log"); ax[1].set_yscale("log"); ax[1].set_xlabel("steps below addition n"); ax[1].set_ylabel("error in x_n = E^n(−1)")
ax[1].set_title("B. Below addition the coupling runs like 2/|s|; two loops are visible")
ax[1].legend(fontsize=8.5)
save(fig, "f7_beta_and_running.png")

# ---------- f8: renormalons rotated by 90 degrees ----------
fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
ac = np.array(RC["abs_c"]); k = ac[:, 0]; lc = ac[:, 1] + ac[:, 0]*math.log10(2)   # coefficients in Z = -2/x, the Fatou variable
ref = np.array([math.lgamma(kk - 1)/math.log(10) - (kk - 1)*math.log10(2*math.pi) for kk in k])
off = np.median((lc - ref)[k > 60])
ax[0].plot(k, lc, color=BLUE, label="log10 2^k |c_k|  (β coefficients in the natural variable Z = −2/x)")
ax[0].plot(k, ref + off, color=ORANGE, lw=1.3, ls="--", label="(k − 2)! / (2π)^(k−1), one constant shift")
ax[0].set_xlabel("order k"); ax[0].set_ylabel("log10 2^k |c_k|")
ax[0].set_title("A. The beta series diverges factorially, at the 2π rate"); ax[0].legend(fontsize=8.5, loc="upper left")
pol = np.array(RC["pade_poles_generator"]); pol = pol[np.abs(pol[:, 0]) < 3]
ax[1].scatter(pol[:, 0], pol[:, 1], s=40, color=BLUE, edgecolor="white", linewidth=1.5, zorder=3, label="height flow of e^x − 1: Padé poles of the Borel transform")
for m_ in [1, 2, 3]:
    ax[1].scatter([0, 0], [2*math.pi*m_, -2*math.pi*m_], s=110, facecolor="none", edgecolor=BLUE, lw=1.2, zorder=2)
    ax[1].scatter([2*m_*1.6], [0], s=80, marker="D", color=ORANGE, edgecolor="white", linewidth=1.5, zorder=3)
ax[1].scatter([], [], s=80, marker="D", color=ORANGE, label="QFT renormalons (cartoon): on the real axis")
ax[1].scatter([], [], s=110, facecolor="none", edgecolor=BLUE, label="expected lattice ±2πi, ±4πi, …")
ax[1].axhline(0, color=INK2, lw=1); ax[1].axvline(0, color=INK2, lw=1)
ax[1].set_xlim(-3, 13); ax[1].set_ylim(-21, 21); ax[1].set_xlabel("Re (Borel plane)"); ax[1].set_ylabel("Im")
ax[1].set_title("B. Same mechanism, singularities turned 90°"); ax[1].legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=1)
save(fig, "f8_renormalons_rotated.png")
