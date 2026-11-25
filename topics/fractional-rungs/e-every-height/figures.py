import sys, os, math, json
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, ".."))
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from regular_iteration import Es, B
D = json.load(open(os.path.join(HERE, "e_every_height.json")))
FIG = os.path.join(HERE, "figures"); os.makedirs(FIG, exist_ok=True)
BLUE, RED, INK, MUTE = "#1d2bd8", "#cf3519", "#15131a", "#8a8794"
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})
def e_s(s): return 1 + float(Es(s, np.array([1.0]))[0])

# ---- 1. e_s against height ----
ss = np.linspace(-4, 2.75, 220); es = np.array([e_s(s) for s in ss])
fig, ax = plt.subplots(figsize=(8, 5.2))
ax.plot(ss, es, color=BLUE, lw=2.4)
ax.set_yscale("log"); ax.set_ylim(1, 40)
ax.axhline(1, color=MUTE, ls=":", lw=1)
marks = [(-1, "1 + log 2\n= 1.693"), (0, "2\n(addition)"), (0.5, "2.271\n(half rung)"), (1, "e = 2.718\n(our operation)"), (2, "e^(e-1)\n= 5.575")]
for s, lab in marks:
    v = e_s(s); ax.plot([s], [v], "o", color=RED, ms=7, zorder=5)
    off = {-1: (-20, 12), 0: (-30, 14), 0.5: (6, -34), 1: (-80, 14), 2: (-78, 4)}[s]
    ax.annotate(lab, (s, v), textcoords="offset points", xytext=off, fontsize=10, color=INK)
ax.annotate("far below, e_s  ->  1 like 1 + 2/|s|", (-3.9, 1.05), fontsize=10, color=MUTE)
ax.annotate("above 2 it climbs like a tower:\ne_3 = 97.0", (2.0, 14), fontsize=10, color=MUTE)
ax.set_xlabel("rung height s   (0 = addition, 1 = a+b+ab, 2 = the next rung up)")
ax.set_ylabel("e_s = 1 + E^s(1)   (log scale)")
ax.set_title("e at every rung height: one smooth curve through 2, e and e^(e-1)")
fig.tight_layout(); fig.savefig(os.path.join(FIG, "f1_e_by_height.png"), dpi=150); plt.close(fig)

# ---- 2. convergence of the n-fold build ----
fig, axs = plt.subplots(1, 2, figsize=(11, 4.6))
ax = axs[0]
cols = {"-1.0": "#3a9a5b", "0.5": BLUE, "1.0": RED, "2.0": "#8a3fd1"}
for s, c in cols.items():
    b = D["build"][s]; ns = [int(n) for n in b["values"]]; v = [b["values"][str(n)] for n in ns]
    ax.plot(ns, v, "o-", color=c, label=f"s = {float(s):g}"); ax.axhline(b["target"], color=c, ls="--", lw=1)
ax.set_xscale("log"); ax.set_xlabel("n  (slices of the unit; the n-cube)")
ax.set_ylabel("1 + (1/n combined n times with B_s)")
ax.set_title("Cut the unit finer: the build settles\non 1 + E^s(1) (dashed)")
ax.legend(frameon=False, fontsize=9)
ax = axs[1]
for s, c in cols.items():
    b = D["build"][s]; ns = np.array([int(n) for n in b["values"]]); v = np.array([b["values"][str(n)] for n in ns])
    ax.loglog(ns, np.abs(v - b["target"]), "o", color=c, label=f"s = {float(s):g}")
    ax.loglog(ns, abs(b["pred_n_err"])/ns, "-", color=c, lw=1)
ax.set_xlabel("n"); ax.set_ylabel("|build - (1 + E^s(1))|")
ax.set_title("Error falls like 1/n, exactly as predicted\n(lines: (s/2n)·slope of E^s at 1)")
ax.legend(frameon=False, fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "f2_build_convergence.png"), dpi=150); plt.close(fig)

# ---- 3. the natural log at each height ----
fig, axs = plt.subplots(1, 2, figsize=(11.5, 5))
ax = axs[0]
S = [-1, -0.5, 0, 0.5, 1, 1.5, 2]
PAL = ["#2b3a9e", "#4f7bd6", "#5aa9c9", "#3a9a5b", "#cf3519", "#b5473a", "#7a1f6e"]
for i, s in enumerate(S):
    X = np.linspace(0.05, 7, 500)
    with np.errstate(all="ignore"): Y = Es(-s, X - 1)
    c = PAL[i]; ax.plot(X, Y, color=c, lw=2.6 if s == 1 else 1.6, label=f"ln_{s:g}" + ("  = ordinary ln" if s == 1 else "  = X - 1" if s == 0 else ""))
    es_ = e_s(s)
    if es_ < 7: ax.plot([es_], [1], "o", color=c, ms=6)
ax.axhline(1, color=MUTE, ls=":", lw=1); ax.axhline(0, color=MUTE, lw=.6); ax.axvline(1, color=MUTE, lw=.6)
ax.plot([1], [0], "ko", ms=5)
ax.annotate("every curve passes (1, 0)\nwith slope 1; its bend there is -s", (1.05, -1.25), fontsize=10)
ax.annotate("dots: ln_s(e_s) = 1", (4.2, 1.08), fontsize=10, color=INK)
ax.set_xlim(0, 7); ax.set_ylim(-2, 3)
ax.set_xlabel("X"); ax.set_ylabel("ln_s(X) = E^-s(X - 1)")
ax.set_title("The supernatural logarithm ln_s: one family,\ncontinuous in the height s")
ax.legend(frameon=False, fontsize=8.5, loc="lower right")
ax = axs[1]
X = np.linspace(0.3, 6, 400)
for s, c in [(0.5, BLUE), (2, RED)]:
    with np.errstate(all="ignore"):
        ax.plot(X, Es(-s, X - 1), color=c, lw=2, label=f"ours, ln_{s:g}")
        ax.plot(X, (X**(1 - s) - 1)/(1 - s), color=c, lw=1.4, ls="--", label=f"Tsallis q-log, q = {s:g}")
ax.axvline(1, color=MUTE, lw=.6); ax.axhline(0, color=MUTE, lw=.6)
ax.set_ylim(-3, 2.5); ax.set_xlabel("X")
ax.set_title("Known cousin: Tsallis q-logarithm.\nSame value, slope and bend at 1 when q = s; they part after")
ax.legend(frameon=False, fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "f3_ln_family.png"), dpi=150); plt.close(fig)

# ---- 4. the closed-ish forms: power series in s, Levy formula ----
fig, axs = plt.subplots(1, 2, figsize=(11, 4.6))
ax = axs[0]
T = D["taylor"]; sg = np.linspace(-2.5, 2.6, 300)
ax.plot(ss[(ss > -2.5) & (ss < 2.6)], es[(ss > -2.5) & (ss < 2.6)], color=INK, lw=3, alpha=.35, label="e_s (computed by the flow)")
for K, c in [(2, "#c9b04a"), (4, "#3a9a5b"), (8, BLUE), (16, RED)]:
    ax.plot(sg, sum(T[m]*sg**m for m in range(K)), color=c, lw=1.4, label=f"first {K} terms of the series in s")
ax.set_ylim(1, 8); ax.set_xlabel("height s"); ax.set_ylabel("value")
ax.set_title("e_s as a power series in the height\n(converges for every height: an entire function)")
ax.legend(frameon=False, fontsize=8.5, loc="upper left")
ax = axs[1]
l1 = D["levy_e_half"]; l2 = D["levy2_e_half"]; tgt = e_s(0.5)
n1 = [int(k) for k in l1]; n2 = [int(k) for k in l2]
ax.loglog(n1, [abs(l1[str(n)] - tgt) for n in n1], "o-", color=BLUE, label="Levy: walk in n steps, bend by y/(1 - s·y/2), walk out")
ax.loglog(n2, [abs(l2[str(n)] - tgt) for n in n2], "s-", color=RED, label="same, with the log term of the Abel function")
ax.set_xlabel("steps n"); ax.set_ylabel("|error| at e_1/2")
ax.set_title("A one-line limit formula for any height\n(shown at the half rung)")
ax.legend(frameon=False, fontsize=8.5)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "f4_closed_forms.png"), dpi=150); plt.close(fig)
print("figures written")

# ---- 5. two survey leads: the height flow as a running coupling; Poisson "no event" floor ----
fig, axs = plt.subplots(1, 2, figsize=(11.5, 4.7))
ax = axs[0]
sneg = -np.logspace(-0.3, 2.6, 60)
run = np.array([float(Es(s, np.array([1.0]))[0]) for s in sneg])
ax.loglog(-sneg, run, "o", color=BLUE, ms=4, label="E^s(1), the height flow run downward")
ax.loglog(-sneg, 2/(-sneg), "-", color=RED, lw=1.5, label="2/|s|: one-loop running of a marginal coupling")
ax.set_xlabel("|s|, steps down the rungs  (like RG time)"); ax.set_ylabel("value")
ax.set_title("Below addition, the flow runs like a marginal\ncoupling: speed g²/2 - g³/12 + g⁴/48 - …")
ax.legend(frameon=False, fontsize=9)
ax = axs[1]
lam = np.linspace(0, 8, 200)
for s, c in [(0.5, "#4f7bd6"), (1, RED), (1.5, "#b5473a"), (2, "#7a1f6e")]:
    ax.plot(lam, 1 + Es(s, -lam), color=c, lw=2.4 if s == 1 else 1.6, label=f"height {s:g}" + ("  (e^-λ, true Poisson)" if s == 1 else ""))
ax.axhline(0, color=INK, lw=.8); ax.axhline(1/math.e, color=MUTE, ls=":", lw=1)
ax.annotate("floor 1/e at height 2", (5.2, 0.4), fontsize=9, color=MUTE)
ax.annotate("below height 1 it goes negative:\nnot a probability", (3.2, -0.38), fontsize=9, color=MUTE)
ax.set_xlabel("λ (expected number of events)"); ax.set_ylabel("chance of no event = 1 + E^s(-λ)")
ax.set_title("Only height 1 lets the chance of 'nothing happens'\nfall all the way to 0 without going negative")
ax.legend(frameon=False, fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "f5_survey_leads.png"), dpi=150); plt.close(fig)
print("figure 5 written")
