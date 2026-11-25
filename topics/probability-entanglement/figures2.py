"""Pictures for the DRAFT stress-test section. Run after checks2.py (reads results2.json)."""
import json, math
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon
from lib import circ
R = json.load(open("results2.json"))
BLUE, ORANGE, GREEN, RED, PURPLE, GRAY = "#2a6fdb", "#e8833a", "#2e9e5b", "#d64545", "#8a5cc2", "#888888"
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
def save(fig, name): fig.tight_layout(); fig.savefig(f"figures/{name}.png", dpi=140); plt.close(fig)
def diag(a, lo, hi): a.plot([lo, hi], [lo, hi], color=GRAY, lw=1, ls="--", zorder=0)

# ---------------- fig9: heat and spin
fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
a = ax[0]
for beta, col in [("15.0", RED), ("40.0", ORANGE)]:
    rows = np.array(R["heat"][beta])
    a.scatter(rows[:, 4], rows[:, 3], color=col, s=40, label=f"thermal clock, β={float(beta):.0f}")
    a.scatter(rows[:, 5], rows[:, 3], color=col, marker="x", s=40, label=f"vacuum clock (wrong), β={float(beta):.0f}")
diag(a, 0, .3); a.set_xlabel("clock prediction t/3"); a.set_ylabel("measured CMI (free fermions)")
a.set_title("(a) Heated state: the clock survives,\nbut it is the thermal (sinh) clock"); a.legend(fontsize=8)
a = ax[1]; rot = np.array(R["rot"])
a.scatter(rot[:, 4], rot[:, 3], color=PURPLE, s=50, label="(t_L + t_R)/6: two half-rate clocks")
a.scatter(rot[:, 5], rot[:, 3], color=RED, marker="v", s=40, label="left clock alone, t_L/3")
a.scatter(rot[:, 6], rot[:, 3], color=BLUE, marker="^", s=40, label="right clock alone, t_R/3")
diag(a, 0, .3); a.set_xlabel("prediction"); a.set_ylabel("measured CMI")
a.set_title("(b) Spinning state (left movers 4× hotter):\nx = y_L ∘₋₁ y_R, each clock at half speed"); a.legend(fontsize=8)
a = ax[2]; run = np.array(R["run"])
a.plot(run[:, 0], run[:, 1], color=RED, lw=2)
a.axhline(1, color=GRAY, ls=":"); a.set_xlabel("region size / β (thermal length)"); a.set_ylabel("thermal clock ÷ vacuum clock")
a.set_title("(c) Heat changes the clock's SPEED with scale\n(it runs, like RG); the gluing dial k stays −1")
save(fig, "fig9_heat_and_spin")

# ---------------- fig10: where the geometric reading breaks; heat vs squeezing
fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
a = ax[0]; plat = np.array(R["btz"]["plat"]); ts = R["btz"]["th_star"]
a.plot(plat[:, 0] / math.pi, plat[:, 1], color=RED, lw=2, label="true CMI (with the plateau)")
a.plot(plat[:, 0] / math.pi, plat[:, 2], color=BLUE, ls="--", lw=2, label="clock value t")
a.axvline(ts / math.pi, color=GRAY, ls=":"); a.text(ts / math.pi - .02, a.get_ylim()[1] * .8, "plateau", ha="right", color=GRAY)
a.set_xlabel("A+B+C (units of π; whole circle = 2)"); a.set_ylabel("CMI (units c/3)"); a.set_yscale("log")
a.set_title("(a) Black hole on a circle (BTZ): clock holds\nuntil the region nearly wraps the horizon"); a.legend(fontsize=8)
a = ax[1]
for al, col in [("0.5", PURPLE), ("0.8", GREEN)]:
    rows = np.array(R["cone"][al])
    a.plot(rows[:, 0] / math.pi, rows[:, 1], color=col, lw=2, label=f"true CMI, α={al}")
    a.plot(rows[:, 0] / math.pi, rows[:, 2], color=col, ls="--", lw=1.5, label=f"clock, α={al}")
a.axvline(1, color=GRAY, ls=":"); a.set_xlabel("A+B+C (units of π)"); a.set_ylabel("CMI (units c/3)")
a.set_ylim(0, 4.5); a.set_title("(b) Heavy particle (conical defect): clock holds\nfor exactly half the circle, then the pure-state flip"); a.legend(fontsize=8)
a = ax[2]; wb = np.array(R["window"]["btz"]); wp = np.array(R["window"]["pure"])
a.plot(wb[:, 0], wb[:, 1] / math.pi, "o-", color=RED, lw=2, label="HEAT: black hole, horizon radius r₊ (∝ temperature)")
a.plot(1 - wp[:, 0], wp[:, 1] / math.pi, "s-", color=PURPLE, lw=2, label="SQUEEZE: heavy matter, 1 − α (→1 = threshold)")
a.axhline(2, color=GRAY, ls=":"); a.axhline(1, color=GRAY, ls=":")
a.set_xscale("log"); a.set_xlim(1e-3, 6); a.set_ylim(0.8, 2.1)
a.set_xlabel("knob (heat: r₊; squeeze: 1 − α)"); a.set_ylabel("clock's reach: largest A+B+C (units of π)")
a.set_title("(c) Heat widens the reach toward the whole circle;\nsqueezing keeps it at exactly half"); a.legend(fontsize=8, loc="center right")
save(fig, "fig10_where_the_reading_breaks")

# ---------------- fig11: local quench
q = R["quench"]; rows = np.array(q["rows"]); dens = np.array(q["density"]); N = len(dens); x0 = N // 2
la, lb, lc, vg = q["la"], q["lb"], q["lc"], q["vg"]; start = x0 - lb // 2 - la
fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
a = ax[0]; xs = np.arange(N) - x0
for (s, w, lab, col) in [(start, la, "A", BLUE), (start + la, lb, "B", GRAY), (start + la + lb, lc, "C", GREEN)]:
    a.axvspan(s - x0, s + w - x0, color=col, alpha=.15); a.text(s + w / 2 - x0, dens.max() * 1.02, lab, ha="center", color=col, fontsize=12)
a.plot(xs, dens, color=ORANGE, lw=2); a.set_xlim(-200, 200)
a.set_xlabel("site (0 = where the particle was added)"); a.set_ylabel("extra particle density")
a.set_title("(a) One added fermion splits into two halves\nthat fly apart (time 25: about to enter A and C)")
a = ax[1]
a.plot(rows[:, 0], rows[:, 1], "o-", color=ORANGE, lw=2)
a.axhline(2 * math.log(2), color=RED, ls="--"); a.text(1, 2 * math.log(2) + .05, "2 ln 2 (one shared bit, counted from both sides)", color=RED, fontsize=8)
for edge, lab in [(lb / 2, "halves leave B"), (lb / 2 + la, "halves leave A, C")]:
    a.axvline(edge / vg, color=GRAY, ls=":"); a.text(edge / vg + 1, .4, lab, rotation=90, color=GRAY, fontsize=8)
a.set_xlabel("time"); a.set_ylabel("CMI − ground-state CMI"); a.set_ylim(-.1, 1.7)
a.set_title("(b) The clock reading holds, breaks by 2 ln 2\nwhile the halves sit in A and C, then returns")
a = ax[2]; cases = q["cases"]
a.bar(range(len(cases)), list(cases.values()), color=[RED, BLUE, GRAY]); a.set_xticks(range(len(cases))); a.set_xticklabels(list(cases.keys()))
a.axhline(2 * math.log(2), color=RED, ls="--", lw=1)
for i, v in enumerate(cases.values()): a.text(i, v + .03, f"{v:.3f}", ha="center")
a.set_xlabel("where the (left half, right half) pair sits"); a.set_ylabel("extra CMI")
a.set_title("(c) The extra weight is a point lump:\nonly the A × C placement feels it")
save(fig, "fig11_local_quench")

# ---------------- fig12: dense fermions and tunnelling
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
a = ax[0]; d = np.array(R["dense"])
a.plot(d[:, 3], np.abs(d[:, 1] - d[:, 2]), "o-", color=PURPLE, lw=2)
a.axvline(40, color=GRAY, ls=":"); a.text(41, .022, "block B = 40 sites", color=GRAY, fontsize=8)
a.set_xscale("log"); a.set_xlabel("spacing between holes (sites) as the chain fills from 50% to 99%")
a.set_ylabel("|CMI − clock|"); a.set_title("(a) Squeezed cold fermions: the clock holds until\nthe particle spacing nears the block size")
a = ax[1]; t = np.array(R["tunnel"])
a.plot(t[:, 0], t[:, 1], "o-", color=BLUE, lw=2, label="−ln T for the whole barrier")
a.plot(t[:, 0], t[:, 2], "s--", color=ORANGE, lw=2, label="sum of the two halves' −ln T")
a.set_xlabel("width of each half"); a.set_ylabel("tunnelling exponent")
a.set_title("(b) Tunnelling: exponents add, so 'blocked' glues\nby ∘₋₁ up to a fixed edge term (rel. miss → 2%)"); a.legend(fontsize=8)
save(fig, "fig12_squeezed_fermions_and_tunnelling")

# ---------------- fig13: what never breaks
fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.4))
a = ax[0]; cols = {"heated": RED, "quench": ORANGE, "random qubits": PURPLE}
for nm, col in cols.items():
    g = np.array([(r[1], r[2]) for r in R["glue"] if r[0] == nm])
    a.scatter(g[:, 1], g[:, 0], s=14, color=col, label=f"{nm} ({len(g)})", alpha=.7)
diag(a, 0, 1); a.set_xlabel("x₁ ∘₋₁ x₂  (glued from the two pieces)"); a.set_ylabel("x of the whole")
a.set_title("(a) Gluing never breaks: it is the chain rule\nof conditional mutual information"); a.legend(fontsize=8)
a = ax[1]; names = list(dict.fromkeys(n for n, _ in R["keff"]))
for i, nm in enumerate(names):
    dev = [max(abs(k + 1), 1e-17) for n, k in R["keff"] if n == nm]
    a.scatter(dev, [i] * len(dev), color=BLUE, s=25)
a.set_yticks(range(len(names))); a.set_yticklabels(names, fontsize=8)
a.set_xscale("log"); a.set_xlim(1e-18, 1e-2); a.axvline(1e-6, color=RED, ls="--"); a.text(1.3e-6, len(names) - 1.5, "check\ntolerance", color=RED, fontsize=8)
a.set_xlabel("|k_eff − (−1)|  for each triple of regions (log scale)")
a.set_title("(b) The dial k = −1 in every environment:\nevery miss is rounding error")
save(fig, "fig13_gluing_never_breaks")

# ---------------- fig14: the bath as a bigger cube
fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.6), gridspec_kw={"width_ratios": [1.25, 1, 1]})
a = ax[0]; a.axis("off"); a.set_xlim(0, 12); a.set_ylim(0, 7.2)
def square(ox, oy, s, col, labs):
    a.add_patch(Polygon([(ox, oy), (ox + s, oy), (ox + s, oy + s), (ox, oy + s)], fc=col, ec="k", alpha=.25))
    for (dx, dy), lab in zip([(0, 0), (s, 0), (0, s), (s, s)], labs): a.text(ox + dx, oy + dy, lab, ha="center", va="center", fontsize=8, bbox=dict(fc="white", ec="none", pad=.5))
labs = ["∅", "A", "B", "AB"]
square(0.6, 0.6, 2.2, GRAY, [l + ", off" for l in labs]); square(1.8, 3.6, 2.2, ORANGE, [l + ", on" for l in labs])
for (dx, dy) in [(0, 0), (2.2, 0), (0, 2.2), (2.2, 2.2)]: a.plot([0.6 + dx, 1.8 + dx], [0.6 + dy, 3.6 + dy], color="k", lw=.8, ls=":")
a.text(0.2, 6.6, "system + bath cube", fontsize=10, weight="bold")
a.text(4.3, 5.0, "bath on (weight π):\nA, B independent → ∘₋₁", fontsize=8, color=ORANGE)
a.text(3.4, 1.2, "bath off (1 − π):\nnothing fires", fontsize=8, color=GRAY)
a.annotate("", (8.0, 3.2), (5.9, 3.2), arrowprops=dict(arrowstyle="->", lw=2, color=BLUE))
a.text(6.95, 3.5, "forget the bath\n(sum over on/off)", ha="center", fontsize=8, color=BLUE)
square(8.4, 2.1, 2.2, BLUE, labs); a.text(9.5, 5.0, "system cube", ha="center", fontsize=10, weight="bold")
a.text(9.5, 0.9, "P(A or B) = x_A ∘k′ x_B\nk′ = −1/π", ha="center", fontsize=9, color=BLUE)
a.set_title("(a) The system cube is a face of a bigger cube;\ntracing out the bath changes the dial")
a = ax[1]; T = np.linspace(0.05, 5, 300)
for mu, col in [(0.0, BLUE), (0.8, GREEN), (1.5, RED)]:
    kp = -(1 + np.exp((1.0 - mu) / T)); a.plot(T, kp, color=col, lw=2, label=f"μ = {mu} (Δ = 1)")
a.axhline(-1, color=GRAY, ls=":"); a.axhline(-2, color=GRAY, ls=":")
a.text(0.3, -0.75, "−1: plain union rule (bath always on)", ha="left", fontsize=8, color=GRAY); a.text(4.9, -1.75, "−2: infinitely hot", ha="right", fontsize=8, color=GRAY)
a.set_ylim(-8, -0.5); a.set_xlabel("bath temperature T = 1/β"); a.set_ylabel("dial k′ after forgetting the bath")
a.set_title("(b) k′ = −(1 + exp β(Δ − μ − ΩJ)) = −1/n_Fermi\n(two-state bath; ours, toy)"); a.legend(fontsize=8, loc="lower right")
a = ax[2]
pi_on, Lam, av = 0.3, 0.8, np.array([.5, .4, .7]); x = list(pi_on * Lam * av)
fold = lambda xs, k: __import__("functools").reduce(lambda s, v: circ(s, v, k), xs, 0.0)
bars = [("simulated\n(400k runs)", 0.2459, ORANGE), ("∘k′ with k′ = −1/π", fold(x, -1 / pi_on), BLUE), ("plain ∘₋₁\n(bath ignored)", fold(x, -1), GRAY)]
a.bar(range(3), [b[1] for b in bars], color=[b[2] for b in bars]); a.set_xticks(range(3)); a.set_xticklabels([b[0] for b in bars], fontsize=8)
for i, b in enumerate(bars): a.text(i, b[1] + .005, f"{b[1]:.4f}", ha="center", fontsize=9)
a.set_ylabel("P(at least one of 3 events)")
a.set_title("(c) Check with π = 0.3: the forgotten bath\nglues at k′ = −3.33, not −1")
save(fig, "fig14_bath_as_bigger_cube")
print("ok")
