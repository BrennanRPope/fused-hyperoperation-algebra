"""DRAFT picture for round 21. Run after checks.py. Blue = space (the law turns), orange = time (the law boosts), grey = light / reference."""
import json, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

OR, BL, GR = "#e07b00", "#1f6fc2", "#777777"
R = json.load(open("results.json"))
plt.rcParams.update({"font.size": 10.5, "figure.dpi": 130})
fig, ax = plt.subplots(2, 2, figsize=(16, 13)); ax = ax.ravel()

a = ax[0]
K = R["kinds"]
per = np.array(K["periodic"]); ang = 2 * np.arctan(per)
ang = np.unwrap(np.r_[ang[:12], np.pi, ang[12:]]) if len(ang) < 25 else np.unwrap(ang)
a.plot(np.arange(len(ang)), np.degrees(ang) / 360, "o-", color=BL, ms=4, label="c > 0 (space): turns of the circle, back to start after 24 kicks")
hyp = np.array(K["boost"]); a.plot(np.arange(len(hyp)), 2 * np.arctanh(hyp) / 18, "s-", color=OR, ms=4, label="c < 0 (time): rapidity / 18, grows for ever")
a.plot(np.arange(len(hyp)), hyp, "s--", color=OR, ms=3, alpha=0.5, label="c < 0 (time): chart position, never reaches the edge 1")
par = np.array(K["plain"]); a.plot(np.arange(len(par)), par / 4, "^-", color=GR, ms=4, label="c = 0 (light): plain addition / 4")
a.axhline(1, color=GR, lw=0.7, ls=":")
a.set_xlabel("number of equal kicks composed with the law"); a.set_ylabel("value (scaled to share one axis)")
a.legend(fontsize=8.5, loc="upper left"); a.set_ylim(0, 1.45)
tk = K["trace2_over_det"]
a.text(0.97, 0.05, f"trace²/det of the law's matrix:\n{tk['1.0']:.2f} (< 4, turn), 4 (plain), {tk['-1.0']:.2f} (> 4, boost)",
       transform=a.transAxes, ha="right", va="bottom", fontsize=8.5, color=GR)
a.set_title("A. Along any direction the law either turns,\nboosts or adds plainly: c = q·η(v, v) decides")

b = ax[1]
th = np.linspace(0, 2 * np.pi, 721)
for t0, t1 in zip(th[:-1], th[1:]):
    tm = (t0 + t1) / 2; v = np.array([math.sin(tm), math.cos(tm)])          # (time, space) components
    c = -v[0]**2 + v[1]**2
    col = BL if c > 1e-9 else OR
    b.fill([0, 1.2 * math.cos(t0), 1.2 * math.cos(t1)], [0, 1.2 * math.sin(t0), 1.2 * math.sin(t1)], color=col, alpha=0.35, lw=0)
for sg in (1, -1):
    b.plot([-1.25, 1.25], [-1.25 * sg, 1.25 * sg], color=GR, lw=2.2)
b.text(0.0, 0.85, "time directions:\nthe law boosts (c < 0)", ha="center", color=OR, fontsize=9)
b.text(0.78, 0.0, "space:\nthe law\nturns\n(c > 0)", ha="center", va="center", color=BL, fontsize=9)
b.text(0.62, 1.12, "light: c = 0,\nplain addition", color=GR, fontsize=9, bbox=dict(fc="w", ec="none", alpha=0.8))
b.set_xlim(-1.3, 1.3); b.set_ylim(-1.3, 1.3); b.set_aspect("equal"); b.set_xlabel("space component of v"); b.set_ylabel("time component of v")
b.text(0.02, -0.2, f"with an all-plus dot every direction turns (smallest c {R['light']['min_c_euclid']:.3f} > 0): no light;\n"
       f"with one minus sign the law along the cone is plain addition to {R['light']['null']:.0e}", transform=b.transAxes, fontsize=8.5, color=GR)
b.set_title("B. One minus sign splits directions into turning\n(space), boosting (time) and light between")

c = ax[2]
tt = np.linspace(0, 2 * np.pi, 200)
c.plot(np.cos(tt), np.sin(tt), color=OR, lw=2.5)
for k in range(8):
    p = np.array([math.cos(k * math.pi / 4), math.sin(k * math.pi / 4)]); v = np.array([-p[1], p[0]])
    c.add_patch(FancyArrowPatch(p, p + 0.35 * v, arrowstyle="-|>", mutation_scale=12, color=OR))
c.text(0, 1.55, "two time directions: a circle stays\ntime-like all the way round (a closed time loop)", ha="center", va="bottom", fontsize=8.5, color=OR)
c.set_xlabel("time 1  |  space"); c.set_ylabel("time 2  |  time")
x2 = np.linspace(-1, 1, 200)
c.plot(2.9 + 0.8 * np.cos(tt), np.sin(tt), color=BL, lw=2.5)
c.plot([2.9 - 0.8, 2.9 + 0.8], [0.62, -0.62], alpha=0)
for yy in (-1, 1):
    c.plot(2.9, yy, "o", color="k", ms=6)
c.text(2.9, -1.95, "one time direction (up): a loop must turn\nback in time where it is flat (dots),\nand there it is not time-like", ha="center", va="bottom", fontsize=8.5, color=BL)
c.set_xlim(-1.6, 4.3); c.set_ylim(-2.0, 2.1); c.set_aspect("equal")
lp = R["loops"]
c.text(0.0, -0.17, "closed loops that stay time-like all the way round: none with 0 or 1 time direction\n(4,000 random loops each); the circle above with 2, 3 or 4", transform=c.transAxes, fontsize=8.5, color=GR, va="top")
c.set_title("C. With two or more time directions time can\nloop; with exactly one it is ordered")

d = ax[3]
T = R["table"]
needs = [("space_closes", "space closes\n(rounds 12–13)"), ("time_ordered", "an ordered\ndirection (time)"), ("light", "light ignores\nthe dial"), ("no_time_loops", "no time\nloops"), ("all", "all four")]
M = np.array([[1.0 if r[k] else 0.0 for k, _ in needs] for r in T])
d.imshow(M, cmap="Blues", vmin=-0.3, vmax=1.3, aspect="auto")
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        d.text(j, i, "yes" if M[i, j] else "no", ha="center", va="center", fontsize=9, color="w" if M[i, j] else "k", weight="bold" if j == 4 and M[i, j] else None)
d.set_xticks(range(len(needs))); d.set_xticklabels([n for _, n in needs], fontsize=9)
d.set_yticks(range(5)); d.set_yticklabels([f"{r['dims_time']} time + {r['dims_space']} space" for r in T])
d.set_ylabel("number of minus signs in the dot product (q > 0)")
d.text(0.0, -0.17, "flipping q and the dot product together leaves the law unchanged,\nso '3 time + 1 space at q < 0' is the same as '1 time + 3 space'", transform=d.transAxes, fontsize=8.5, color=GR)
d.set_title("D. Only one time and three space\ndirections meet all four needs")

fig.tight_layout()
fig.savefig("figures/f28_why_lorentzian.png", bbox_inches="tight")
