"""DRAFT figure O8: where gravity's 2v comes from. Run after checks_twov.py."""
import json, math, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "..", "figures")
SPACE, SR, EN = "#888888", "#d9822b", "#2b6cb0"
plt.rcParams.update({"font.size": 11, "axes.titlesize": 12})
fig, axs = plt.subplots(1, 3, figsize=(17, 5.2), gridspec_kw={"width_ratios": [1, 1.1, 1.1]})
ax = axs[0]
cases = ["electric orbit", "GR gyroscope", "GR orbit"]
space = [0, 2, 2]; sr = [1, 1, 1]; en = [0, 0, 3]
ax.bar(cases, space, color=SPACE, label="curved space (1/B)")
ax.bar(cases, sr, bottom=space, color=SR, label="special relativity (Thomas disc, radius v)")
ax.bar(cases, en, bottom=[a + b for a, b in zip(space, sr)], color=EN, label="gravity pulls on energy (3v²)")
for i, t in enumerate([1, 3, 6]): ax.text(i, t + 0.12, str(t), ha="center", fontsize=13, fontweight="bold")
ax.set_ylim(0, 7.3); ax.set_ylabel("turn per lap ÷ (hodograph area / c²), first order"); ax.legend(fontsize=8.5, loc="upper left")
ax.set_title("A. The 6 is 2 + 1 + 3: the '2v' is\n1 + 3 = 4, not a doubling")
ax = axs[1]
rr = np.linspace(2.6, 20, 400)
ax.plot(rr, -1 / rr, color="black", lw=1.5, label="Newton: −M/r")
ax.plot(rr, 1 - 1 / np.sqrt(1 - 2 / rr), color=EN, lw=2.2, label="what the orbit feels: E(1 − 1/√(1 − 2M/r)), per unit E")
ax.set_ylim(-0.9, 0.02); ax.set_xlabel("r / M"); ax.set_ylabel("potential energy per unit energy")
ax.text(6.5, -0.55, "steeper than 1/r: a falling body's\nenergy grows as 1/√(1 − 2M/r),\nand gravity pulls on that energy", color=EN, fontsize=9.5)
ax.legend(fontsize=8.5, loc="lower right")
ax.set_title("B. Inside the clock warp, the pull is steeper\nthan Newton's, which spoils the closure by 3v²")
ax = axs[2]
x = np.linspace(0, 0.2, 300); v2 = x / (1 - 2 * x)
ax.plot(x, 1 - 2 * x, color=SPACE, lw=2.2, label="1 − 2M/r  (space)")
ax.plot(x, 1 - v2, color=SR, lw=2, label="1 − v²  (special relativity)")
ax.plot(x, 1 - 3 * v2, color=EN, lw=2, label="1 − 3v²  (energy pull)")
ax.plot(x, (1 - 2 * x) * (1 - 4 * v2), color="black", lw=1.5, ls="--", label="(1 − 2M/r)(1 − v² − 3v²) = 1 − 6M/r")
ax.axvline(1 / 6, color=EN, lw=0.8, ls=":"); ax.text(1 / 6 + 0.002, 0.8, "r = 6M", color=EN, fontsize=9)
ax.axhline(0, color="black", lw=0.5)
ax.set_xlabel("M/r"); ax.set_ylabel("factor"); ax.set_ylim(-0.4, 1.05); ax.legend(fontsize=8.5, loc="lower left")
ax.set_title("C. Exact at every field strength: the\nspeed parts add (v² + 3v²), space multiplies")
fig.suptitle("Figure O8. Where gravity's '2v' comes from: special relativity's v² (as for the electric orbit) plus 3v² because gravity pulls on energy", y=1.02)
fig.savefig(os.path.join(OUT, "figO8-where-2v-comes-from.png"), dpi=140, bbox_inches="tight"); print("wrote figO8")
