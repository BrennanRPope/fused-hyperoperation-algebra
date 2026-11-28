"""Picture for the information-causality test (DRAFT). Run after checks.py.
Colour key: blue = the law's own reading; orange = textbook physics; grey = reference lines."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BLUE, ORANGE, GREY, RED = "#1f77b4", "#ff7f0e", "#7f7f7f", "#c0392b"
R = json.load(open("results.json"))
fig, axs = plt.subplots(1, 2, figsize=(11.5, 4.8))
ax = axs[0]
ax.fill([0, 1, 0], [0, 0, 1], color="#dddddd", label="local boxes: E1 + E2 <= 1 (CHSH <= 2)")
ph = np.linspace(0, np.pi / 2, 200)
ax.fill(np.r_[0, np.cos(ph)], np.r_[0, np.sin(ph)], color=ORANGE, alpha=0.18, label="allowed by information causality")
ax.plot(np.cos(ph), np.sin(ph), color=BLUE, lw=2.5, label="quantum boxes = e^(i phi): the law's rotation lift")
ax.plot([1], [1], "s", color=RED); ax.text(0.86, 1.04, "PR box, CHSH 4", fontsize=8, color=RED)
s = 1 / np.sqrt(2)
ax.plot([s], [s], "o", color=BLUE); ax.text(s + 0.03, s + 0.03, "2 sqrt 2 (Tsirelson)", fontsize=8, color=BLUE)
ax.plot([0, 1.35], [1.35 * 0 + np.sqrt(2), np.sqrt(2) - 1.35], color=GREY, ls=":", lw=1)
ax.text(0.86, 0.3, "line of\nconstant CHSH\ntouches the circle\nat 45 degrees", fontsize=7, color=GREY)
ax.set_xlim(0, 1.15); ax.set_ylim(0, 1.15); ax.set_aspect("equal")
ax.set_xlabel("bias E1 (Bob's question 0)"); ax.set_ylabel("bias E2 (Bob's question 1)")
ax.legend(fontsize=7, loc="lower left")
ax.set_title("A · Where boxes live (I2, I4)")
ax = axs[1]
cols = {"1.0": RED, "0.75": "#e67e22", "0.72": "#f5b041", str(1 / np.sqrt(2)): BLUE, "0.6": GREY}
for key, ys in R["I3"]["curves"].items():
    r = float(key)
    lab = f"r = {r:.4f}" + (" (on the circle)" if abs(r - 1 / np.sqrt(2)) < 1e-9 else "")
    ax.plot(range(1, len(ys) + 1), ys, "-o", ms=2.5, color=cols.get(key, GREY), label=lab)
ax.axhline(1, color="k", ls="--", lw=1); ax.text(30, 1.05, "1 bit sent", fontsize=8)
ax.set_ylim(0, 2.2); ax.set_xlabel("levels n in the tree (Alice holds 2^n bits)")
ax.set_ylabel("Bob's information about Alice's bits (bits)")
ax.legend(fontsize=7, title="box with E1 = E2 = r", title_fontsize=7)
ax.set_title("B · Information causality (I3): boxes past the circle\nlet one sent bit carry more than one bit")
fig.savefig("figures/ic1_disc.png", bbox_inches="tight")
print("figure written")
