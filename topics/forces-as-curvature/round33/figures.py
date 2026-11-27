"""DRAFT round 33 picture: f40_weak_pattern.png, from results.json."""
import json, os, math
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
R = json.load(open(os.path.join(HERE, "results.json")))
GREY, BLUE, ORANGE, RED, PURPLE, GREEN = "#777777", "#1f6fc5", "#e08a1e", "#c8372d", "#7b4fa8", "#2e9e4f"
fig = plt.figure(figsize=(20, 13.5))
gs = fig.add_gridspec(2, 2, hspace=0.28, wspace=0.16)
fig.suptitle("Round 33 (DRAFT): taking the hidden 2-sphere as a premise, does the weak force's left-handed pattern follow?", fontsize=16)

# A: sectors
ax = fig.add_subplot(gs[0, 0])
for r in R["W1"]["rows"]:
    w, n = r["weight"], r["left"] + r["right"]
    col = BLUE if r["left"] else ORANGE
    ax.bar(w, n if r["left"] else -n, width=0.38, color=col, alpha=0.85)
    lab = {0.0: "singlet", 0.5: "doublet", 1.0: "triplet"}.get(r["spin"], f"spin {r['spin']:g}")
    ax.text(w, (n + 0.3) if r["left"] else -(n + 0.3), lab, ha="center", va="bottom" if r["left"] else "top", fontsize=9, color=col, rotation=90)
ax.axhline(0, color=GREY, lw=1)
ax.set_xlabel("weight w along the hidden 3-sphere's circles  (flux through the 2-sphere = 2w)")
ax.set_ylabel("massless modes: left-handed up, right-handed down")
ax.set_title("A. Under the premise each weight gives one multiplet, and its size fixes its weight (W1)")
ax.text(-3.2, 7.6, "size and weight are tied: |w| = spin + ½\ndoublets always have |w| = 1, singlets always ½", color=PURPLE, fontsize=11)
ax.set_ylim(-9.5, 9.5)

# B: the tie against the measured hypercharges
ax = fig.add_subplot(gs[0, 1])
SM = [("Q", 1.0, 1/6), ("L", 1.0, 1/2), ("u^c", 0.5, 2/3), ("d^c", 0.5, 1/3), ("e^c", 0.5, 1.0)]
xs = np.linspace(0, 1.7, 10)
for cc in (1/6, 1/2, 1.0):
    ax.plot(xs, cc * xs, color=PURPLE, lw=1.2, ls="--", alpha=0.7)
ax.text(0.05, 1.17, "premise: |Y| = |c| × (spin + ½), one line through 0 for every field", color=PURPLE, fontsize=11, va="top", bbox=dict(fc="white", ec="none", alpha=0.9))
for name, x, y in SM:
    ax.plot(x, y, "o", ms=13, color=RED if x == 1.0 else ORANGE)
    ax.annotate(name, (x, y), (x + 0.04, y + 0.03), fontsize=12, color=RED if x == 1.0 else ORANGE)
ax.annotate("", (1.0, 0.5), (1.0, 1 / 6), arrowprops=dict(arrowstyle="<->", color=RED, lw=1.5))
ax.text(1.05, 0.30, "two doublets,\ntwo hypercharges:\nno one line can hold both", color=RED, fontsize=10)
ax.set_xlim(0, 1.7); ax.set_ylim(0, 1.22)
ax.set_xticks([0.5, 1.0]); ax.set_xticklabels(["singlets\n(spin + ½ = ½)", "doublets\n(spin + ½ = 1)"])
ax.set_ylabel("measured |hypercharge| of each Standard Model field")
ax.set_title("B. The measured hypercharges break the tie: no line fits (W2)")

# C: anomalies
ax = fig.add_subplot(gs[1, 0])
labels = ["Tr Y", "Tr Y³", "Tr T²·Y", "odd number\nof doublets"]
sm = [0, 0, 0, 0]; r9 = [2, 2, 0.5, 1]
xx = np.arange(4)
ax.bar(xx - 0.18, sm, 0.34, color=GREEN, label="Standard Model generation: all cancel, 4 doublets")
ax.bar(xx + 0.18, r9, 0.34, color=RED, label="premise, round 9's sector alone (c = 1): none cancel, 1 doublet")
for x in xx: ax.text(x - 0.18, 0.05, "0", ha="center", color=GREEN, fontsize=12)
ax.set_xticks(xx); ax.set_xticklabels(labels, fontsize=11)
ax.set_ylabel("anomaly (must be 0 for a consistent theory)")
ax.set_title("C. The premise's matter can't cancel its anomalies: every hypercharge has one sign (W3)")
ax.text(-0.4, 2.55, "under the premise every term is positive,\nso the sums vanish only with no matter at all", ha="left", color=RED, fontsize=11)
ax.legend(loc="upper right", fontsize=10); ax.set_ylim(0, 3.5)

# D: running
ax = fig.add_subplot(gs[1, 1])
MZ = 91.1876; a2, aY = R["W4"]["a2inv"], R["W4"]["aYinv"]; bY, b2 = 41 / 6, -19 / 6
lmu = np.linspace(math.log10(MZ), 21.5, 300); t = np.log(10**lmu / MZ)
ax.plot(lmu, aY - bY / (2 * math.pi) * t, color=PURPLE, lw=2.5, label="1/α of hypercharge")
ax.plot(lmu, a2 + (-b2) / (2 * math.pi) * t, color=BLUE, lw=2.5, label="1/α of the weak SU(2)")
ax.axvline(math.log10(R["W4"]["MPl"]), color=GREY, ls="--"); ax.text(math.log10(R["W4"]["MPl"]) - 0.2, 82, "Planck mass", rotation=90, color=GREY, ha="right")
ms = math.log10(R["W4"]["mu_star"])
ax.plot(ms, aY - bY / (2 * math.pi) * math.log(R["W4"]["mu_star"] / MZ), "o", ms=13, color=RED)
ax.annotate("equal (sin²θ = ½, what the law's\nround 3-sphere would set) only at\n~5 × 10²⁰ GeV, beyond Planck", (ms, 51), (8.5, 58), color=RED, fontsize=11, bbox=dict(fc="white", ec="none", alpha=0.9),
            arrowprops=dict(arrowstyle="->", color=RED))
ax.text(2.3, 40, "measured at the Z: sin²θ = 0.231", color=GREY, fontsize=11)
ax.set_xlabel("log₁₀ (energy / GeV)"); ax.set_ylabel("inverse coupling (one-loop Standard Model running)")
ax.set_title("D. Even granting hypercharge = circle weight, the couplings miss (W4)")
ax.legend(loc="center left", fontsize=10)
os.makedirs(os.path.join(HERE, "figures"), exist_ok=True)
fig.savefig(os.path.join(HERE, "figures", "f40_weak_pattern.png"), dpi=120, bbox_inches="tight")
print("saved")
