"""Pictures for the chain-rule test (DRAFT). Run after checks.py.
Colour key: blue = the law's own reading; orange = textbook physics; grey = reference lines."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BLUE, ORANGE, GREY, RED = "#1f77b4", "#ff7f0e", "#7f7f7f", "#c0392b"
R = json.load(open("results.json"))
fig, axs = plt.subplots(1, 2, figsize=(11.5, 4.6))
ax = axs[0]
al = [0.5, 1.0, 2.0, 3.0]
law = [max(R["C2"]["law_chain_error"][str(a)], 1e-17) for a in al]
plain = [max(R["C2"]["plain_chain_error"][str(a)], 1e-17) for a in al]
x = np.arange(len(al))
ax.bar(x - 0.18, law, 0.36, color=BLUE, label="law's chain rule at its own dial (lift, average, unlift)")
ax.bar(x + 0.18, plain, 0.36, color=GREY, label="plain chain rule (Shannon's: weights = probabilities)")
ax.set_yscale("log"); ax.set_ylim(1e-17, 3)
ax.set_xticks(x); ax.set_xticklabels([f"Renyi {a}\nk = {(1 - a) * np.log(2):+.2f}" for a in al], fontsize=8)
ax.set_ylabel("chain-rule error (bits)")
ax.text(1, 3e-15, "Shannon:\nboth agree", ha="center", fontsize=8, color=ORANGE)
ax.legend(fontsize=7, loc="upper left")
ax.set_title("A · One chain rule per dial (C2): the law holds them all;\nonly k = 0 uses plain probabilities")
ax = axs[1]
E = np.logspace(-3, -0.05, 200)
def info(E, a):
    P1, P2 = (1 + E) / 2, (1 - E) / 2
    if a == np.inf: return 1 + np.log2(P1)
    if a == 1: return ((1 + E) * np.log1p(E) + (1 - E) * np.log1p(-E)) / (2 * np.log(2))
    return 1 - np.log2(P1 ** a + P2 ** a) / (1 - a)
for a, col in [(0.5, "#9cc6ea"), (1, ORANGE), (2, BLUE), (10, "#0b4f8a"), (np.inf, RED)]:
    lab = "guessing (alpha = inf): slope 1, local 2" if a == np.inf else (f"Shannon: slope 2" if a == 1 else f"Renyi {a}: slope 2")
    ax.loglog(E, [info(e, a) for e in E], color=col, label=lab)
ax.set_xlabel("bias E of a coin"); ax.set_ylabel("information (bits)")
ax.legend(fontsize=7); ax.set_title("B · Every finite dial is smooth (C3): same circle,\nsame 2 sqrt 2; only guessing differs")
fig.savefig("figures/cr1_dial.png", bbox_inches="tight")
print("figure written")
