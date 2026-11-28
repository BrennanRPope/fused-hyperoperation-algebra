"""Picture for the forces-circle test (DRAFT). Run after checks.py.
Colour key: blue = the law's own reading; orange = textbook physics; grey = reference lines."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BLUE, ORANGE, GREY, RED = "#1f77b4", "#ff7f0e", "#7f7f7f", "#c0392b"
fig, axs = plt.subplots(1, 2, figsize=(11.5, 4.8))
ax = axs[0]
t = np.linspace(0, np.pi / 2, 300)
for av, col, lab in [(1, GREY, "a = 1: total variation, guessing (CHSH 2)"), (2, BLUE, "a = 2: Shannon, any smooth measure (2 sqrt 2)"),
                     (4, "#e67e22", "a = 4 (3.36)"), (60, RED, "a -> infinity (PR box, 4)")]:
    c, s = np.cos(t), np.sin(t)
    rr = (np.abs(c) ** av + np.abs(s) ** av) ** (-1 / av)
    ax.plot(rr * c, rr * s, color=col, lw=2.2 if av == 2 else 1.3, label=lab)
ax.set_aspect("equal"); ax.set_xlim(0, 1.1); ax.set_ylim(0, 1.1)
ax.set_xlabel("bias E1"); ax.set_ylabel("bias E2")
ax.legend(fontsize=7, loc="lower left")
ax.set_title("A · Boxes allowed when information grows like |E|^a (F1)")
ax = axs[1]
a = np.linspace(0.8, 12, 400)
ax.plot(a, 2 ** (2 - 1 / a), color=BLUE)
for y, txt in [(2, "local 2"), (2 * np.sqrt(2), "quantum 2 sqrt 2"), (4, "PR box 4")]:
    ax.axhline(y, color=GREY, ls=":", lw=1); ax.text(9, y + 0.04, txt, fontsize=8, color=GREY)
ax.plot([1, 1], [2, 2], "o", color=GREY); ax.text(1.15, 1.85, "guessing, total variation", fontsize=8)
ax.plot([2], [2 * np.sqrt(2)], "o", color=ORANGE, ms=8)
ax.text(2.2, 2.65, "Shannon, Renyi-2, chi^2:\nsmooth + flip-symmetric", fontsize=8, color=ORANGE)
ax.set_xlabel("exponent a: information ~ |bias|^a near a fair coin"); ax.set_ylabel("CHSH ceiling")
ax.set_title("B · The ceiling is set by one number (F1, F2)")
fig.savefig("figures/fc1_exponent.png", bbox_inches="tight")
print("figure written")
