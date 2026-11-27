"""Figures for the arXiv draft. Run: python3 paper_figures.py (writes figures/*.pdf)."""
import os
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "figures"); os.makedirs(OUT, exist_ok=True)
GREEN, RED, GREY = "#2f855a", "#c53030", "#718096"
plt.rcParams.update({"font.size": 10, "axes.titlesize": 10, "mathtext.fontset": "cm"})

# Fig. 1: residual of the Planck-Einstein condition under a shift of the energy zero
J2 = np.array([[0., -1.], [1., 0.]]); Jstd = np.kron(np.eye(2), J2)
P = np.eye(4)[[0, 2, 1, 3]]; JB = P @ Jstd @ P.T
H = P @ np.diag([1., 1., 2., 2.]) @ P.T            # H commutes with JB, not with Jstd
cs = np.linspace(-3, 3, 241); I4 = np.eye(4)
def resid(B, Jshift):
    return [np.linalg.norm((B + c*Jshift) @ (B + c*Jshift) + (H + c*I4) @ (H + c*I4)) for c in cs]
fig, a = plt.subplots(figsize=(4.6, 3.1))
a.plot(cs, resid(JB @ H, JB), color=GREEN, lw=2.2, label=r"$J_H = J$ (shared quarter turn)")
a.plot(cs, resid(JB @ H, Jstd), color=RED, lw=2.0, ls="--", label=r"$J_H \neq J$")
a.set_xlabel(r"shift of the energy zero $c$")
a.set_ylabel(r"$\|(B+cJ)^2+(H+c)^2\|$")
a.legend(fontsize=8.5, frameon=False); a.set_ylim(-1, 29)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig1-energy-zero.pdf")); plt.close(fig)

# Fig. 2: joint parameters versus products of local observables for two two-level systems
labels = ["real tensor\nproduct", "complex (balanced)\nproduct", "quaternionic"]
joint = [10, 16, 28]; local = [9, 16, 36]
x = np.arange(3)
fig, a = plt.subplots(figsize=(4.6, 3.1))
a.bar(x - .18, joint, .34, color=GREY, label="joint observables")
a.bar(x + .18, local, .34, color=[RED, GREEN, RED], label="products of local observables")
for xi, j, l in zip(x, joint, local):
    a.text(xi - .18, j + .6, str(j), ha="center", fontsize=9); a.text(xi + .18, l + .6, str(l), ha="center", fontsize=9)
a.set_xticks(x); a.set_xticklabels(labels, fontsize=8.5); a.set_ylim(0, 44)
a.legend(fontsize=8.5, frameon=False, loc="upper left"); a.set_ylabel("real dimension")
fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig2-local-tomography.pdf")); plt.close(fig)
print("paper figures written")
