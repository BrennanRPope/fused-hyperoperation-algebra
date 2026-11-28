import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
out = "/mnt/project-files/math-theory/bridges-study/c2-one-body/figures/c2e_one_body_test.png"
fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
# panel a: axis moments |M (alpha)^l| for M=1, alpha=0.6
M, al = 1.0, 0.6
l = np.arange(6)
exact = M*al**l
rules = np.where(l <= 1, exact, 0)
w = 0.38
ax[0].bar(l - w/2, exact, w, color="#2f6db5", label="moved mass (river's exact shift by i·a)")
ax[0].bar(l + w/2, rules, w, color="#e08a2c", label="Part 8 moving-body rules (1PN)")
for li in l[2:]:
    ax[0].text(li + w/2, 0.01, "none", ha="center", va="bottom", fontsize=8, color="#e08a2c", rotation=90)
ax[0].set_xticks(l); ax[0].set_xlabel("moment l (0 = mass, 1 = spin/dipole, 2+ = shape)")
ax[0].set_ylabel("size of moment  M·α^l   (M=1, α=0.6)")
ax[0].set_title("a. The moved mass needs every moment;\nthe 1PN rules stop after the first step", fontsize=10)
ax[0].legend(fontsize=8, frameon=False)
# panel b: turn per loop for each shell, beta*OmE(r)/2pi
u = np.arcsinh(al/M)
rp = np.sqrt(M**2 + al**2) + M
rm = M - np.sqrt(M**2 + al**2)
beta = 4*np.pi*(rp**2 - al**2)/(rp - rm)
r = np.linspace(rp, 8, 300)
turn = beta*al/(r**2 - al**2)/(2*np.pi)
ax[1].plot(r, turn, color="#2f6db5")
ax[1].axhline(np.tanh(u), ls="--", color="gray", lw=0.8)
ax[1].plot([rp], [np.tanh(u)], "o", color="#c0392b")
ax[1].annotate("horizon shell: turn per loop = tanh u\n(the Euclidean Kerr helix)", (rp, np.tanh(u)),
               xytext=(rp+1.2, np.tanh(u)*0.8), fontsize=8, arrowprops=dict(arrowstyle="->", lw=0.7))
ax[1].text(5.2, 0.12, "far shells turn less;\nonly the horizon\nsets the loop length β", fontsize=8, color="#555")
ax[1].set_xlabel("shell radius r  (M=1, α=0.6)"); ax[1].set_ylabel("turns per imaginary-time loop")
ax[1].set_ylim(0, 0.7)
ax[1].set_title("b. The river's rigid shells, run in imaginary time,\ngive the helix at the horizon", fontsize=10)
fig.suptitle("DRAFT · One spinning body: what the river's rules give", fontsize=11)
fig.tight_layout()
fig.savefig(out, dpi=150)
print("saved", out, "tanh u =", round(np.tanh(u),4))
