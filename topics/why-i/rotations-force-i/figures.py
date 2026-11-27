"""Pictures for topics/why-i/rotations-force-i (Nov 26, 2026). Run: python3 figures.py"""
import numpy as np, matplotlib, os
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Circle, FancyBboxPatch
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "figures"); os.makedirs(OUT, exist_ok=True)
BLUE, ORANGE, GREEN, RED, GREY, PURPLE = "#2b6cb0", "#dd6b20", "#2f855a", "#c53030", "#718096", "#6b46c1"
plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "figure.dpi": 130})
def save(fig, name):
    fig.tight_layout(); fig.savefig(os.path.join(OUT, name + ".png"), bbox_inches="tight"); plt.close(fig)
def arrow(ax, x0, y0, x1, y1, col, lw=2, style="-|>"):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle=style, mutation_scale=14, color=col, lw=lw))

# r1. a rotation needs a plane
fig, ax = plt.subplots(1, 3, figsize=(13, 4.4))
a = ax[0]; a.set_aspect("equal"); a.axis("off"); a.set_xlim(-1.6, 1.6); a.set_ylim(-1.4, 1.4)
a.plot([-1.5, 1.5], [0, 0], color=GREY, lw=1.5)
arrow(a, 0, 0, 1, 0, BLUE); arrow(a, 0, 0, -1, 0, RED, lw=1.5)
a.text(0, 0.35, "a single energy level on a line:\nthe state can only stay, stretch or flip", ha="center", fontsize=9)
a.text(0, -0.55, "no steady turn possible", ha="center", color=RED, fontsize=10)
a.set_title("1. one real level: frozen")
a = ax[1]; a.set_aspect("equal"); a.axis("off"); a.set_xlim(-1.6, 1.6); a.set_ylim(-1.4, 1.4)
a.add_patch(Circle((0, 0), 1, fill=False, color=GREY, lw=1.2))
for th, al in [(0.0, .25), (0.5, .45), (1.0, .7), (1.5, 1.0)]:
    arrow(a, 0, 0, np.cos(th), np.sin(th), BLUE, lw=2*al + .5)
a.annotate("", xy=(np.cos(1.9)*1.2, np.sin(1.9)*1.2), xytext=(1.2, 0.05), arrowprops=dict(arrowstyle="->", color=PURPLE, connectionstyle="arc3,rad=0.35"))
a.text(0.55, 1.15, "turns at rate E/hbar", color=PURPLE, fontsize=9)
a.text(0, -1.3, "a level that is a plane can turn steadily", ha="center", fontsize=9)
a.set_title("2. a plane: the state turns")
a = ax[2]; a.set_aspect("equal"); a.axis("off"); a.set_xlim(-1.6, 1.6); a.set_ylim(-1.4, 1.4)
arrow(a, 0, 0, 1, 0, BLUE); arrow(a, 0, 0, 0, 1, GREEN); arrow(a, 0, 0, -1, 0, ORANGE)
a.text(1.02, -0.18, "u", color=BLUE); a.text(0.05, 1.02, "Ju", color=GREEN); a.text(-1.25, -0.18, "JJu = -u", color=ORANGE)
a.text(0, -0.8, "the quarter turn J of that plane:\nJ squared is -1, the job of i", ha="center", fontsize=9)
a.set_title("3. the turn is a J")
fig.suptitle("Steady turning needs every energy level to be a plane, and each plane brings a quarter turn J (known: Aste 2019, Vaughan 2026)", fontsize=11)
save(fig, "r1-a-rotation-needs-a-plane")

# r2. the energy-zero lever: shared J vs each Hamiltonian with its own J
J2 = np.array([[0., -1.], [1., 0.]]); Jstd = np.kron(np.eye(2), J2)
P = np.eye(4)[[0, 2, 1, 3]]; JB = P @ Jstd @ P.T; HB = P @ np.diag([1., 1., 2., 2.]) @ P.T
cs = np.linspace(-3, 3, 121)
def resid(Bfun, H, Jshift):
    return [np.linalg.norm((Bfun(H) + c*Jshift) @ (Bfun(H) + c*Jshift) + (H + c*np.eye(4)) @ (H + c*np.eye(4))) for c in cs]
fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.6))
a = ax[0]; a.axis("off"); a.set_xlim(0, 10); a.set_ylim(0, 6)
for x, lab, col in [(1.6, "Hamiltonian H:\nturns with J_H", BLUE), (5, "shift the zero:\nH + c adds a\ncommon turn c J", PURPLE), (8.4, "steady turning (A3)\nfor every c forces\nJ_H = J", GREEN)]:
    a.add_patch(FancyBboxPatch((x - 1.45, 2.2), 2.9, 1.6, boxstyle="round,pad=0.05", fc="white", ec=col, lw=2))
    a.text(x, 3.0, lab, ha="center", va="center", fontsize=9)
arrow(a, 3.1, 3.0, 3.5, 3.0, "k"); arrow(a, 6.5, 3.0, 6.9, 3.0, "k")
a.text(5, 1.2, "the leftover J_H - J squares to zero,\nand a rotation generator whose square is zero is zero", ha="center", fontsize=9, color=GREY)
a.text(5, 5.2, "the proof in one line (ours, assembled from known pieces)", ha="center", fontsize=10, weight="bold")
a = ax[1]
a.plot(cs, resid(lambda H: JB @ H, HB, JB), color=GREEN, lw=2.4, label="shared turn: H's own J is also the shift's J")
a.plot(cs, resid(lambda H: JB @ H, HB, Jstd), color=RED, lw=2.2, ls="--", label="H turns with its own J, the shift with another")
a.set_xlabel("shift of the energy zero c"); a.set_ylabel("how badly steady turning fails")
a.set_title("check C2/T3: with two different quarter turns,\nsteady turning breaks as soon as the zero moves"); a.legend(fontsize=8.5)
save(fig, "r2-the-energy-zero-lever")

# r3. composites: booking the shift to either part
fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.6))
a = ax[0]; a.axis("off"); a.set_xlim(0, 10); a.set_ylim(0, 6)
a.add_patch(FancyBboxPatch((0.4, 3.4), 3.2, 1.6, boxstyle="round,pad=0.05", fc="white", ec=BLUE, lw=2)); a.text(2, 4.2, "part A\nturns with J_A", ha="center", va="center")
a.add_patch(FancyBboxPatch((6.4, 3.4), 3.2, 1.6, boxstyle="round,pad=0.05", fc="white", ec=ORANGE, lw=2)); a.text(8, 4.2, "part B\nturns with J_B", ha="center", va="center")
a.add_patch(FancyBboxPatch((2.6, 0.5), 4.8, 1.8, boxstyle="round,pad=0.05", fc="white", ec=GREEN, lw=2))
a.text(5, 1.4, "add c to the TOTAL energy:\nbooked to A it turns by J_A, booked to B by J_B;\nsame physics, so J_A x 1 = 1 x J_B", ha="center", va="center", fontsize=9)
arrow(a, 2.2, 3.3, 3.8, 2.4, BLUE); arrow(a, 7.8, 3.3, 6.2, 2.4, ORANGE)
a.text(5, 5.5, "one shared i for the whole world", ha="center", weight="bold")
a = ax[1]
labels = ["real product\n(two rebits)", "shared turn\n(two qubits)", "quaternionic\n(two quaterbits)"]
joint = [10, 16, 28]; seen = [9, 16, 36]
x = np.arange(3)
a.bar(x - .18, joint, .34, color=GREY, label="joint parameters")
a.bar(x + .18, seen, .34, color=[RED, GREEN, RED], label="products of local observables")
for xi, j, s_ in zip(x, joint, seen): a.text(xi - .18, j + .6, str(j), ha="center"); a.text(xi + .18, s_ + .6, str(s_), ha="center")
a.set_xticks(x); a.set_xticklabels(labels); a.legend(fontsize=8.5, loc="upper left"); a.set_ylim(0, 48)
a.text(0, 13.5, "1 joint turn\nhidden", ha="center", fontsize=8.5, color=RED); a.text(2.18, 40.5, "too many: no\ncomposite fits", ha="center", fontsize=8.5, color=RED)
a.set_title("the shared turn makes the whole equal what the parts see\n(local tomography: check T6; counts as on the e-i-pi page)")
save(fig, "r3-composites-share-one-i")

# r4. real, complex, quaternion: how many quarter turns commute with everything
fig, ax = plt.subplots(figsize=(12, 4.2)); ax.axis("off"); ax.set_xlim(0, 12); ax.set_ylim(0, 5)
cols = [(2, "real numbers", "0 shared quarter turns", "fails A4: no turn commutes with\nevery observable (check T4)", RED),
        (6, "complex numbers", "1 shared quarter turn J", "meets A1-A6:\nordinary quantum mechanics", GREEN),
        (10, "quaternions", "3 quarter turns (i, j, k)", "quaternion-linear dynamics: fails A4 (C3b);\nwith J = R_i instead: meets A1-A4\nbut hides observables (fails A5, C3)", RED)]
for x, t, n, v, col in cols:
    ax.add_patch(FancyBboxPatch((x - 1.8, 0.6), 3.6, 3.6, boxstyle="round,pad=0.05", fc="white", ec=col, lw=2.2))
    ax.text(x, 3.7, t, ha="center", weight="bold", fontsize=11); ax.text(x, 2.9, n, ha="center", fontsize=10)
    ax.text(x, 1.5, v, ha="center", fontsize=9, color=col)
ax.set_title("The only number systems that can carry the turn are R, C and H (Frobenius 1878); the assumptions leave only C")
save(fig, "r4-real-complex-quaternion")

# r5. scorecard of assumptions
rows = [("A1 steps are rotations", "Brennan's premise; Hardy / Masanes-Muller continuous reversibility", "steps could stretch: heat equation, not quantum", BLUE),
        ("A2 energy is conserved", "Noether; Stone", "energy leaks", BLUE),
        ("A3 E = hbar x turning rate", "Planck-Einstein (1900-05)", "rotations unrelated to energy: no J at all", BLUE),
        ("A4 only energy differences matter", "classical mechanics", "each Hamiltonian has its own J: real QM (C2, T4);\nquaternionic QM (C3b)", ORANGE),
        ("A5 nothing hidden", "Masanes-Muller 'all effects allowed'", "a quaternionic variant survives (C3)", ORANGE),
        ("A6 a total shift has no owner", "same as A4, for two parts", "plain real product: hidden joint turn (C4)", ORANGE)]
fig, ax = plt.subplots(figsize=(13.5, 4.2)); ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
for x, h in [(0.02, "assumption"), (0.3, "where it comes from"), (0.62, "what goes wrong without it")]:
    ax.text(x, 0.94, h, weight="bold")
for i, (a_, b_, c_, col) in enumerate(rows):
    y = 0.8 - i*0.14
    ax.add_patch(plt.Rectangle((0.0, y - 0.05), 0.008, 0.1, color=col))
    ax.text(0.02, y, a_, va="center", fontsize=9.5, color=col, weight="bold"); ax.text(0.3, y, b_, va="center", fontsize=9.2); ax.text(0.62, y, c_, va="center", fontsize=9.2)
ax.set_title("What the proof assumes (blue: gives a quarter turn; orange: makes it one shared i)")
save(fig, "r5-assumptions")

# r6. where this sits in the 2021-2026 debate
items = [(2021, "Renou et al.: real QM with the plain real product can be falsified", RED),
         (2022, "Chen et al.; Li et al.: experiments side with complex QM", RED),
         (2025, "Hoffreumon-Woods: real QM with a new combination rule matches everything", BLUE),
         (2026, "Barrios Hita et al. (PRL): same claim, quotient-space construction", BLUE),
         (2026, "Bang-Cho-Baek: those real theories carry a hidden J; Vaughan: evolution needs (at least) a J", GREEN),
         (2026, "this page: rotations + E = hbar w + 'only energy differences matter' force that one shared J", PURPLE)]
fig, ax = plt.subplots(figsize=(12.5, 4.4)); ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
ax.plot([0.07, 0.07], [0.05, 0.92], color="k", lw=1)
for i, (yr, t, col) in enumerate(items):
    y = 0.86 - i*0.15
    ax.plot(0.07, y, "o", color=col, ms=8); ax.text(0.0, y, str(yr), va="center", weight="bold", color=col); ax.text(0.1, y, t, va="center", fontsize=9.8, color=col)
ax.set_title("The debate: real theories that match experiment carry a shared quarter turn; the open question is whether to call it i")
save(fig, "r6-the-debate")
print("figures written")
