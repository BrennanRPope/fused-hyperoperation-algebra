"""Figures for the spin / i tie (DRAFT). Run after checks.py:  python3 figures.py
Colour key: blue = the law, orange = textbook physics (Bartnik-McKinnon, Dirac, Einstein), grey = reference."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from scipy.linalg import expm

BLUE, ORANGE, GREY, RED, GREEN = "#1f6fb4", "#e07b22", "#888888", "#c0392b", "#2e8b57"
R = json.load(open("results.json"))
rng = np.random.default_rng(7)


def box(ax, x, y, w, h, text, color, fs=9.5):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02", fc=color, ec="k", lw=0.8, alpha=0.9))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs, wrap=True)


def arrow(ax, x1, y1, x2, y2, label=None):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1), arrowprops=dict(arrowstyle="->", lw=1.4))
    if label:
        ax.text((x1 + x2) / 2 + 0.01, (y1 + y2) / 2, label, fontsize=8.5, va="center")


# ---------------------------------------------------------------- t1
fig, ax = plt.subplots(figsize=(11, 6))
ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
ax.text(0.25, 0.95, "Why quantum needs i", ha="center", fontsize=13, weight="bold")
ax.text(0.75, 0.95, "Can the twist bubble spin?", ha="center", fontsize=13, weight="bold")
box(ax, 0.07, 0.78, 0.36, 0.1, "shift every energy by a constant c", "#d6e6f5")
box(ax, 0.07, 0.58, 0.36, 0.1, "the state turns by an extra common phase e^(-Jct/ħ)", "#d6e6f5")
box(ax, 0.07, 0.38, 0.36, 0.1, "invisible only if J commutes with everything\n→ one shared J with J² = -1: the i", "#d6e6f5")
arrow(ax, 0.25, 0.78, 0.25, 0.68); arrow(ax, 0.25, 0.58, 0.25, 0.48)
box(ax, 0.57, 0.78, 0.36, 0.1, "rotate the round twist field (f = 1 at its node)", "#fbe3cc")
box(ax, 0.57, 0.58, 0.36, 0.1, "the rotation is a constant gauge turn W = -T₃ (S4)", "#fbe3cc")
box(ax, 0.57, 0.38, 0.36, 0.1, "angular momentum = electric charge of that turn\n→ spinning needs charge (BHSV 1997; vdB–R 2002)", "#fbe3cc")
arrow(ax, 0.75, 0.78, 0.75, 0.68); arrow(ax, 0.75, 0.58, 0.75, 0.48)
box(ax, 0.03, 0.05, 0.94, 0.2,
    "Shared root (textbook): a symmetry that acts only up to an internal U(1) turn\n"
    "hands its conserved charge to that U(1).  Energy ↔ phase rate (E = ħω);  angular momentum ↔ electric charge.\n"
    "Sharper (our reading): a constant shift is a convention exactly when its generator commutes (S1, S3).",
    "#eeeeee", fs=10)
arrow(ax, 0.25, 0.38, 0.4, 0.24); arrow(ax, 0.75, 0.38, 0.6, 0.24)
fig.savefig("figures/t1-the-tie-at-a-glance.png", dpi=140, bbox_inches="tight"); plt.close(fig)

# ---------------------------------------------------------------- t2
n = 3
def realify(M): return np.block([[M.real, -M.imag], [M.imag, M.real]])
Hc = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n)); Hc = (Hc + Hc.conj().T) / 2
H = realify(-1j * Hc); J = realify(1j * np.eye(n))
G = rng.normal(size=(2 * n, 2 * n)); G = G - G.T; G *= np.linalg.norm(J) / np.linalg.norm(G)
cs = np.linspace(0, 3, 61)
dc = [np.abs(expm(H + c * J) - expm(H) @ expm(c * J)).max() + 1e-17 for c in cs]
dn = [np.abs(expm(H + c * G) - expm(H) @ expm(c * G)).max() + 1e-17 for c in cs]
fig, axs = plt.subplots(1, 2, figsize=(12, 4.6))
axs[0].semilogy(cs, dn, color=GREY, lw=2, label="shift generator G does not commute: the shift shows")
axs[0].semilogy(cs, dc, color=BLUE, lw=2, label="shift generator = the shared J: invisible (machine zero)")
axs[0].set_xlabel("size of the constant shift c"); axs[0].set_ylabel("|e^(H + cG) − e^H e^(cG)|  (how visible)")
axs[0].set_title("A. The shift lemma: a constant shift is a convention\nexactly when its generator commutes", fontsize=10.5)
axs[0].legend(fontsize=8.5, loc="center right"); axs[0].set_ylim(1e-17, 1e2)
ax = axs[1]; ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
ax.set_title("B. The two questions use the lemma in opposite directions", fontsize=10.5)
box(ax, 0.02, 0.58, 0.45, 0.32, "Quantum (why-i, A4)\nDemand: the energy zero is a convention\n⇒ the shift generator J is central\n⇒ J² = −1 shared by all: the i", "#d6e6f5", fs=9)
box(ax, 0.53, 0.58, 0.45, 0.32, "Twist bubble (BHSV 1997)\nIn the non-Abelian field the potential\nzero c₀τ_z does NOT commute\n⇒ c₀ is physical (Abelian: convention)", "#fbe3cc", fs=9)
box(ax, 0.53, 0.12, 0.45, 0.32, "Staticity theorem: a regular bubble is\nstatic unless Tr(φ_∞ Q) ≠ 0\n⇒ to spin it needs a physical\npotential zero AND a charge", "#fbe3cc", fs=9)
box(ax, 0.02, 0.12, 0.45, 0.32, "Same algebra, contrapositive:\nwhere the shift is a convention, nothing\nis gained; where it is physical, the\nbubble can store spin in it", "#eeeeee", fs=9)
arrow(ax, 0.75, 0.58, 0.75, 0.44)
fig.tight_layout(); fig.savefig("figures/t2-the-shift-lemma.png", dpi=140, bbox_inches="tight"); plt.close(fig)

# ---------------------------------------------------------------- t3
fig, axs = plt.subplots(1, 2, figsize=(12, 4.6))
f = np.linspace(0, 2, 401)
axs[0].plot(f, np.abs(1 - f), color=BLUE, lw=2.2, label="|D(x̂·T)| / |D(x̂·T)| at f=0  (= |w| = |1 − f|, S3)")
pts = [(0, 3), (1, 1), (2, 3)]
ax2 = axs[0].twinx()
ax2.scatter([0, 1, 2], [3, 1, 3], color=BLUE, s=90, zorder=5, marker="s")
ax2.scatter([0.4, 0.6, 1.4, 1.7], [0, 0, 0, 0], color=GREY, s=50, marker="s")
for x_, y_, t_, off in [(0, 3, "vacuum: 3", (10, -4)), (1, 1, "the law: 1", (10, -4)), (2, 3, "vacuum: 3", (-62, -4))]:
    ax2.annotate(t_, (x_, y_), textcoords="offset points", xytext=off, fontsize=9, color=BLUE)
ax2.annotate("elsewhere: 0", (1.4, 0), textcoords="offset points", xytext=(4, 8), fontsize=9, color=GREY)
ax2.set_ylim(-0.4, 3.6); ax2.set_ylabel("generators whose constant shift is a convention (squares)")
axs[0].set_xlabel("position on the twist line f  (f = 1 − w)"); axs[0].set_ylabel("how visible a shift along x̂ is (line)")
axs[0].set_title("\n\nA. Along the law's twist line, exactly one quarter turn\nsurvives at f = 1: the hedgehog J(x) = 2x̂·T, J² = −1", fontsize=10.5)
axs[0].legend(fontsize=8, loc="lower left", bbox_to_anchor=(0.0, 1.13))
labels = ["real QM", "complex QM", "quaternionic", "f = 0, 2\n(vacua)", "f = 1\n(the law)", "other f"]
vals = [0, 1, 3, 3, 1, 0]
cols = [GREY, ORANGE, GREY, BLUE, BLUE, BLUE]
axs[1].bar(range(6), vals, color=cols, alpha=0.85)
axs[1].axvline(2.5, color="k", lw=0.8, ls="--")
axs[1].text(1, 3.25, "why-i: shared quarter turns", ha="center", fontsize=9)
axs[1].text(4, 3.25, "twist line: surviving quarter turns", ha="center", fontsize=9)
axs[1].set_xticks(range(6)); axs[1].set_xticklabels(labels, fontsize=8.5); axs[1].set_ylim(0, 3.6)
axs[1].set_ylabel("count")
axs[1].set_title("B. The same 0 / 1 / 3 pattern (Frobenius) on both sides.\nOnly 1 gives a single shared i; the law's f = 1 is that case", fontsize=10.5)
fig.tight_layout(); fig.savefig("figures/t3-one-quarter-turn-survives-at-f1.png", dpi=140, bbox_inches="tight"); plt.close(fig)

# ---------------------------------------------------------------- t4
fig, axs = plt.subplots(1, 2, figsize=(12, 4.6))
# slice of the round field: arrows of A^3_i in the z = 0 plane before and after rotation
xs = np.linspace(-1.5, 1.5, 13); X, Y = np.meshgrid(xs, xs)
r = np.sqrt(X ** 2 + Y ** 2) + 1e-9; w = np.cos(1.7 * r) * np.exp(-0.2 * r)
# component along T_3 of A_i: eps_{3ij} x_j (1-w)/r^2 -> (y, -x)(1-w)/r^2
Ux, Uy = Y * (1 - w) / r ** 2, -X * (1 - w) / r ** 2
axs[0].quiver(X, Y, Ux, Uy, color=BLUE, scale=25)
th = np.linspace(0, 2 * np.pi, 200)
for rr in [0.6, 1.2]:
    axs[0].plot(rr * np.cos(th), rr * np.sin(th), color=GREY, lw=0.6, ls=":")
axs[0].set_aspect("equal"); axs[0].set_title("A. The T₃ part of the round twist field in the z = 0 plane.\n"
                                            "Turning space about z = turning the colours by −T₃:\n"
                                            "residual 1.8e-10 against a field change ≥ 0.06 (S4)", fontsize=10)
axs[0].set_xlabel("x"); axs[0].set_ylabel("y")
ax = axs[1]; ax.set_xlim(-1.6, 1.6); ax.set_ylim(-1.6, 1.6); ax.set_aspect("equal"); ax.axis("off")
ax.add_patch(plt.Circle((0, 0), 1.2, fc="none", ec=ORANGE, lw=2))
ax.add_patch(plt.Circle((0, 0), 0.45, fc="#fbe3cc", ec=ORANGE))
ax.text(0, 0, "bubble", ha="center", va="center", fontsize=9)
for a_ in np.linspace(0, 2 * np.pi, 12, endpoint=False):
    ax.annotate("", xy=(1.5 * np.cos(a_), 1.5 * np.sin(a_)), xytext=(1.2 * np.cos(a_), 1.2 * np.sin(a_)),
                arrowprops=dict(arrowstyle="->", color=ORANGE))
ax.text(0, 1.38, "J = ∮ ⟨W, E⟩ dS = charge in the W direction", ha="center", fontsize=10)
ax.text(0, -1.45, "static, purely magnetic (E = 0): J = 0 (rigid spin is a gauge turn)\n"
                  "spinning: needs electric flux Q, with J/Q fixed by the winding\n(Volkov–Wöhnert 2003; van der Bij–Radu 2002)", ha="center", fontsize=9)
ax.set_title("B. Gauss's law turns the bubble's angular momentum\ninto electric flux through a far sphere (textbook)", fontsize=10)
fig.tight_layout(); fig.savefig("figures/t4-rotation-is-a-gauge-turn.png", dpi=140, bbox_inches="tight"); plt.close(fig)

# ---------------------------------------------------------------- t5
fig, axs = plt.subplots(1, 3, figsize=(15, 4.6))
d = [float(k_) for k_ in R["S5"]["J"]]; Jv = list(R["S5"]["J"].values())
axs[0].plot(d, Jv, "o", color=ORANGE, ms=9, label="computed ∫ r × (E × B) d³x")
axs[0].axhline(R["S5"]["target"], color=GREY, ls="--", label="eg/4π")
axs[0].set_ylim(0, 0.1); axs[0].set_xlabel("separation of charge and monopole"); axs[0].set_ylabel("field angular momentum")
axs[0].set_title("A. A charge beside a monopole spins the field\nby eg/4π, at any distance (Thomson; S5)", fontsize=10); axs[0].legend(fontsize=8.5)
nn = np.arange(1, 6)
axs[1].scatter(4 * np.pi * nn / 2, np.ones(5) * 1, color=ORANGE, s=80, label="rotation: J = n·ħ/2 (Saha 1936, Wilson 1949)")
axs[1].scatter(2 * np.pi * nn, np.ones(5) * 0, color=BLUE, s=80, label="phase: e^(ieg/ħ) single-valued (Dirac 1931)")
for v in 2 * np.pi * nn: axs[1].plot([v, v], [0, 1], color=GREY, lw=0.6, ls=":")
axs[1].set_yticks([0, 1]); axs[1].set_yticklabels(["from i", "from spin"]); axs[1].set_ylim(-0.6, 1.8)
axs[1].set_xlabel("allowed eg (ħ = 1)"); axs[1].legend(fontsize=8, loc="upper left")
axs[1].set_title("B. Quantized rotation and a single-valued complex phase\ngive the same charge lattice (textbook)", fontsize=10)
# Berry monopole
th_ = np.linspace(0, np.pi, 30); ph_ = np.linspace(0, 2 * np.pi, 60)
TH, PH = np.meshgrid(th_, ph_)
axs[2].remove(); ax3 = fig.add_subplot(1, 3, 3, projection="3d")
ax3.plot_surface(np.sin(TH) * np.cos(PH), np.sin(TH) * np.sin(PH), np.cos(TH), color=BLUE, alpha=0.18, lw=0)
for t0 in [0.5, 1.0, 1.5, 2.0, 2.6]:
    for p0 in np.linspace(0, 2 * np.pi, 8, endpoint=False):
        v = np.array([np.sin(t0) * np.cos(p0), np.sin(t0) * np.sin(p0), np.cos(t0)])
        ax3.quiver(*v, *(0.3 * v), color=BLUE, lw=1)
ax3.set_title(f"C. A spin following the direction over the sphere:\nChern number {R['S7']['chern_complex']:+.0f}, a unit monopole.\n"
              "Real amplitudes can't: the gap closes, phases only 0 or π (S7)", fontsize=9.5)
ax3.set_axis_off()
fig.tight_layout(); fig.savefig("figures/t5-rotation-phase-charge.png", dpi=140, bbox_inches="tight"); plt.close(fig)

# ---------------------------------------------------------------- t6
fig, axs = plt.subplots(1, 3, figsize=(15, 4.6))
Hq = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4)); Hq = (Hq + Hq.conj().T) / 2
cs = np.linspace(-20, 20, 81)
P = [np.abs(expm(-1j * (Hq + c * np.eye(4)) * 1.7))[0, 1] ** 2 for c in cs]
axs[0].plot(cs, P, color=BLUE, lw=2.2)
axs[0].set_ylim(0, max(P) * 2); axs[0].set_xlabel("constant added to every energy, c")
axs[0].set_ylabel("a transition probability")
axs[0].set_title("A. Quantum: a constant energy is a pure phase.\nNothing measurable moves (S6): why-i's assumption A4", fontsize=10)
t = np.linspace(0, 2, 100)
for rho, ls in [(0.01, ":"), (0.04, "--"), (0.16, "-")]:
    axs[1].plot(t, np.exp(np.sqrt(8 * np.pi * rho / 3) * t), color=ORANGE, ls=ls, lw=2, label=f"vacuum energy ρ = {rho}")
axs[1].set_xlabel("time"); axs[1].set_ylabel("size of the universe a(t)")
axs[1].set_title("B. Gravity: the same constant gravitates.\nIt sets the expansion (Λ = 8πGρ): A4 fails", fontsize=10); axs[1].legend(fontsize=8.5)
N = 40; tau = 2.0 / N
zr = np.cumprod(np.full(N, 1 + np.sqrt(8 * np.pi * 0.04 / 3) * tau))
zi = np.cumprod(np.full(N, 1 - 1j * 2.5 * tau)) / np.abs(1 - 1j * 2.5 * tau) ** np.arange(1, N + 1)
axs[2].plot(zi.real, zi.imag, "o-", color=BLUE, ms=3, label="k imaginary: (1 + kτ)ⁿ turns (quantum phase)")
axs[2].plot(zr, 0 * zr, "o-", color=ORANGE, ms=3, label="k real = √(Λ/3): (1 + kτ)ⁿ grows (de Sitter, round 16)")
axs[2].add_patch(plt.Circle((0, 0), 1, fc="none", ec=GREY, ls=":"))
axs[2].set_aspect("equal"); axs[2].set_xlim(-1.4, 3.6); axs[2].set_ylim(-1.6, 1.6)
axs[2].legend(fontsize=7.8, loc="lower right")
axs[2].set_title("C. One fused shift a ↦ 1 + k·a carries both.\nThe constant is invisible as a turn, visible as growth", fontsize=10)
fig.tight_layout(); fig.savefig("figures/t6-quantum-gravity-edge.png", dpi=140, bbox_inches="tight"); plt.close(fig)

# ---------------------------------------------------------------- t7
rows = [
    ("Rotation of a gauge soliton is a gauge turn; J = surface flux = charge", "Textbook", "Forgács–Manton 1980; Volkov–Wöhnert 2003; van der Bij–Radu 2002"),
    ("Spinning Bartnik–McKinnon bubble must be charged", "Textbook (perturbative)", "BHSV 1997; Volkov–Gal'tsov 1999 review"),
    ("Exact rotating regular bubble exists", "Open", "conjectured absent: van der Bij–Radu 2002"),
    ("Rotation quantum = phase single-valued (same charge lattice)", "Textbook", "Saha 1936; Wilson 1949; Dirac 1931"),
    ("Spin following a direction carries a monopole; needs complex phase", "Textbook", "Wu–Yang 1975; Berry 1984"),
    ("Shift lemma links A4 (forces i) and c₀ (lets the bubble spin)", "Ours (reading)", "not found in print in my search"),
    ("Law's f = 1 leaves exactly one quarter turn: 3 / 1 / 3 / 0", "Ours (computed)", "S3; Abelian projection is 't Hooft 1974"),
    ("Same fused shift: imaginary k = phase, real k = Λ growth", "Ours (reading)", "rounds 15–16; Wick rotation is textbook"),
    ("Energy gravitates, so A4 fails under gravity", "Textbook", "the cosmological constant problem (Weinberg 1989)"),
    ("A tie that fixes ħ, Λ or a new number", "Not found", "the tie stops here"),
]
fig, ax = plt.subplots(figsize=(13, 5.6)); ax.axis("off")
colmap = {"Textbook": "#fbe3cc", "Textbook (perturbative)": "#fbe3cc", "Open": "#f5d0d0", "Ours (reading)": "#d6e6f5",
          "Ours (computed)": "#d6e6f5", "Not found": "#eeeeee"}
tab = ax.table(cellText=[[a_, b_, c_] for a_, b_, c_ in rows], colLabels=["Claim", "Status", "Source"],
               cellColours=[["white", colmap[b_], "white"] for _, b_, _ in rows], loc="center", colWidths=[0.5, 0.16, 0.34])
tab.auto_set_font_size(False); tab.set_fontsize(9); tab.scale(1, 1.6)
ax.set_title("t7. What is textbook, what is ours, what is open", fontsize=12)
fig.savefig("figures/t7-scorecard.png", dpi=140, bbox_inches="tight"); plt.close(fig)
print("wrote figures t1–t7")
