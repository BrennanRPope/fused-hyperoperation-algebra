"""Pictures for the thermodynamics page. Run: python3 figures.py  -> figures/*.png"""
import itertools, math, random
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle, Polygon

HOT, COLD, INK, MUTE, PALE_H, PALE_C = "#c8321f", "#1f63a8", "#16181d", "#6b7280", "#f7c8bd", "#c3dcf3"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11, "axes.edgecolor": INK,
                     "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK})
OUT = "figures/"
random.seed(3); np.random.seed(3)


def save(fig, name):
    fig.savefig(OUT + name, dpi=160, bbox_inches="tight", facecolor="white")
    plt.close(fig)


# ---------- F1: the square three ways ----------
fig, axs = plt.subplots(1, 3, figsize=(13.5, 4.6))
# (a) partition functions
ax = axs[0]; zA, zB = 0.8, 0.55
ax.add_patch(Rectangle((0, 0), 1, 1, fc="white", ec=INK, lw=1.5))
ax.add_patch(Rectangle((1, 0), zA, 1, fc=PALE_H, ec=INK, lw=1.5))
ax.add_patch(Rectangle((0, 1), 1, zB, fc=PALE_C, ec=INK, lw=1.5))
ax.add_patch(Rectangle((1, 1), zA, zB, fc=HOT, ec=INK, lw=1.5, alpha=.85))
ax.text(.5, .5, "both in\nground state\n1", ha="center", va="center")
ax.text(1 + zA / 2, .5, "A excited\nz_A", ha="center", va="center")
ax.text(.5, 1 + zB / 2, "B excited\nz_B", ha="center", va="center")
ax.text(1 + zA / 2, 1 + zB / 2, "both\nexcited\nz_A·z_B", ha="center", va="center", color="white", weight="bold")
ax.set_title("Two systems: Z_AB = (1 + z_A)(1 + z_B)\nthe corner of mult-add is 'both excited'", fontsize=11)
ax.set_xlim(-.05, 1.9); ax.set_ylim(-.05, 1.62); ax.set_aspect("equal"); ax.axis("off")
# (b) engines in series: the waste box shrinks
ax = axs[1]; e1, e2 = .35, .3
ax.add_patch(Rectangle((0, 0), 1, 1, fc=HOT, ec=INK, lw=1.5, alpha=.85))
ax.add_patch(Rectangle((e1, 0), 1 - e1, 1, fc=PALE_H, ec=INK, lw=1.5))
ax.add_patch(Rectangle((e1, e2), 1 - e1, 1 - e2, fc=PALE_C, ec=INK, lw=1.5))
ax.text(e1 / 2, .5, "work\nstage 1\nη₁", ha="center", va="center", color="white", weight="bold")
ax.text(e1 + (1 - e1) / 2, e2 / 2, "work stage 2\nη₂(1 − η₁)", ha="center", va="center", fontsize=10)
ax.text(e1 + (1 - e1) / 2, e2 + (1 - e2) / 2, "waste heat\n(1 − η₁)(1 − η₂)", ha="center", va="center")
ax.set_title("Two engines in series: η = η₁ ⊞₋₁ η₂\nthe waste box multiplies (lift 1 − η)", fontsize=11)
ax.set_xlim(-.05, 1.05); ax.set_ylim(-.05, 1.05); ax.set_aspect("equal"); ax.axis("off")
# (c) fridges in cascade: the rejected-heat box grows
ax = axs[2]; u1, u2 = .45, .35
ax.add_patch(Rectangle((0, 0), 1, 1, fc=PALE_C, ec=INK, lw=1.5))
ax.add_patch(Rectangle((1, 0), u1, 1, fc=PALE_H, ec=INK, lw=1.5))
ax.add_patch(Rectangle((0, 1), 1 + u1, u2 * 1, fc=HOT, ec=INK, lw=1.5, alpha=.85))
ax.plot([1, 1], [1, 1 + u2], color=INK, lw=1)
ax.text(.5, .5, "cold heat\npumped\nQ_c = 1", ha="center", va="center")
ax.text(1 + u1 / 2, .5, "work 1\nu₁", ha="center", va="center")
ax.text((1 + u1) / 2, 1 + u2 / 2, "work 2: u₂(1 + u₁)", ha="center", va="center", color="white", weight="bold")
ax.set_title("Two fridges in cascade: u = u₁ ⊞ u₂ (k = +1)\nrunning an engine backwards is the inverse", fontsize=11)
ax.set_xlim(-.05, 1.6); ax.set_ylim(-.05, 1.5); ax.set_aspect("equal"); ax.axis("off")
save(fig, "f1_three_squares.png")

# ---------- F2: fermion cube and boson staircase ----------
fig = plt.figure(figsize=(12.5, 5))
ax = fig.add_subplot(1, 2, 1)
P = lambda bits: np.array([bits[0] * 1.0 + bits[2] * .45, bits[1] * 1.0 + bits[2] * .38])
for bits in itertools.product((0, 1), repeat=3):
    for d in range(3):
        if bits[d] == 0:
            nb = list(bits); nb[d] = 1
            p, q = P(bits), P(nb)
            ax.plot([p[0], q[0]], [p[1], q[1]], color=MUTE, lw=1.2, zorder=1)
labels = {1: "z", 2: "z²", 3: "z³"}
for bits in itertools.product((0, 1), repeat=3):
    p = P(bits); nocc = sum(bits)
    col = ["white", PALE_H, HOT, "#7a1a0e"][nocc]
    ax.add_patch(Circle(p, .12, fc=col, ec=INK, lw=1.3, zorder=3))
    word = "".join("●" if bb else "○" for bb in bits)
    ax.text(p[0], p[1] - .22, word, ha="center", va="top", fontsize=10)
ax.text(.73, -.62, "8 corners = 8 ways to fill 3 fermion modes (Pauli: each 0 or 1)\n"
        "Ξ = (1 + z₁)(1 + z₂)(1 + z₃) = sum over corners", ha="center", fontsize=10.5)
ax.set_title("Fermions are the mult-add cube", fontsize=12)
ax.set_xlim(-.4, 1.85); ax.set_ylim(-.85, 1.65); ax.set_aspect("equal"); ax.axis("off")
ax = fig.add_subplot(1, 2, 2)
zb = .62; x0 = 0
for n in range(9):
    w = zb ** n
    ax.add_patch(Rectangle((x0, 0), w, w, fc=PALE_C if n else "white", ec=INK, lw=1.2))
    if n < 4:
        ax.text(x0 + w / 2, w / 2, "1" if n == 0 else f"z^{n}", ha="center", va="center", fontsize=10)
    x0 += w
ax.annotate("", xy=(1 / (1 - zb) + .02, .12), xytext=(x0 - .2, .12),
            arrowprops=dict(arrowstyle="->", color=MUTE))
ax.text(1.3, -.18, "a boson mode holds 0, 1, 2, … quanta: Ξ_B = 1 + z + z² + … = 1/(1 − z)\n"
        "Ξ_B(z) · Ξ_F(−z) = 1: the boson is the way back (inverse) from a fermion of weight −z",
        ha="center", va="top", fontsize=10.5)
ax.set_title("Bosons are the inverse of a negative fermion", fontsize=12)
ax.set_xlim(-.05, 2.75); ax.set_ylim(-.62, 1.1); ax.set_aspect("equal"); ax.axis("off")
save(fig, "f2_fermion_cube_boson_stairs.png")

# ---------- F3: Fermi function between two walls, and negative temperature ----------
fig, axs = plt.subplots(1, 2, figsize=(12.5, 4.4))
ax = axs[0]
xs = np.linspace(-6, 6, 400)
ax.axhspan(-1.25, -1, color=PALE_C); ax.axhspan(1, 1.25, color=PALE_H)
ax.plot(xs, np.tanh(xs / 2), color=INK, lw=2)
ax.axhline(1, color=HOT, lw=2); ax.axhline(-1, color=COLD, lw=2); ax.axhline(0, color=MUTE, lw=.8, ls=":")
ax.text(-5.8, 1.07, "wall: n = 0 (empty), T → 0⁺", color=HOT, fontsize=10)
ax.text(-5.8, -1.18, "wall: n = 1 (full), T → 0⁻", color=COLD, fontsize=10)
ax.text(1.2, .1, "identity: n = ½, T = ∞", color=MUTE, fontsize=10)
ax.set_xlabel("(ε − μ)/kT"); ax.set_ylabel("t = 1 − 2n = tanh(x/2)")
ax.set_title("Fermi function: energies add, t adds like velocities\n(t₁ + t₂)/(1 + t₁t₂), walls ±1", fontsize=11)
ax.set_ylim(-1.25, 1.25)
ax = axs[1]
th = np.linspace(-1, 1, 9)
ax.add_patch(Rectangle((-1, -.06), 2, .12, fc="#eef0f3", ec="none"))
ax.plot([-1, 1], [0, 0], color=INK, lw=1.5)
for tv, lab in [(-1, "T = 0⁻"), (-.6, "T < 0"), (0, "T = ±∞"), (.6, "T > 0"), (1, "T = 0⁺")]:
    col = HOT if tv > 0 else (COLD if tv < 0 else MUTE)
    ax.add_patch(Circle((tv, 0), .045 if abs(tv) < 1 else .06, fc=col, ec=INK, zorder=3))
    ax.text(tv, .16, lab, ha="center", fontsize=10, color=col)
ax.text(0, -.32, "hotter than infinity is not past a wall: it is the other side of the identity.\n"
        "The only walls are absolute zero from either side (Ramsey 1956, read on the segment).",
        ha="center", va="top", fontsize=10)
ax.set_xlim(-1.25, 1.25); ax.set_ylim(-.7, .4); ax.set_aspect("equal"); ax.axis("off")
ax.set_title("Two-level system on the velocity segment", fontsize=11)
save(fig, "f3_fermi_walls.png")

# ---------- F4: Ising series, parallel and the lens ----------
fig, axs = plt.subplots(1, 2, figsize=(12.5, 4.8))
ax = axs[0]
def spin(xy, col="white"):
    ax.add_patch(Circle(xy, .09, fc=col, ec=INK, lw=1.4, zorder=3))
spin((0, 1)); spin((1, 1), "#ddd"); spin((2, 1))
ax.plot([0, 2], [1, 1], color=HOT, lw=3, zorder=1)
ax.text(.5, 1.12, "K₁", ha="center"); ax.text(1.5, 1.12, "K₂", ha="center")
ax.text(1, .72, "series: sum out the middle spin\ntanh K = tanh K₁ · tanh K₂  (parity, k = −2)", ha="center", fontsize=10)
spin((0, 0)); spin((2, 0))
arc = np.linspace(0, math.pi, 50)
ax.plot(1 + np.cos(arc), .0 + .28 * np.sin(arc), color=COLD, lw=3, zorder=1)
ax.plot(1 + np.cos(arc), .0 - .28 * np.sin(arc), color=COLD, lw=3, zorder=1)
ax.text(1, .34, "K₁", ha="center"); ax.text(1, -.43, "K₂", ha="center")
ax.text(1, -.62, "parallel: couplings add\ntanh K = velocity sum of tanh K₁, tanh K₂", ha="center", va="top", fontsize=10)
ax.set_xlim(-.3, 2.3); ax.set_ylim(-1.05, 1.35); ax.set_aspect("equal"); ax.axis("off")
ax.set_title("One Ising bond, two ways to combine", fontsize=11)
ax = axs[1]
t = np.linspace(0, 1, 300)
ax.plot(t, (1 - t) / (1 + t), color=INK, lw=2, label="lens M(t) = (1 − t)/(1 + t)  (Kramers–Wannier)")
ax.plot(t, t, color=MUTE, ls="--", lw=1, label="t* = t")
tc = math.sqrt(2) - 1
ax.add_patch(Circle((tc, tc), .018, fc=HOT, ec=INK, zorder=4))
ax.annotate("self-dual point t = √2 − 1\n= tanh K_c, the 2D Ising critical point", xy=(tc, tc), xytext=(.48, .62),
            arrowprops=dict(arrowstyle="->", color=INK), fontsize=10)
ax.text(.03, .05, "hot (t → 0) ↔ cold (t* → 1)", color=MUTE, fontsize=10)
ax.set_xlabel("t = tanh K"); ax.set_ylabel("dual t* = e^(−2K)")
ax.set_xlim(0, 1); ax.set_ylim(0, 1.02); ax.set_aspect("equal")
ax.legend(loc="upper right", fontsize=9, frameon=False)
ax.set_title("Duality swaps series and parallel; its fixed point is T_c\n(the same lens gives the KdV soliton weight)", fontsize=11)
save(fig, "f4_ising_lens.png")

# ---------- F5: Lee-Yang zeros ----------
def Z_poly(n, J):
    coeffs = [0.0] * (n + 1)
    for S in itertools.product((0, 1), repeat=n):
        w = 1.0
        for (i, j), Jij in J.items():
            if S[i] != S[j]:
                w *= math.exp(-2 * Jij)
        coeffs[sum(S)] += w
    return coeffs

fig, axs = plt.subplots(1, 2, figsize=(12.5, 5.2))
ax = axs[0]
circ = np.linspace(0, 2 * math.pi, 400)
ax.plot(np.cos(circ), np.sin(circ), color=MUTE, lw=1)
for trial in range(25):
    n = random.randint(4, 8)
    J = {(i, j): random.uniform(.05, 1.2) for i in range(n) for j in range(i + 1, n) if random.random() < .5}
    J.update({(i, i + 1): random.uniform(.05, 1.2) for i in range(n - 1)})
    r = np.roots(Z_poly(n, J)[::-1]); ax.plot(r.real, r.imag, "o", ms=4, color=HOT, alpha=.8)
for trial in range(8):
    n = random.randint(4, 7)
    J = {(i, i + 1): -random.uniform(.3, 1.2) for i in range(n - 1)}
    r = np.roots(Z_poly(n, J)[::-1]); ax.plot(r.real, r.imag, "s", ms=4, color=COLD, alpha=.8)
ax.plot([], [], "o", color=HOT, label="ferromagnets (25 random)"); ax.plot([], [], "s", color=COLD, label="antiferromagnets")
ax.legend(loc="lower left", fontsize=9, frameon=False)
ax.set_aspect("equal"); ax.set_xlim(-2.6, 2.6); ax.set_ylim(-2.1, 2.1)
ax.set_title("Zeros of Z in the field variable z\nferromagnets: every zero on the unit circle (Lee–Yang 1952)", fontsize=11)
ax.set_xlabel("Re z"); ax.set_ylabel("Im z")
ax = axs[1]
for a_, col in [(1.0, INK), (.6, PALE_H), (.3, HOT), (0.0, "#7a1a0e")]:
    sq = [(0, 0, "1"), (1, 0, f"{a_:g}·zᵢ"), (0, 1, f"{a_:g}·zⱼ"), (1, 1, "zᵢzⱼ")]
ax.axis("off")
x0 = 0
for a_, title in [(1.0, "T = ∞: a = 1\nplain mult-add square"), (.45, "ferromagnet: 0 < a < 1\nedges shrink"), (0.0, "T = 0: a = 0\nonly the diagonal")]:
    s0 = 1.0
    cells = [((0, 0), "white", "1"), ((1, 0), PALE_H if a_ > 0 else "#f3f3f3", f"{a_:g} zᵢ" if a_ else "0"),
             ((0, 1), PALE_H if a_ > 0 else "#f3f3f3", f"{a_:g} zⱼ" if a_ else "0"), ((1, 1), HOT, "zᵢzⱼ")]
    for (cx, cy), fc, txt in cells:
        ax.add_patch(Rectangle((x0 + cx * .62, cy * .62), .6, .6, fc=fc, ec=INK, lw=1.2, alpha=.9 if fc == HOT else 1))
        ax.text(x0 + cx * .62 + .3, cy * .62 + .3, txt, ha="center", va="center", fontsize=10,
                color="white" if fc == HOT else INK)
    ax.text(x0 + .61, -.2, title, ha="center", va="top", fontsize=10)
    x0 += 1.6
ax.text(2.2, 1.55, "each bond's Lee–Yang polynomial 1 + a zᵢ + a zⱼ + zᵢzⱼ\n"
        "is a mult-add square with its two edge cells weighted by a = e^(−2βJ)", ha="center", fontsize=10.5)
ax.text(2.2, -.75, "Asano (1970): glue the squares spin by spin, keeping only the cube's all-0 and all-1 corners;\n"
        "the zeros stay on the circle. Checked against brute force on 15 random graphs.", ha="center", va="top", fontsize=10)
ax.set_xlim(-.2, 4.6); ax.set_ylim(-1.15, 1.85); ax.set_aspect("equal")
save(fig, "f5_lee_yang.png")

# ---------- F6: Stokes pile / radiation shields on the border ----------
fig, ax = plt.subplots(figsize=(8.5, 4.4))
n = np.arange(0, 21); R = .3
T_border = 1 - n * R / (1 + (n - 1) * R); T_border[0] = 1
T_absorb = (1 - .3) ** n
ax.plot(n, T_border, "o-", color=HOT, ms=4, label="non-absorbing plates, R = 0.3 (Stokes): on the border, falls like 1/n")
ax.plot(n, T_absorb, "s-", color=COLD, ms=4, label="plates that absorb 30%: two walls, falls exponentially")
ax.set_xlabel("number of plates (or radiation shields)"); ax.set_ylabel("fraction transmitted")
ax.legend(fontsize=9, frameon=False, loc="upper right"); ax.set_ylim(0, 1.1)
ax.set_title("Same first step, different countries", fontsize=11)
save(fig, "f6_stokes_border.png")
print("figures written")
