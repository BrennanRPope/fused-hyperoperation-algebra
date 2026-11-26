"""Labeled pictures for the velocity-lines DRAFT. Writes figures/*.png."""
import math, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle
from scipy.integrate import solve_ivp

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "figures"); os.makedirs(OUT, exist_ok=True)
BLUE, ORANGE, RED, GREEN, GRAY, PURPLE = "#2a6fdb", "#e8833a", "#d23c3c", "#2e9e5b", "#777777", "#7b4fc9"
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})

def save(fig, name):
    fig.savefig(os.path.join(OUT, name), dpi=150, bbox_inches="tight"); plt.close(fig); print("wrote", name)

def orbit(p, x0, v0, T, n=4000):
    def f(t, y):
        r = math.hypot(y[0], y[1]); a = -1/r**p/r
        return [y[2], y[3], a*y[0], a*y[1]]
    s = solve_ivp(f, (0, T), list(x0) + list(v0), method="DOP853", rtol=1e-12, atol=1e-13, t_eval=np.linspace(0, T, n))
    return s.t, s.y[0] + 1j*s.y[1], s.y[2] + 1j*s.y[3]

def dial_add(u, v, k):
    return (u + v)/(1 - k*np.conj(u)*v)

def to_sphere(V, k):
    w = V*math.sqrt(k); w2 = np.abs(w)**2
    return np.stack([2*w.real, 2*w.imag, w2 - 1], 1)/(w2 + 1)[:, None]

# ---------------------------------------------------------------- fig 1: seen from the start velocity, it runs straight
t, X, V = orbit(2, (1, 0), (0.1, 1.15), 11.6, n=3000)
E = 0.5*abs(V[0])**2 - 1/abs(X[0]); k = -1/(2*E)
fig, axs = plt.subplots(1, 3, figsize=(15, 5))
ax = axs[0]; ax.set_aspect("equal")
ax.plot(X.real, X.imag, color=BLUE, lw=2); ax.plot(0, 0, "o", color=RED); ax.plot(X[0].real, X[0].imag, "o", color="k")
ax.text(X[0].real + 0.05, X[0].imag + 0.05, "start"); ax.set_title("1. An orbit under 1/r² (positions)")
ax = axs[1]; ax.set_aspect("equal")
ax.plot(V.real, V.imag, color=BLUE, lw=2, label="velocity v(t): a circle")
ax.plot(V[0].real, V[0].imag, "o", color="k"); ax.text(V[0].real + 0.04, V[0].imag, "v₀")
ax.plot(0, 0, "+", color="k"); ax.text(0.03, -0.1, "v = 0")
ax.legend(frameon=False, loc="lower left"); ax.set_title("2. Its velocities (hodograph)")
ax = axs[2]; ax.set_aspect("equal")
for kk, col, lab in ((k, BLUE, f"dial κ = −1/(2E) = {k:.2f}: one straight line"), (1.7*k, ORANGE, "dial ×1.7: bends"),
                     (0.0, GRAY, "flat dial κ = 0 (plain subtraction): bends")):
    Z = dial_add(-V[0], V, kk); Z = Z[np.abs(Z) < 6]
    ax.plot(Z.real, Z.imag, color=col, lw=2.5 if col == BLUE else 1.5, label=lab)
ax.plot(0, 0, "o", color="k"); ax.text(0.1, 0.25, "v₀ moved to 0")
ax.set_xlim(-3, 3); ax.set_ylim(-3, 3)
ax.legend(frameon=False, fontsize=9, loc="lower left")
ax.set_title("3. Every velocity measured from v₀ with the\ntheory's dial addition: (−v₀) ⊕ v(t)")
save(fig, "fig1_runs_straight.png")

# ---------------------------------------------------------------- fig 2: the dial's clock
t, X, V = orbit(2, (1, 0), (0.1, 1.15), 11.6, n=3000)
E = 0.5*abs(V[0])**2 - 1/abs(X[0]); k = -1/(2*E); a = -1/(2*E)
Z = dial_add(-V[0], V, k); dist = 2*np.arctan(math.sqrt(k)*np.abs(Z))
half = np.argmax(dist > 3.1)
fig, axs = plt.subplots(1, 2, figsize=(13, 4.8))
ax = axs[0]
tau = np.concatenate([[0], np.cumsum(np.diff(t)/np.abs(X[1:]))])*math.sqrt(-2*E)
ax.plot(t[:half], dist[:half], color=BLUE, lw=2.5, label="dial distance travelled by the velocity")
ax.plot(t[:half], tau[:half], "--", color=RED, lw=2, label="the clock ∫ √(−2E) dt / r")
ax.set_xlabel("ordinary time t"); ax.set_ylabel("radians on the dial sphere"); ax.legend(frameon=False)
ax.set_title("Steady pace, if time is counted as dt/r\n(this clock is Kepler's eccentric anomaly)")
ax = axs[1]
lam = 2/(1 + k*np.abs(V)**2)
ax.plot(t, np.abs(X), color=BLUE, lw=2.5, label="distance from the sun r")
ax.plot(t[::60], a*lam[::60], "o", color=ORANGE, ms=5, label="a × λ(v), the dial chart's scale at v")
ax.set_xlabel("time t"); ax.legend(frameon=False)
ax.set_title("Position is read off the velocity geometry:\nr = a × 2/(1 + κ|v|²)")
save(fig, "fig2_dial_clock.png")

# ---------------------------------------------------------------- fig 3: one dial per energy
E0 = -0.3; k0 = -1/(2*E0)
fig = plt.figure(figsize=(13, 5.4))
ax = fig.add_subplot(1, 2, 1); ax.set_aspect("equal")
ax3 = fig.add_subplot(1, 2, 2, projection="3d")
uu, vv = np.mgrid[0:2*np.pi:40j, 0:np.pi:20j]
ax3.plot_wireframe(np.cos(uu)*np.sin(vv), np.sin(uu)*np.sin(vv), np.cos(vv), color=GRAY, lw=0.3, alpha=0.5)
cols = [BLUE, GREEN, ORANGE, PURPLE]
for (h_, rot), col in zip(((0.35, 0.0), (0.7, 1.2), (1.0, 2.4), (1.15, 4.0)), cols):
    v0m = math.sqrt(2*(E0 + 1)); vr = math.sqrt(v0m**2 - h_**2)
    t, X, V = orbit(2, (1, 0), (vr, h_), 30, n=40000)
    V = V*np.exp(1j*rot)
    ax.plot(V.real, V.imag, color=col, lw=2, label=f"angular momentum {h_}")
    S = to_sphere(V, k0); ax3.plot(S[:, 0], S[:, 1], S[:, 2], color=col, lw=2.5)
ax.plot(0, 0, "+", color="k"); ax.legend(frameon=False, fontsize=9, loc="lower left")
ax.set_title(f"Four orbits, same energy E = {E0}:\nfour different circles of velocity")
ax3.set_box_aspect((1, 1, 1)); ax3.set_xticks([]); ax3.set_yticks([]); ax3.set_zticks([])
ax3.set_title(f"On the one dial sphere κ = −1/(2E) = {k0:.2f}\nall four are great circles (straight lines)")
save(fig, "fig3_one_dial_per_energy.png")

# ---------------------------------------------------------------- fig 4: the dial family, Kepler and relativity side by side
fig, axs = plt.subplots(1, 3, figsize=(15, 5), subplot_kw=dict(projection="3d"))
th = np.linspace(0, 2*np.pi, 60)
titles = ["bound orbits, E < 0\nvelocity geometry: sphere (d = +1)\nsymmetry: rotations in 4D, SO(4)",
          "parabolic, E = 0\nflat (d = 0)\nsymmetry: Galilean-like E(3)",
          "escaping, E > 0\nhyperboloid (d = −1)\nsymmetry: Lorentz group SO(3,1)"]
for ax, ttl, kind in zip(axs, titles, ("s", "f", "h")):
    if kind == "s":
        ax.plot_wireframe(np.cos(uu)*np.sin(vv), np.sin(uu)*np.sin(vv), np.cos(vv), color=BLUE, lw=0.3, alpha=0.6)
    elif kind == "f":
        xx, yy = np.meshgrid(np.linspace(-1.5, 1.5, 12), np.linspace(-1.5, 1.5, 12))
        ax.plot_wireframe(xx, yy, 0*xx, color=GRAY, lw=0.4)
    else:
        rr_, tt_ = np.meshgrid(np.linspace(0, 1.4, 12), th)
        ax.plot_wireframe(np.sinh(rr_)*np.cos(tt_), np.sinh(rr_)*np.sin(tt_), np.cosh(rr_), color=ORANGE, lw=0.4)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([]); ax.set_title(ttl, fontsize=10)
fig.suptitle("Orbits' energy sets the dial (Fock, Bargmann, Moser, Milnor). Same family as relativity's boosts: rotations · Galilean · Lorentz.\n"
             "Bracket of the Runge–Lenz vectors: {A_x, A_y} = −2E·L_z, so −2E plays the dial", fontsize=11)
save(fig, "fig4_dial_family.png")

# ---------------------------------------------------------------- fig 5: dimensions: flux law vs straight velocity lines
fig, axs = plt.subplots(2, 3, figsize=(13, 8))
for j, (n, p, col) in enumerate(((2, 1.0, PURPLE), (3, 2.0, BLUE), (4, 3.0, ORANGE))):
    T = 30 if p < 3 else 2.4
    t, X, V = orbit(p, (1, 0), (0, 0.9 if p < 3 else 1.02), T, n=6000)
    ax = axs[0, j]; ax.set_aspect("equal"); ax.plot(X.real, X.imag, color=col, lw=1); ax.plot(0, 0, "o", color=RED)
    ax.set_title(f"{n} dimensions: flux law gives pull 1/r^{n-1}"); ax.set_xticks([]); ax.set_yticks([])
    ax = axs[1, j]; ax.set_aspect("equal"); ax.plot(V.real, V.imag, color=col, lw=1)
    ax.set_title("velocities: " + ("a circle ✓" if n == 3 else "not a circle ✗")); ax.set_xticks([]); ax.set_yticks([])
fig.suptitle("Two rules for gravity: 'conserved flux' depends on the dimension, 'velocities run straight' always wants 1/r².\nThey agree only in 3 dimensions")
save(fig, "fig5_dimensions.png")

# ---------------------------------------------------------------- fig 6: hydrogen: the same sphere, quantised
fig, axs = plt.subplots(1, 2, figsize=(13, 5))
ax = axs[0]
for n in range(1, 5):
    En = -1/(2*n*n)
    ax.hlines(En, 0, 1, color=BLUE, lw=3)
    xs = np.linspace(1.15, 1.15 + 0.08*(n*n - 1), n*n)
    ax.plot(xs, [En]*len(xs), "o", color=ORANGE, ms=5)
    ax.text(-0.05, En, f"n = {n}", ha="right", va="center")
    ax.text(1.2 + 0.08*n*n, En, f"{n*n} states, dial κ = {n*n}", va="center", fontsize=10)
ax.set_xlim(-0.4, 4.2); ax.set_ylim(-0.56, 0.02); ax.set_xticks([])
ax.set_ylabel("energy (atomic units)"); ax.set_title("Hydrogen: level n has n² states and dial κ = −1/(2E) = n²")
ax = axs[1]; ax.set_aspect("equal"); ax.axis("off")
for n, col in zip(range(1, 5), (BLUE, GREEN, ORANGE, PURPLE)):
    ax.add_patch(Circle((0, 0), 1/n, fill=False, lw=2.5, ec=col)); ax.text(1/n*0.72, 1/n*0.72, f"n={n}", color=col)
ax.set_xlim(-1.2, 1.2); ax.set_ylim(-1.2, 1.2)
ax.set_title("The dial sphere in momentum space shrinks to radius 1/n\n(Fock 1935: hydrogen = waves on this sphere)")
save(fig, "fig6_hydrogen.png")

# ---------------------------------------------------------------- fig 7: map
fig, ax = plt.subplots(figsize=(12, 6.4)); ax.axis("off"); ax.set_xlim(0, 12); ax.set_ylim(0, 7.4)
def box(x, y, w, h, txt, fc, ec):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.1", fc=fc, ec=ec, lw=1.8))
    ax.text(x + w/2, y + h/2, txt, ha="center", va="center", fontsize=10)
box(0.2, 5.6, 5.4, 1.5, "POSTULATE (new, stated in the theory's words)\nfree fall = velocity runs straight in one dial\ngeometry per energy, at a steady pace", "#fde8e8", RED)
box(6.4, 5.6, 5.4, 1.5, "FROM THE THEORY\ndial addition u ⊕ v = (u+v)/(1 − κ ū v)\nand its sphere / flat / hyperboloid family", "#e6f4ea", GREEN)
box(0.2, 3.1, 3.6, 1.6, "pull ∝ 1/r²\n(no dimension, no flux assumed)", "#fff6e0", ORANGE)
box(4.2, 3.1, 3.6, 1.6, "dial = −1/(2E)\nclock = dt/r, r = a·λ(v)\nKepler's equation", "#fff6e0", ORANGE)
box(8.2, 3.1, 3.6, 1.6, "symmetry L, A close\nwith coefficient −2E\nSO(4) / E(3) / SO(3,1)", "#fff6e0", ORANGE)
box(0.2, 0.5, 5.4, 1.6, "+ Gauss (conserved flux in n dimensions)\n⇒ n = 3", "#eef3fb", BLUE)
box(6.4, 0.5, 5.4, 1.6, "quantise the same sphere ⇒ hydrogen:\nn² states, κ = n² (Fock 1935)", "#eef3fb", BLUE)
for x0 in (2.9, 9.1):
    for x1 in (2.0, 6.0, 10.0):
        ax.annotate("", xy=(x1, 4.8), xytext=(x0, 5.5), arrowprops=dict(arrowstyle="->", lw=1.2, color=GRAY))
ax.annotate("", xy=(2.9, 2.2), xytext=(2.0, 3.0), arrowprops=dict(arrowstyle="->", lw=1.4, color=GRAY))
ax.annotate("", xy=(9.1, 2.2), xytext=(10.0, 3.0), arrowprops=dict(arrowstyle="->", lw=1.4, color=GRAY))
ax.set_title("What is assumed and what follows")
save(fig, "fig7_map.png")
