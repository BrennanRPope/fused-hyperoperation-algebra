"""Pictures for the Lee-Yang checks (Nov 27): fig9 the ring's zeros as a pendulum, fig10 the 2D dial gas on the pendulum dial.
Run from this folder: python3 ly_figures.py -> ../../figures/fig9, fig10."""
import sys, math, cmath, pickle, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.special import ellipk
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
BLUE, ORANGE, GREEN, GREY, RED, PURPLE = "#2f6fb5", "#e07b1f", "#2a9d55", "#9a9a9a", "#c8382f", "#7a4fb5"
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("bo", [BLUE, "#7a7a7a", ORANGE])
W = dict(fc="white", ec="none", alpha=0.85, pad=1.5)
def save(fig, name): fig.savefig(f"../../figures/{name}.png", dpi=150, bbox_inches="tight"); plt.close(fig)
def ring_zeta(k, n):
    out = []
    for j in range(n):
        c = 4*math.cos(math.pi*(2*j+1)/(2*n))**2; s = c*(1 - 1/k) - 2; d = cmath.sqrt(s*s - 4)
        out += [(s + d)/2, (s - d)/2]
    return np.array(out)

# ---------------- fig9
fig, axs = plt.subplots(1, 3, figsize=(16, 5.2), gridspec_kw=dict(width_ratios=[1.05, 1.1, 1.25], wspace=0.28))
ax = axs[0]; ax.set_aspect("equal"); ax.axis("off")
k, n = 3.0, 41; m = 1 - 1/k; th0 = 2*math.asin(math.sqrt(m))
ax.add_patch(plt.Circle((0, 0), 1, fill=False, ls=":", color=GREY))
# pendulum drawn hanging from the centre: bottom = zeta = -1 (left), top = zeta = +1 (right); rotate so the wall is at the bottom
rot = lambda zeta: zeta * (-1j)     # zeta = -1 -> +i? we want zeta=-1 at the bottom: multiply by +1j*(-1) ...
rot = lambda zeta: zeta * 1j       # zeta = -1 -> -i (bottom), zeta = +1 -> +i (top)
a = np.linspace(-th0, th0, 200); arc = rot(-np.exp(-1j*a))
ax.plot(arc.real, arc.imag, color=cmap(m), lw=7, alpha=0.35, solid_capstyle="round")
zs = rot(ring_zeta(k, n)); ax.plot(zs.real, zs.imag, "o", color=GREEN, ms=4.5)
end = rot(-np.exp(-1j*th0)); ax.plot([0, end.real], [0, end.imag], color="k", lw=1.5); ax.plot(end.real, end.imag, "o", color=cmap(m), ms=11)
ax.plot(0, -1, "s", color=RED, ms=7); ax.text(0, -1.17, "ζ = −1: the wall\n= the pendulum at rest", ha="center", va="top", color=RED, fontsize=9)
ax.plot(0, 1, "^", color=PURPLE, ms=8); ax.text(0, 1.14, "ζ = +1: zero field\n= the top", ha="center", color=PURPLE, fontsize=9)
ax.text(end.real - 0.12, end.imag + 0.02, "Yang–Lee edge =\nturning point θ₀,\nsin²(θ₀/2) = 1 − 1/k", fontsize=9, color="#555", ha="right")
ax.set_xlim(-2.3, 1.5); ax.set_ylim(-1.7, 1.5)
ax.set_title("(a) A ring's Lee–Yang zeros (green, k = 3, 41 sites) on the\ncircle |ζ| = 1, drawn as a pendulum hanging from the wall", fontsize=10)
ax = axs[1]
for k in [1.5, 3.0, 30.0]:
    m = 1 - 1/k; th0 = 2*math.asin(math.sqrt(m))
    nz = 40001; ang = np.array([2*math.asin(math.sqrt(m)*abs(math.cos(math.pi*(2*j+1)/(2*nz)))) for j in range(nz)])
    h, e = np.histogram(ang, bins=60, range=(0, th0), density=True); c = (e[1:] + e[:-1])/2
    ax.plot(c, h, "o", color=cmap(m), ms=3.5)
    t = np.linspace(0, th0*0.9995, 400); dwell = 1/(2*np.sqrt(m - np.sin(t/2)**2))
    pred = np.cos(t/2)*dwell; pred /= np.trapezoid(pred, t)
    ax.plot(t, pred, color=cmap(m), lw=1.6, label=f"k = {k:g}  (m = {m:.2f}, θ₀ = {math.degrees(th0):.0f}°)")
    dn = dwell/np.trapezoid(dwell, t); ax.plot(t, dn, color=cmap(m), lw=1, ls=":")
ax.set_ylim(0, 2.9); ax.set_xlabel("angle from the wall, θ = π − arg ζ"); ax.set_ylabel("zeros per unit angle")
from matplotlib.lines import Line2D
h0, l0 = ax.get_legend_handles_labels()
h0 += [Line2D([], [], marker="o", ls="", color="k", ms=3.5), Line2D([], [], color="k", lw=1.6), Line2D([], [], color="k", lw=1, ls=":")]
l0 += ["zeros of a 40,001-site ring", "cos(θ/2) × the pendulum's dwell time", "dwell time alone"]
ax.legend(h0, l0, fontsize=8.5, frameon=False, loc="upper left")
ax.set_title("(b) The zeros crowd toward the edge exactly as a\npendulum lingers near its turning point (× cos(θ/2))", fontsize=10)
ax = axs[2]; ax.set_aspect("equal")
ax.add_patch(plt.Circle((0, 0), 1, fill=False, ls=":", color=GREY))
for k, colr, lab in [(3.0, cmap(2/3), "k = 3: swinging (m = 2/3)"), (-1.0, PURPLE, "k = −1: over the top (m = 2)"),
                     (0.5, "#555", "k = 1/2: Wick image (m = −1)")]:
    zs = ring_zeta(k, 61); ax.plot(zs.real, zs.imag, "o", color=colr, ms=4, label=lab)
ax.set_xlim(-3.2, 3.2); ax.set_ylim(-1.6, 1.6); ax.axhline(0, color=GREY, lw=0.6); ax.axvline(0, color=GREY, lw=0.6)
ax.set_xlabel("Re ζ"); ax.set_ylabel("Im ζ"); ax.legend(fontsize=8.5, loc="lower left", frameon=True)
ax.text(1.1, 0.25, "spills past the top\nonto the positive axis", color=PURPLE, fontsize=8.5)
ax.text(-3.1, 0.25, "leaves the circle\nonto the negative axis", color="#555", fontsize=8.5)
ax.set_title("(c) The ring's three regimes are the pendulum's three\nenergy bands, with m = 1 − 1/k", fontsize=10)
save(fig, "fig9-ring-zeros-as-a-pendulum")

# ---------------- fig10
fig, axs = plt.subplots(1, 2, figsize=(14, 5.4), gridspec_kw=dict(width_ratios=[1, 1.15], wspace=0.25))
ax = axs[0]; ax.set_aspect("equal")
phi = np.linspace(0, 2*np.pi, 800); A = np.cos(phi) + 2; r1 = A + np.sqrt(A*A - 1); r2 = 1/r1
for r in (r1, r2): ax.plot(r*np.cos(phi), r*np.sin(phi), color=ORANGE, lw=2)
Fz = pickle.load(open('../../../lee-yang/fisher_zeros.pkl', 'rb'))
for Lf, mk in [(4, "s"), (5, "^"), (6, "o")]:
    z = np.array(Fz[Lf]); ax.plot(z.real, z.imag, mk, color=GREEN, ms=5, mfc="none" if Lf < 6 else GREEN, label=f"{Lf}×{Lf} torus")
kc = 3 + 2*math.sqrt(2); ax.plot([kc, 1/kc], [0, 0], "*", color=RED, ms=13)
ax.text(kc, -0.55, "k_c = 3 + 2√2", color=RED, ha="center", fontsize=9)
ax.set_xlim(-2.5, 7.2); ax.set_ylim(-4.3, 4.3); ax.axhline(0, color=GREY, lw=0.6)
ax.legend(fontsize=8.5, loc="upper left", frameon=False); ax.set_xlabel("Re k"); ax.set_ylabel("Im k")
ax.text(1.0, 3.05, "orange: |k_O| = 1, the swap line\nbetween swinging and going over", color=ORANGE, fontsize=9, bbox=W)
ax.set_title("(a) The 2D dial gas's zero-field (Fisher) zeros in the dial k\nclose in on the curve where Onsager's modulus has size 1", fontsize=10)
ax = axs[1]
kk = np.linspace(1.02, 40, 4000); kO = 4*kk/(kk - 1)**2
ax.plot(kk, kO, color="k", lw=1.5, label="Onsager's modulus k_O = 4k/(k − 1)²")
ax.fill_between(kk, 0, 3, where=kk > kc, color=BLUE, alpha=0.08); ax.fill_between(kk, 0, 3, where=kk < kc, color=PURPLE, alpha=0.08)
ax.axhline(1, color=ORANGE, lw=1.5, ls=":"); ax.axvline(kc, color=RED, lw=1.2)
ax.set_xscale("log"); ax.set_ylim(0, 3); ax.set_xlabel("dial k (log scale)"); ax.set_ylabel("k_O")
ax2 = ax.twinx(); ax2.spines["right"].set_visible(True)
m = np.minimum(kO, 1/kO)**2; ax2.plot(kk, ellipk(m), color=GREEN, lw=2, label="K: time to the top")
ax2.set_ylabel("K(m), the pendulum's time to the top", color=GREEN); ax2.set_ylim(1.5, 6)
ax.text(12, 2.6, "swinging\n(ordered)", color=BLUE, fontsize=9.5); ax.text(1.4, 2.6, "over the top\n(disordered)", color=PURPLE, fontsize=9.5)
ax.text(kc*1.07, 0.15, "k_c: the separatrix", color=RED, fontsize=9.5)
ax.legend(loc="center right", fontsize=8.5, frameon=False); ax2.legend(loc="upper right", bbox_to_anchor=(1, 0.8), fontsize=8.5, frameon=False)
ax.set_title("(b) Along the real dial: the transition is where the modulus hits 1.\nThe log peak of K there is Onsager's log in the specific heat (known)", fontsize=10)
save(fig, "fig10-dial-gas-on-the-pendulum-dial")
print("done")
