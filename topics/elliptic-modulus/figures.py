"""Pictures for the elliptic-modulus page (DRAFT). Every curve is computed from a formula or from heights.py.
Run: python3 figures.py -> figures/fig1..fig5 (PNG). Colour key: blue = round side, orange = relativity side,
green = ours, grey = known background."""
import numpy as np, math, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.special import ellipj, ellipk
from heights import height, G
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
BLUE, ORANGE, GREEN, GREY, RED, PURPLE = "#2f6fb5", "#e07b1f", "#2a9d55", "#9a9a9a", "#c8382f", "#7a4fb5"
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("bo", [BLUE, "#7a7a7a", ORANGE])
def save(fig, name): fig.savefig(f"figures/{name}.png", dpi=150, bbox_inches="tight"); plt.close(fig)

# ---------------- fig1: the pendulum dial
fig = plt.figure(figsize=(15.5, 4.6)); gs = fig.add_gridspec(1, 3, width_ratios=[1.5, 1.3, 1.0], wspace=0.28)
ax = fig.add_subplot(gs[0]); ax.set_aspect("equal"); ax.axis("off")
for i, (th0, lab) in enumerate([(30, "30°"), (90, "90°"), (150, "150°"), (180, "180°")]):
    cx = i*2.9; mm = math.sin(math.radians(th0)/2)**2; kap = 1 - 2*mm; col = cmap(mm)
    ax.add_patch(plt.Circle((cx, 0), 1, fill=False, ls=":", color=GREY))
    a = np.radians(np.linspace(-th0, th0, 60)); ax.plot(cx + np.sin(a), -np.cos(a), color=col, lw=3, alpha=0.5)
    r = math.radians(th0); ax.plot([cx, cx + math.sin(r)], [0, -math.cos(r)], color="k", lw=1.5)
    ax.plot(cx + math.sin(r), -math.cos(r), "o", color=col, ms=10)
    ax.text(cx, -1.45, "swing %s\nm = %.3f\nκ = %+.2f" % (lab, mm, kap), ha="center", va="top", fontsize=9)
ax.set_xlim(-1.3, 10.1); ax.set_ylim(-2.6, 1.3)
ax.set_title("(a) Swing amplitude θ₀ sets Jacobi's modulus m = sin²(θ₀/2)\nand the curvature reading κ = 1 − 2m = cos θ₀", fontsize=10)
ax = fig.add_subplot(gs[1])
for mm in [0, 0.15, 0.3, 0.5, 0.7, 0.85, 0.97, 1.0]:
    K = ellipk(mm) if mm < 1 else 99; uu = np.linspace(0, min(K*0.985, 3.0), 400)
    g = np.tan(uu) if mm == 0 else (np.tanh(uu) if mm == 1 else G(uu, mm))
    ax.plot(uu, g, color=cmap(mm), lw=2 if mm in (0, 0.5, 1) else 1.2, label="m = %g" % mm if mm in (0, 0.5, 1) else None)
    if mm < 1 and K < 3: ax.axvline(K, color=cmap(mm), ls=":", lw=0.8)
ax.set_ylim(0, 4); ax.set_xlim(0, 3); ax.set_xlabel("u (the law's own coordinate: time on the pendulum)")
ax.set_ylabel("slope reading t = G_m(u) = sn dn / cn")
ax.text(1.62, 3.6, "tan (m = 0): round side", color=BLUE); ax.text(1.75, 0.75, "tanh (m = 1): relativity side", color=ORANGE)
ax.text(1.88, 2.4, "dotted: the wall u = K(m),\nwhere the pendulum turns back", fontsize=8.5, color="#555")
ax.set_title("(b) One family from tan to tanh. Each wall is a turning point;\nat the separatrix it runs off to infinity (light speed)", fontsize=10)
ax = fig.add_subplot(gs[2])
mm = np.linspace(0, 0.9999, 400); ax.plot(mm, ellipk(mm), color=GREEN, lw=2, label="K(m), distance to the wall")
ax.plot(mm[200:], np.log(4/np.sqrt(1 - mm[200:])), "--", color=GREY, label="ln(4/√(1−m))")
ax.axhline(np.pi/2, color=BLUE, ls=":", lw=1); ax.text(0.02, 1.62, "π/2: the circle's quarter turn", color=BLUE, fontsize=8.5)
ax.axvline(0.5, color=PURPLE, ls=":"); ax.text(0.52, 3.3, "m = 1/2: 90° swing,\nκ = 0, square torus", color=PURPLE, fontsize=8.5)
ax.set_xlabel("m"); ax.set_ylim(1.4, 6); ax.legend(loc="upper left", fontsize=8.5, frameon=False, bbox_to_anchor=(0, 0.85))
ax.set_title("(c) The wall moves out as the swing grows\nand reaches infinity at the separatrix", fontsize=10)
save(fig, "fig1-pendulum-dial")

# ---------------- fig2: periods tell the country
fig, axs = plt.subplots(1, 5, figsize=(15, 3.6))
for ax, mm in zip(axs, [0.0, 0.2, 0.5, 0.8, 1.0]):
    if 0 < mm < 1:
        K, Kp = ellipk(mm), ellipk(1 - mm)
        for i in range(-2, 3):
            for j in range(-2, 3):
                ax.plot(2*i*K, 2*j*Kp, "o", color=GREEN, ms=5)                 # zeros
                ax.plot(2*i*K + K, 2*j*Kp, "x", color=RED, ms=6, mew=2)        # poles (walls)
                ax.plot(2*i*K, 2*j*Kp + Kp, "x", color=RED, ms=6, mew=2)
        ax.add_patch(plt.Rectangle((0, 0), 2*K, 2*Kp, fill=False, ec="k", lw=1.2))
        ax.set_title("m = %g: torus, 2K × 2K'\n(K = %.2f, K' = %.2f)%s" % (mm, K, Kp, "  square" if mm == 0.5 else ""), fontsize=9.5)
    elif mm == 0:
        for i in range(-2, 3):
            ax.axvline(i*np.pi, color=GREEN, lw=1); ax.axvline(i*np.pi + np.pi/2, color=RED, lw=1, ls="--")
        ax.set_title("m = 0: tan, one REAL period π\nturning country (circle)", fontsize=9.5, color=BLUE)
    else:
        for j in range(-2, 3):
            ax.axhline(j*np.pi, color=GREEN, lw=1); ax.axhline(j*np.pi + np.pi/2, color=RED, lw=1, ls="--")
        ax.set_title("m = 1: tanh, one IMAGINARY period iπ\ntwo-wall country (relativity)", fontsize=9.5, color=ORANGE)
    ax.set_xlim(-5, 5); ax.set_ylim(-5, 5); ax.set_aspect("equal"); ax.set_xlabel("Re u"); ax.set_ylabel("Im u")
axs[0].text(-4.8, -4.7, "green: zeros   red: poles (walls)", fontsize=8)
fig.suptitle("Fig 2. Where G_m = sn dn / cn is zero (green) and infinite (red) in the complex u-plane. Our conic laws have one period; "
             "its direction is the country.\nBetween them every m gives two periods, a torus. At m = 1/2 the torus is square and Wick rotation u → iu is a symmetry.", fontsize=10, y=1.08)
save(fig, "fig2-periods-and-countries")

# ---------------- fig3: the fingerprint
fig, axs = plt.subplots(1, 2, figsize=(12.5, 4.6))
ax = axs[0]
for sgn in (-1, 1):
    k = sgn*np.linspace(0.18, 1.7, 300); ax.plot(k, 1.2 - 0.9/k**2, color=GREEN, lw=2)
ax.axhline(1.2, color=GREY, ls="--"); ax.text(-1.68, 1.25, "6/5: the sine chart (only class the slope chart misses)", fontsize=8.5, color="#555")
ax.plot([-1, 1], [0.3, 0.3], "o", color="k", ms=8, zorder=5)
ax.annotate("tanh, κ = −1\n(relativity) 3/10", (-1, 0.3), xytext=(-1.65, -0.9), color=ORANGE, arrowprops=dict(arrowstyle="->", color=ORANGE))
ax.annotate("tan, κ = +1\n(circle) 3/10", (1, 0.3), xytext=(1.1, -0.9), color=BLUE, arrowprops=dict(arrowstyle="->", color=BLUE))
ax.axvline(0, color=PURPLE, ls=":"); ax.text(0.04, -2.6, "κ = 0: 90° swing,\nthe square torus;\ncontact jumps 3 → 5", color=PURPLE, fontsize=8.5)
ax.axhline(0, color=GREY, lw=0.6); ax.plot([-np.sqrt(3)/2, np.sqrt(3)/2], [0, 0], "s", color=RED, ms=6)
ax.text(0.25, -0.45, "R = 0 at κ = ±√3/2:\nthe 30° and 150° swings", color=RED, fontsize=8.5)
ax.set_ylim(-3, 1.6); ax.set_xlim(-1.7, 1.7); ax.set_xlabel("curvature reading κ = 1 − 2m = cos θ₀"); ax.set_ylabel("fingerprint R (iterative residue)")
ax.set_title("(a) R = 6/5 − 9/(10κ²). Our conic laws all sit at one value, 3/10;\nthe modulus moves the fingerprint, symmetric under Wick rotation κ → −κ", fontsize=10)
ax = axs[1]
kk, ee = np.meshgrid(np.linspace(-1.6, 1.6, 401), np.linspace(-1.2, 2.6, 401))
with np.errstate(divide="ignore", invalid="ignore"): R = 1.2 - 0.9*ee/kk**2
im = ax.contourf(kk, ee, np.clip(R, -3, 1.6), levels=np.linspace(-3, 1.6, 24), cmap="PuOr_r", alpha=0.8); plt.colorbar(im, ax=ax, label="R")
k = np.linspace(-1.6, 1.6, 200)
ax.plot(k, k**2, color="k", lw=2.2, label="ε = κ²: our conic laws (perfect squares), R = 3/10")
ax.plot(k, 0*k, color=GREY, lw=2, ls="--", label="ε = 0: sine charts (arcsin), R = 6/5")
for c in (0.4, 2.5): ax.plot(k, c*k**2, color="k", lw=0.7, ls=":")
ax.plot(0, 0, "o", color="k"); ax.text(0.05, -0.15, "plain addition", fontsize=8)
ax.plot(1 - 2*np.array([0.1, 0.3, 0.7, 0.9]), [1, 1, 1, 1], "o", color=GREEN, ms=6, label="Jacobi slope charts (ε = 1, κ = 1 − 2m)")
ax.set_ylim(-1.2, 2.6); ax.set_xlim(-1.6, 1.6); ax.set_xlabel("κ (t² coefficient / 2)"); ax.set_ylabel("ε (t⁴ coefficient)")
ax.legend(fontsize=8, loc="upper center", frameon=True)
ax.set_title("(b) Odd laws: log' = 1/√(1 + 2κt² + εt⁴). Rescaling slides along ε = cκ² (dotted)\nand keeps R; this is the k = 0 twin of the (k, d) residue map", fontsize=10)
save(fig, "fig3-fingerprint")

# ---------------- fig4: heights
fig, axs = plt.subplots(1, 3, figsize=(15, 4.5))
x = np.linspace(0.0, 1.0, 41)
for ax, mm in zip(axs[:2], (0.2, 0.8)):
    for sv, ls in [(-1, ":"), (-0.5, "--"), (0, "-"), (0.5, "--"), (1, "-"), (2, ":")]:
        xx = x if sv <= 1 else x[x < (0.85 if mm == 0.2 else 1.0)]
        y = xx if sv == 0 else height(xx, mm, sv)
        ax.plot(xx, y - xx, ls, color=cmap(mm) if sv != 0 else GREY, lw=2 if sv in (1,) else 1.3, label="height %g" % sv)
    ax.set_xlabel("x"); ax.set_ylabel("G^s(x) − x"); ax.legend(fontsize=8, frameon=False, ncol=2)
    ax.set_title("(%s) Heights of the m = %g law (%s side)\nhalf steps compose to one step to 10⁻¹¹" % ("a" if mm == 0.2 else "b", mm, "round" if mm < .5 else "relativity"), fontsize=10)
ax = axs[2]
def mprime(R0, s):
    with np.errstate(invalid="ignore", divide="ignore"):
        k2 = 9/(12 - 10*R0/s); return (1 - np.sqrt(k2))/2
ss = np.logspace(-2, 2, 400)
for m0, col in [(0.0, BLUE), (0.2, cmap(0.25)), (0.4, PURPLE)]:
    R0 = 1.2 - 0.9/(1 - 2*m0)**2; mp_ = mprime(R0, ss)
    ax.semilogx(ss, mp_, color=col, lw=2, label="start m = %g or %g (R = %.2f)" % (m0, 1 - m0, R0)); ax.plot(1, m0, "o", color=col)
ax.axhline(math.sin(math.radians(15))**2, color=RED, ls="--"); ax.text(0.012, 0.1, "m = sin²15°: the 30° swing (R = 0)", color=RED, fontsize=8.5)
ax.axhline(0.5, color=PURPLE, ls=":"); ax.text(0.011, 0.47, "m = 1/2", color=PURPLE, fontsize=8.5)
ax.set_ylim(-0.4, 0.55); ax.set_xlabel("height s (dot: s = 1, the law itself)"); ax.set_ylabel("modulus m' that height s matches (formally)")
ax.legend(fontsize=8, frameon=False, loc="lower left")
ax.set_title("(c) Heights move the modulus: G_m^s has R/s, which is the\nfingerprint of the elliptic law at m'. Every climb heads to the 30° swing", fontsize=10)
save(fig, "fig4-heights-move-the-modulus")

# ---------------- fig5: theta = the lift copied round the torus
fig, axs = plt.subplots(1, 3, figsize=(15, 4.4))
ax = axs[0]
z = np.linspace(-14, 1.5, 3000); q = 0.35
lift = 1 + z
n = np.arange(-30, 31)
z = np.linspace(-14, -0.35, 3000); lift = 1 + z
theta = np.array([np.sum(q**(n*(n - 1)/2)*zz**n) if zz != 0 else 1 for zz in z])
ax.plot(z, np.arcsinh(lift), color=GREY, lw=2, label="the lift 1 + z: one wall at z = −1")
ax.plot(z, np.arcsinh(theta), color=GREEN, lw=1.8, label="θ(z) = Σ q^(n(n−1)/2) zⁿ, q = 0.35")
for j in range(0, 3): ax.axvline(-q**(-j), color=RED, ls=":", lw=1)
ax.text(-13.5, -3.6, "red: the wall −1 and its copies −1/q, −1/q², ...", color=RED, fontsize=8.5)
ax.axhline(0, color="k", lw=0.5); ax.set_ylim(-4, 5); ax.set_ylabel("arcsinh of value (to fit both)"); ax.set_xlabel("z"); ax.legend(fontsize=8, frameon=False, loc="upper left")
ax.set_title("(a) Theta = product of lifts (1 + z qʲ):\nmult-add's one wall copied round the torus", fontsize=10)
ax = axs[1]; ax.axis("off")
ax.set_title("(b) The corners of an infinite cube, one edge per step of the orbit z, zq, zq², ...", fontsize=10)
rows = [("{}", "1"), ("{0}", "z"), ("{1}", "z·q"), ("{0,1}", "z²·q"), ("{2}", "z·q²"), ("{0,2}", "z²·q²"), ("{1,2}", "z²·q³"), ("{0,1,2}", "z³·q³")]
for i, (S, w) in enumerate(rows):
    ax.text(0.05, 0.9 - i*0.105, "corner %s" % S, fontsize=10); ax.text(0.45, 0.9 - i*0.105, "weight %s" % w, fontsize=10, color=GREEN)
ax.text(0.0, 0.02, "corner weight = zⁿ q^(sum of positions): no pair terms.\nGenus one is pure mult-add; pair weights need two or more waves.", fontsize=9, color="#444")
ax = axs[2]
xx = np.linspace(-3, 3, 600)
for qv, col in [(0.05, BLUE), (1e-3, PURPLE), (1e-6, GREEN), (1e-14, ORANGE)]:
    nn = np.arange(-30, 31)[:, None]; w = qv**(nn*(nn - 1)/2)*np.exp(2*nn*xx[None, :])
    t0, t1, t2 = w.sum(0), (2*nn*w).sum(0), (4*nn*nn*w).sum(0)
    ax.plot(xx, 2*(t2/t0 - (t1/t0)**2), color=col, lw=1.6, label="q = %g" % qv)
ax.plot(xx, 2/np.cosh(xx)**2, "k:", lw=2, label="2 sech²x (one soliton)")
ax.set_xlabel("x"); ax.set_ylabel("u = 2 ∂²ₓ log τ"); ax.legend(fontsize=8, frameon=False)
ax.set_title("(c) The same τ as a KdV wave: cnoidal for q > 0,\nthe soliton when only corners n = 0, 1 survive", fontsize=10)
save(fig, "fig5-theta-lift-round-the-torus")
print("figures written")
