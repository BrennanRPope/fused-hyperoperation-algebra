"""Pictures for the named-functions page. Run: python3 -I figures.py"""
import sys, os
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np, mpmath as mp
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm
from flows import E_k, E_kd
mp.mp.dps = 25
OUT = os.path.join(HERE, "figures")
BLUE, ORANGE, GREY = "#1f5fbf", "#e07b00", "#888888"
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})
E = E_k(1)
def vec(fl, s, xs, n=150):
    out = []
    for x in xs:
        try: out.append(float(fl(s, x, n=n)))
        except Exception: out.append(np.nan)
    return np.array(out)
HEIGHTS = [-1, 0, 0.5, 1, 1.5, 2]
cmap = cm.viridis
col = lambda s: cmap((s + 1)/3.2)

# ---------- fig1: the hypernatural pair ----------
fig, ax = plt.subplots(1, 2, figsize=(12, 5.4))
X = np.linspace(0.15, 4.0, 160)
for s in HEIGHTS:
    y = vec(E, -s, X - 1)
    ax[0].plot(X, y, color=col(s), lw=2.2 if s == 1 else 1.6, label=f"height s = {s:g}" + ("  (ln)" if s == 1 else "  (X − 1)" if s == 0 else ""))
ax[0].axhline(0, color=GREY, lw=.6); ax[0].axvline(1, color=GREY, lw=.6)
ax[0].plot([1], [0], "ko"); ax[0].set_ylim(-2.5, 3.2)
ax[0].set_title("Hypernatural logarithm ln⁽ˢ⁾(X)\nevery height passes through (1, 0) with slope 1; its bend there is −s")
ax[0].set_xlabel("X"); ax[0].legend(fontsize=9, loc="lower right")
T = np.linspace(-2.2, 1.6, 160)
for s in HEIGHTS:
    y = 1 + vec(E, s, T)
    ax[1].plot(T, y, color=col(s), lw=2.2 if s == 1 else 1.6, label=f"s = {s:g}" + ("  (eᵗ)" if s == 1 else "  (1 + t)" if s == 0 else ""))
for s in [0, 0.5, 1, 1.5, 2]:
    es = 1 + float(E(s, 1, n=150)); ax[1].plot([1], [es], "o", color=col(s))
    ax[1].annotate(f"e⁽{s:g}⁾ = {es:.3f}", (1, es), xytext=(6, -4), textcoords="offset points", fontsize=8, color=col(s))
ax[1].axhline(1, color=GREY, lw=.6); ax[1].axvline(0, color=GREY, lw=.6); ax[1].set_ylim(-0.2, 7)
ax[1].set_title("Its inverse, the hypernatural exponential exp⁽ˢ⁾(t) = 1 + Eˢ(t)\nvalue at t = 1 is e at that height")
ax[1].set_xlabel("t"); ax[1].legend(fontsize=9, loc="upper left")
fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig1-hypernatural-pair.png"), dpi=130); plt.close(fig)

# ---------- fig2: shapes ----------
panels = [("Cube  (k = 1, d = 0)", E_kd(1, 0), "yes: the original"),
          ("Simplex shares  (k = −1)", E_kd(-1, 0), "yes: same flow, other half"),
          ("Ellipse  (k = 1, d = 1)", E_kd(1, 1), "yes: its own flow"),
          ("Parabola  (k = 1, d = 1/4)", E_kd(1, mp.mpf(1)/4), "yes: exact, height = slide"),
          ("Circle  (k = 0, d = 1)", E_kd(0, 1), "yes: odd flow (iterated tan)"),
          ("Relativity  (k = 0, d = −1)", E_kd(0, -1), "yes: odd flow (iterated tanh)"),
          ("Plain addition  (k = 0, d = 0)", None, "no: every height is the same"),
          ("ℓp balls (p dial)", "lp", "no new function: height = p dial")]
fig, axs = plt.subplots(2, 4, figsize=(16, 8.2))
xs = np.linspace(-0.9, 0.9, 91)
for ax, (title, fl, verdict) in zip(axs.flat, panels):
    if fl is None:
        for s in HEIGHTS: ax.plot(xs, xs, color=col(s), lw=1.6)
        ax.set_ylim(-1.2, 1.2)
    elif fl == "lp":
        th = np.linspace(0, 2*np.pi, 400)
        for s, c in zip([0, 0.5, 1, 1.5, 2, 3.3], [col(-1), col(0), col(0.5), col(1), col(1.5), col(2)]):
            p = 2.0**s
            x = np.sign(np.cos(th))*np.abs(np.cos(th))**(2/p); y = np.sign(np.sin(th))*np.abs(np.sin(th))**(2/p)
            ax.plot(x, y, color=c, lw=1.6, label=f"rung {s:g}: p = {p:.3g}")
        ax.set_aspect("equal"); ax.legend(fontsize=7, loc="center")
    else:
        for s in HEIGHTS:
            y = vec(fl, -s, xs, n=120)
            y[np.abs(y) > 3] = np.nan
            ax.plot(xs, y, color=col(s), lw=2.0 if s == 1 else 1.4)
        ax.axhline(0, color=GREY, lw=.5); ax.axvline(0, color=GREY, lw=.5); ax.set_ylim(-2.2, 2.2)
    good = verdict.startswith("yes")
    ax.set_title(title, fontsize=11)
    ax.text(0.03, 0.95 if fl != "lp" else 0.02, verdict, transform=ax.transAxes, fontsize=9.5, va="top" if fl != "lp" else "bottom",
            color=BLUE if good else ORANGE, fontweight="bold")
fig.suptitle("Does each shape have a hypernatural logarithm?  Curves: ln⁽ˢ⁾ for that shape's law at heights s = −1 … 2 (yellow = 2, purple = −1).\nHeight 1 (thick) is the shape's own logarithm; height 0 is the straight line.", fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.92]); fig.savefig(os.path.join(OUT, "fig2-shapes.png"), dpi=120); plt.close(fig)

# ---------- fig3: the (k, d) map ----------
fig, ax = plt.subplots(figsize=(8.5, 6.5))
K = np.linspace(-2.2, 2.2, 400); D = np.linspace(-1.5, 2, 400)
KK, DD = np.meshgrid(K, D); Delta = KK**2 - 4*DD
ax.contourf(KK, DD, np.sign(Delta), levels=[-1.5, 0, 1.5], colors=["#dbe7fb", "#fde9d2"], alpha=0.9)
ax.plot(K, K**2/4, color="magenta", lw=2.2, label="parabola class Δ = 0: exact flow, height is a slide")
ax.axvline(0, color=BLUE, lw=2, ls="--", label="k = 0: odd flow (iterated tan / tanh)")
ax.plot([0], [0], "o", ms=12, mfc="white", mec="k", mew=2, label="k = d = 0, plain addition: no flow")
for (k, d, nm) in [(1, 0, "cube"), (-1, 0, "simplex shares"), (0, 1, "circle"), (0, -1, "relativity"), (1, 1, "ellipse")]:
    ax.plot([k], [d], "ko"); ax.annotate(nm, (k, d), xytext=(6, 5), textcoords="offset points")
ax.text(-2.05, 1.7, "Δ < 0: round side\nflow like eˣ − 1, no closed form", fontsize=9)
ax.text(0.9, -1.3, "Δ > 0: hyperbola side\nflow like eˣ − 1, no closed form", fontsize=9)
ax.set_xlabel("homothety dial k"); ax.set_ylabel("curvature dial d")
ax.set_title("Every conic on the plane dial has a hypernatural logarithm,\nexcept the single point of plain addition")
ax.legend(fontsize=8.5, loc="upper right"); fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig3-kd-map.png"), dpi=130); plt.close(fig)

# ---------- fig4: two halves ----------
fig, ax = plt.subplots(1, 2, figsize=(12, 4.8))
for s in [0.5, 1.5]:
    Xr = np.linspace(1.0005, 3, 80); Xl = np.linspace(0.12, 0.9995, 80)
    ax[0].plot(Xr, vec(E, -s, Xr - 1), color=GREY, lw=2)
    ax[0].plot(Xl, vec(E, -s, Xl - 1), color=ORANGE, lw=2)
    ax[0].annotate(f"s = {s}", (2.7, float(E(-s, 1.7, n=100))), fontsize=9)
ax[0].axvline(1, color="k", lw=.6); ax[0].set_xlabel("X")
ax[0].set_title("One function, two halves\ngrey: growth side X > 1 (the cube's corner sums)\norange: share side X < 1 (stick-breaking, k = −1)", fontsize=10)
h = np.linspace(-0.25, 0.25, 101); h = h[h != 0]
for s in [0.5, 1.5]:
    y = vec(E, -s, h)
    ax[1].plot(h, (y - h)/h**2, color=col(s), lw=2, label=f"s = {s}: bend/2 → −s/2 = {-s/2}")
ax[1].axvline(0, color="k", lw=.6); ax[1].set_xlabel("X − 1")
ax[1].set_title("(ln⁽ˢ⁾ − (X − 1)) / (X − 1)²: both halves agree at 1\n(they join smoothly; Baker 1958: not analytically, at fractional s)", fontsize=10)
ax[1].legend(fontsize=9); fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig4-two-halves.png"), dpi=130); plt.close(fig)

# ---------- fig5: Bell numbers at height s ----------
fig, ax = plt.subplots(1, 2, figsize=(12, 4.8))
S = np.linspace(-1, 3.2, 120)
for n, c in zip(range(2, 7), [col(-1), col(0), col(0.7), col(1.4), col(2)]):
    vals = [float(E.coeffs(mp.mpf(s))[n]*mp.factorial(n)) for s in S]
    ax[0].plot(S, np.sign(vals)*np.log10(1 + np.abs(vals)), color=c, lw=2, label=f"n = {n}")
    for si in range(-1, 4):
        v = float(E.coeffs(mp.mpf(si))[n]*mp.factorial(n))
        ax[0].plot([si], [np.sign(v)*np.log10(1 + abs(v))], "o", color=c)
for si, lab in [(1, "sets\n(all 1s)"), (2, "Bell\nnumbers"), (3, "partitions of\npartitions")]:
    ax[0].annotate(lab, (si, 4.1), ha="center", fontsize=8.5)
ax[0].set_ylim(-3, 4.8); ax[0].set_xlabel("height s"); ax[0].set_ylabel("sign · log₁₀(1 + |bₙ(s)|)")
ax[0].set_title("bₙ(s) = n-th derivative of exp⁽ˢ⁾ at 0\npolynomials in s through the counting numbers", fontsize=10)
ax[0].legend(fontsize=8, loc="lower right")
# nested groupings picture for n = 3
ax[1].axis("off")
ax[1].set_title("What the height counts (n = 3 items)", fontsize=10)
rows = [("s = 1: one way (a set)", 1), ("s = 2: 5 ways to group  (Bell)", 5), ("s = 3: 12 ways to group, then group the groups", 12)]
for i, (lab, cnt) in enumerate(rows):
    ax[1].text(0.02, 0.85 - i*0.3, lab, fontsize=10, transform=ax[1].transAxes)
    for j in range(cnt):
        ax[1].add_patch(plt.Circle((0.06 + j*0.075, 0.74 - i*0.3), 0.025, transform=ax[1].transAxes, color=col(i)))
ax[1].text(0.02, 0.05, "Fractional heights interpolate these counts with Comtet's triangle (OEIS A008826).", fontsize=9, transform=ax[1].transAxes)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig5-bell-at-height.png"), dpi=130); plt.close(fig)

# ---------- fig6: names map ----------
fig, ax = plt.subplots(figsize=(13, 6.4)); ax.axis("off")
def box(x, y, txt, known, w=0.28):
    ax.text(x, y, txt, ha="left", va="center", fontsize=10.5, transform=ax.transAxes,
            bbox=dict(boxstyle="round,pad=0.4", fc="#fff4e6" if known else "#e8f0fc", ec=ORANGE if known else BLUE))
ax.text(0.02, 0.95, "Height goes up top (superscript)", fontsize=12, fontweight="bold", transform=ax.transAxes)
ax.text(0.52, 0.95, "Dials go below (subscript)", fontsize=12, fontweight="bold", transform=ax.transAxes)
H = [("Eˢ  height flow of eˣ − 1 (regular iteration)", True), ("Φ, Ψ  Abel function and superfunction of E", True),
     ("itlog E  speed of the height flow (Julia, Écalle)", True), ("a +⁽ˢ⁾ b  addition at height s", False),
     ("X ×⁽ˢ⁾ Y  multiplication at height s", False), ("ln⁽ˢ⁾  hypernatural logarithm", False),
     ("exp⁽ˢ⁾  hypernatural exponential", False), ("e⁽ˢ⁾ = exp⁽ˢ⁾(1)  e at height s", False),
     ("bₙ(s)  Bell numbers at height s (Bell 1938, Comtet)", True)]
Dl = [("a ∘ₖ b  circle operation", True), ("ℓₖ lift, φₖ flip, σ_c slide", False), ("Lₖ  logarithm of ∘ₖ", True),
      ("Eₖ = (e^{kt} − 1)/k  exponential of ∘ₖ (E₁ = E)", True), ("L_{k,d}, E_{k,d}  log and exp of each conic", True),
      ("C_κ, S_κ, T_κ  Cayley–Klein trig, κ = −Δ/4 a reading", True), ("[n]_q, q = 1 + k  q-integer", True),
      ("Uₙ(k, d)  Lucas sequence (powers of ε)", True), ("∘₋₁  probabilistic sum (noisy-OR)", True)]
for i, (t, kn) in enumerate(H): box(0.02, 0.86 - i*0.095, t, kn)
for i, (t, kn) in enumerate(Dl): box(0.52, 0.86 - i*0.095, t, kn)
ax.text(0.02, 0.0, "orange = existing name in the literature    blue = our name or symbol (marked ours)", fontsize=10, transform=ax.transAxes)
fig.savefig(os.path.join(OUT, "fig6-names-map.png"), dpi=130, bbox_inches="tight"); plt.close(fig)

# ---------- fig7: Cayley-Klein exponential of the plane number ----------
fig, axs = plt.subplots(1, 3, figsize=(14, 4.6))
t = np.linspace(0, 9, 600)
for ax, (k, d, nm) in zip(axs, [(-0.3, 1, "round side (Δ < 0): spiral"), (-1, 0.25, "parabola (Δ = 0): no turn, only a slide"), (-0.3, -0.5, "hyperbola side (Δ > 0): runs off")]):
    kap = d - k*k/4
    if kap > 0: C, S = np.cos(np.sqrt(kap)*t), np.sin(np.sqrt(kap)*t)/np.sqrt(kap)
    elif kap < 0: C, S = np.cosh(np.sqrt(-kap)*t), np.sinh(np.sqrt(-kap)*t)/np.sqrt(-kap)
    else: C, S = np.ones_like(t), t
    g = np.exp(k*t/2)
    ax.plot(C*g, S*g, color=BLUE, lw=2); ax.plot(C, S, color=GREY, lw=1, ls="--")
    ax.set_title(f"{nm}\nk = {k}, d = {d}, curvature −Δ/4 = {kap:.3g}", fontsize=10)
    ax.set_xlabel("C_κ part"); ax.set_ylabel("S_κ part"); ax.set_aspect("equal", adjustable="datalim")
    ax.set_xlim(-3, 3); ax.set_ylim(-3, 3)
fig.suptitle("exp(t·ε) = e^{kt/2} (C_κ(t) + S_κ(t)·ε₀): blue: the true path (k < 0 shrinks it); grey: the curvature part alone (turn, slide or stretch)", fontsize=10.5)
fig.tight_layout(rect=[0, 0, 1, 0.9]); fig.savefig(os.path.join(OUT, "fig7-cayley-klein-exp.png"), dpi=130); plt.close(fig)
print("figures written")
