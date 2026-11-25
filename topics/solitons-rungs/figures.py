"""Figures for the solitons-on-the-rungs page. Writes figures/*.svg."""
import sys, math, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, "/mnt/project-files/math-theory/fractional-rungs")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from regular_iteration import Es

OUT = "/mnt/project-files/math-theory/solitons-rungs/figures/"
INK = "#1d2a33"; SEA = "#0f6e7a"; SAND = "#d98c2b"; ROSE = "#b5404f"; GREY = "#9aa7ad"; FOAM = "#e6f1f2"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11, "axes.edgecolor": INK,
                     "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK, "text.color": INK,
                     "axes.spines.top": False, "axes.spines.right": False, "svg.fonttype": "none"})
def save(fig, name):
    fig.savefig(OUT + name, format="svg", bbox_inches="tight", facecolor="white")
    plt.close(fig)
def Aij(a, b): return ((a - b)/(a + b))**2
T = lambda n: n*(n - 1)//2

# 1. the theta lattice and its q = 0 face -------------------------------------------------
fig, axs = plt.subplots(1, 2, figsize=(9.2, 4.4))
A = Aij(1.0, 1.6)
for ax, q, title in [(axs[0], 0.35, "q = 0.35: the whole lattice carries weight"),
                     (axs[1], 0.0, "q = 0: only the unit square survives")]:
    ax.add_patch(plt.Rectangle((-0.18, -0.18), 1.36, 1.36, fc=FOAM, ec=SEA, lw=1.5, zorder=0))
    for n1 in range(-2, 4):
        for n2 in range(-2, 4):
            wgt = (q**T(n1) if T(n1) else 1.0) * (q**T(n2) if T(n2) else 1.0) * A**(n1*n2)
            if q == 0 and (T(n1) or T(n2)): wgt = 0
            r = 520*min(wgt, 1.0)**0.5
            on_cube = n1 in (0, 1) and n2 in (0, 1)
            ax.scatter([n1], [n2], s=max(r, 0.0), color=SEA if on_cube else SAND, zorder=2, alpha=0.9)
            ax.scatter([n1], [n2], s=6, color=INK, zorder=3)
    for (n1, n2, lab) in [(0, 0, "1"), (1, 0, "w₁"), (0, 1, "w₂"), (1, 1, "A·w₁w₂")]:
        ax.annotate(lab, (n1, n2), xytext=(8, -16), textcoords="offset points", color=INK, fontsize=10)
    ax.set_xlim(-2.6, 3.6); ax.set_ylim(-2.6, 3.6); ax.set_aspect("equal")
    ax.set_xticks(range(-2, 4)); ax.set_yticks(range(-2, 4))
    ax.set_xlabel("copies of wave 1  (n₁)"); ax.set_ylabel("copies of wave 2  (n₂)")
    ax.set_title(title, fontsize=11, loc="left")
fig.text(0.5, -0.04, "Dot area = weight q^C(n₁,2) · q^C(n₂,2) · A^(n₁n₂).  The shaded square is the circle-operation cube.",
         ha="center", fontsize=10, color=GREY)
save(fig, "lattice.svg")

# 2. q glues copies of one soliton into a cnoidal wave ------------------------------------
fig, axs = plt.subplots(4, 1, figsize=(9.2, 6.6), sharex=True)
eta = np.linspace(-12, 12, 2400)
for ax, q in zip(axs, [1e-7, 1e-2, 0.12, 0.45]):
    L = -math.log(q)
    tot = np.zeros_like(eta)
    for m in range(-80, 81):
        cp = 0.5/np.cosh(np.clip((eta - m*L)/2, -350, 350))**2
        tot += cp
        if abs(m*L) < 14: ax.plot(eta, cp, color=SAND, lw=0.9, alpha=0.8)
    ax.plot(eta, tot, color=SEA, lw=2.2)
    ax.set_ylim(0, max(tot.max()*1.15, 0.6)); ax.set_yticks([])
    ax.text(-11.8, ax.get_ylim()[1]*0.78, f"q = {q:g}" + ("   (one soliton)" if q < 1e-6 else f"   spacing −log q = {L:.2f}"),
            fontsize=10)
axs[-1].set_xlabel("position  η")
fig.text(0.5, -0.02, "Thin: copies of the q = 0 soliton.  Thick: their exact sum, the periodic (cnoidal) KdV wave.",
         ha="center", fontsize=10, color=GREY)
save(fig, "copies.svg")

# 3. what gluing looks like: the annulus and its seam --------------------------------------
fig, axs = plt.subplots(1, 2, figsize=(9.2, 4.6))
th = np.linspace(0, 2*np.pi, 400)
for ax, q, title in [(axs[0], 0.55*np.exp(0.7j), "q ≠ 0: the annulus between |q| and 1, rims glued"),
                     (axs[1], 0.04*np.exp(0.7j), "q → 0: the hole shrinks, the plane minus 0 is left")]:
    rq = abs(q)
    ax.fill(np.cos(th), np.sin(th), color=FOAM, zorder=0)
    ax.fill(rq*np.cos(th), rq*np.sin(th), color="white", zorder=1)
    ax.plot(np.cos(th), np.sin(th), color=SEA, lw=2); ax.plot(rq*np.cos(th), rq*np.sin(th), color=SEA, lw=2, ls="--")
    w = 0.78*np.exp(0.25j)
    pts = [w*q**m for m in range(-2, 7)]
    xs = [p.real for p in pts]; ys = [p.imag for p in pts]
    ax.plot(xs, ys, color=SAND, lw=1, alpha=0.7)
    ax.scatter(xs, ys, color=SAND, s=24, zorder=3)
    ax.scatter([w.real], [w.imag], color=ROSE, s=60, zorder=4)
    ax.annotate("w", (w.real, w.imag), xytext=(6, 6), textcoords="offset points", color=ROSE)
    ax.annotate("w·q", (pts[3].real, pts[3].imag), xytext=(6, -12), textcoords="offset points", color=INK, fontsize=9)
    ax.set_xlim(-1.6, 1.6); ax.set_ylim(-1.6, 1.6); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(title, fontsize=10.5, loc="left")
fig.text(0.5, 0.0, "Orange: the orbit w·qᵐ.  Multiplying by q is the circle operation by q−1, so gluing is a rung-1 move.",
         ha="center", fontsize=10, color=GREY)
save(fig, "annulus.svg")

# 4. the rung flow applied to a soliton -----------------------------------------------------
etas = np.linspace(-14, 22, 1801); h = etas[1] - etas[0]
fig, axs = plt.subplots(1, 2, figsize=(9.2, 4.0), gridspec_kw={"width_ratios": [1.25, 1]})
cols = {0.6: SAND, 0.85: "#c9a35a", 1.0: SEA, 1.3: "#6c8f95", 2.0: GREY}
for tv, col in cols.items():
    F = Es(-tv, np.exp(etas)); p = 2*np.gradient(np.gradient(F, h), h)
    m = (etas > -10) & (etas < 16)
    axs[0].plot(etas[m], p[m], color=col, lw=2.6 if tv == 1.0 else 1.4, label=f"t = {tv:g}")
axs[0].axhline(0, color=GREY, lw=0.6); axs[0].set_ylim(-0.12, 0.7)
axs[0].legend(frameon=False, fontsize=9, loc="upper right")
axs[0].set_xlabel("η"); axs[0].set_title("2·d²/dη² of E⁻ᵗ(e^η)", fontsize=11, loc="left")
import json
flow = json.load(open("/mnt/project-files/math-theory/solitons-rungs/results.json"))["flow"]
ts = sorted(float(k) for k in flow); mis = [flow[str(t) if str(t) in flow else f"{t}"]["resid"] for t in ts]
axs[1].semilogy(ts, mis, color=INK, lw=1.8)
axs[1].scatter([1.0], [flow["1.0"]["resid"]], color=SEA, s=50, zorder=3)
axs[1].set_xlabel("rung step t"); axs[1].set_title("mismatch with a KdV soliton", fontsize=11, loc="left")
axs[1].annotate("t = 1", (1.0, flow["1.0"]["resid"]), xytext=(10, 4), textcoords="offset points", color=SEA)
save(fig, "flow.svg")

# 5. two roads from rung 0 to rung 1 --------------------------------------------------------
x = np.linspace(0, 1, 300)
fig, axs = plt.subplots(1, 2, figsize=(9.2, 4.2))
for ax, rung, sol, title, meets in [
        (axs[0], x*(3*x - 1)/2, x**3, "three bodies", [0, 0.5, 1]),
        (axs[1], x*(2*x - 1)*(3*x - 1)/2, x**6, "four bodies", [0, 1])]:
    ax.plot(x, rung, color=SAND, lw=2.2, label="rung flow, rung s")
    ax.plot(x, sol, color=SEA, lw=2.2, label="soliton weight A")
    ax.scatter(meets, [m**3 if title == "three bodies" else m**6 for m in meets], color=INK, s=28, zorder=3)
    ax.axhline(0, color=GREY, lw=0.6)
    ax.set_xlabel("two-body weight (s = A)"); ax.set_title(f"corner weight, {title}", fontsize=11, loc="left")
axs[0].legend(frameon=False, fontsize=9, loc="upper left")
axs[0].annotate("rung 0", (0, 0), xytext=(4, 10), textcoords="offset points", fontsize=9)
axs[0].annotate("rung 1", (1, 1), xytext=(-40, -4), textcoords="offset points", fontsize=9)
save(fig, "roads.svg")

# 6. the soliton paths are the creases of a max ---------------------------------------------
k1, k2 = 1.0, 1.6; A = Aij(k1, k2); ph = math.log(A)
X, Tt = np.meshgrid(np.linspace(-22, 22, 700), np.linspace(-7, 7, 420))
e1 = k1*X - k1**3*Tt; e2 = k2*X - k2**3*Tt
H = np.stack([0*X, e1, e2, e1 + e2 + ph])
lt = np.log(np.exp(H - H.max(0)).sum(0)) + H.max(0)
u = 2*np.gradient(np.gradient(lt, X[0], axis=1), X[0], axis=1)
fig, ax = plt.subplots(figsize=(9.2, 4.6))
ax.imshow(u, extent=[-22, 22, -7, 7], origin="lower", aspect="auto", cmap="Blues", vmin=0, vmax=1.4)
dom = H.argmax(0)
ax.contour(X, Tt, dom, levels=[0.5, 1.5, 2.5], colors=[ROSE], linewidths=1.1, linestyles="--")
ax.set_xlabel("position x"); ax.set_ylabel("time t")
ax.set_title("Two KdV solitons (blue) and the creases of max(0, η₁, η₂, η₁+η₂+log A) (red)", fontsize=11, loc="left")
save(fig, "tropical.svg")
print("figures written")
