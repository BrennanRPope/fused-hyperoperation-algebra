"""Pictures for the Klein-bottle flow write-up. Run after ensemble.py and upward.py: python3 figures.py"""
import glob, json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle, Circle
from flow import Grid, project

os.makedirs("figures", exist_ok=True)
RUNS = "runs_N128"
COL = {"T2": "#7a7a7a", "Trect": "#1f77b4", "K": "#d62728"}
NAME = {"T2": "square torus (double size)", "Trect": "torus (rectangle, straight glue)", "K": "Klein bottle (rectangle, flipped glue)"}
plt.rcParams.update({"font.size": 11, "axes.titlesize": 12})


def load_runs(d=RUNS):
    runs = {s: [] for s in COL}
    for fn in sorted(glob.glob(f"{d}/*_s*.json")):
        base = os.path.basename(fn)[:-5]
        s, seed = base.split("_s")
        r = json.load(open(fn)); r["seed"] = int(seed)
        runs[s].append(r)
    for s in runs:
        runs[s].sort(key=lambda r: r["seed"])
    return runs


def vlim(w):
    m = np.percentile(np.abs(w), 99.5)
    return -m, m


def half(w):
    """The fundamental rectangle [0,pi] x [0,2pi] of a cover field (axis 0 = x)."""
    return w[: w.shape[0] // 2, :]


def show(ax, w, rect=True, title=None, cmap="RdBu_r", lim=None):
    f = half(w) if rect else w
    lo, hi = lim if lim is not None else vlim(w)
    ext = [0, np.pi if rect else 2 * np.pi, 0, 2 * np.pi]
    ax.imshow(f.T, origin="lower", extent=ext, cmap=cmap, vmin=lo, vmax=hi, interpolation="bilinear")
    ax.set_xticks([]); ax.set_yticks([])
    if title:
        ax.set_title(title)


def glue_marks(ax, surface, rect=True):
    """Arrows on the edges of the fundamental rectangle: single arrows = top/bottom glued straight,
    double arrows = left/right edges, drawn opposite on K to show the flip."""
    W = np.pi if rect else 2 * np.pi
    kw = dict(arrowstyle="-|>", mutation_scale=14, lw=2)
    # top & bottom: straight on both
    for yy in (0, 2 * np.pi):
        ax.add_patch(FancyArrowPatch((0.3 * W, yy), (0.7 * W, yy), color="k", **kw))
    if not rect:
        return
    up = ((0, 2.2), (0, 4.1))
    ax.add_patch(FancyArrowPatch(*up, color="#2ca02c", **kw))
    if surface == "K":
        ax.add_patch(FancyArrowPatch((W, 4.1), (W, 2.2), color="#2ca02c", **kw))
    else:
        ax.add_patch(FancyArrowPatch((W, 2.2), (W, 4.1), color="#2ca02c", **kw))


# ---------------------------------------------------------------- Figure 1: the two gluings
def fig1():
    fig, axs = plt.subplots(1, 3, figsize=(13, 5.2), gridspec_kw=dict(width_ratios=[1, 1, 1.6]))
    for ax, s in zip(axs[:2], ["Trect", "K"]):
        ax.add_patch(Rectangle((0, 0), np.pi, 2 * np.pi, fc="#f4f4f4", ec="k", lw=1.5))
        glue_marks(ax, s)
        # a small spinning vortex walking out of the right edge
        ax.add_patch(Circle((2.7, 3.15), 0.32, fc="#ffbbbb", ec="#d62728", lw=2))
        ax.annotate("", xy=(2.7 + 0.32, 3.15 + 0.15), xytext=(2.7 + 0.25, 3.15 - 0.25),
                    arrowprops=dict(arrowstyle="->", color="#d62728", lw=2, connectionstyle="arc3,rad=0.6"))
        yy = 3.15 if s == "Trect" else 2 * np.pi - 3.15
        c = "#ffbbbb" if s == "Trect" else "#bbccff"
        ec = "#d62728" if s == "Trect" else "#1f4fbf"
        ax.add_patch(Circle((0.45, yy), 0.32, fc=c, ec=ec, lw=2, ls="--"))
        rad = 0.6 if s == "Trect" else -0.6
        ax.annotate("", xy=(0.45 + 0.32, yy + 0.15), xytext=(0.45 + 0.25, yy - 0.25),
                    arrowprops=dict(arrowstyle="->", color=ec, lw=2, connectionstyle=f"arc3,rad={rad}"))
        ax.text(np.pi / 2, -0.75, "comes back spinning the same way" if s == "Trect" else "comes back spinning the OTHER way",
                ha="center", fontsize=10, color=ec)
        ax.set_xlim(-0.5, np.pi + 0.5); ax.set_ylim(-1.1, 2 * np.pi + 0.4); ax.set_aspect("equal"); ax.axis("off")
        ax.set_title("Torus: green edges glued straight" if s == "Trect" else "Klein bottle: green edges glued flipped")
    ax = axs[2]
    ax.add_patch(Rectangle((0, 0), np.pi, 2 * np.pi, fc="#ffecec", ec="k", lw=1.5))
    ax.add_patch(Rectangle((np.pi, 0), np.pi, 2 * np.pi, fc="#ecf0ff", ec="k", lw=1.5, ls="--"))
    ax.text(np.pi / 2, 5.6, "the rectangle", ha="center")
    ax.text(1.5 * np.pi, 5.6, "its mirror copy\n(shifted, flipped)", ha="center")
    ax.add_patch(Circle((1.0, 1.6), 0.35, fc="#ffbbbb", ec="#d62728", lw=2))
    ax.add_patch(Circle((1.0 + np.pi, 2 * np.pi - 1.6), 0.35, fc="#bbccff", ec="#1f4fbf", lw=2))
    ax.annotate("", xy=(1.0 + np.pi - 0.3, 2 * np.pi - 1.75), xytext=(1.25, 1.75),
                arrowprops=dict(arrowstyle="->", lw=1.5, connectionstyle="arc3,rad=-0.25"))
    ax.text(np.pi + 0.15, 2.0, "glide: move right\nby half the width,\nflip top to bottom;\nspin reverses", fontsize=10)
    ax.text(np.pi, -0.75, "Simulation trick: run on the double-size square and keep only flows that\nmatch their glide image "
            "(velocity copied, spin reversed). That IS flow on the Klein bottle.", ha="center", fontsize=10)
    ax.set_xlim(-0.3, 2 * np.pi + 0.3); ax.set_ylim(-1.4, 2 * np.pi + 0.4); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title("How the Klein bottle is simulated")
    fig.tight_layout()
    fig.savefig("figures/fig1-two-gluings.png", dpi=130); plt.close(fig)


# ---------------------------------------------------------------- Figure 2: upward flow
def fig2():
    d = np.load("upward.npz")
    ts = [0, 6, 12, 20, 40]
    fig = plt.figure(figsize=(15, 8.6))
    gs = fig.add_gridspec(2, 6, width_ratios=[1, 1, 1, 1, 1, 2.0])
    g = Grid(d["K_w_0"].shape[0])
    xs = g.x[: g.N // 2]
    for j, s in enumerate(["Trect", "K"]):
        lim = vlim(d["K_w_0"]) if s == "K" else (-6, 6)
        for i, t in enumerate(ts):
            ax = fig.add_subplot(gs[j, i])
            show(ax, d[f"{s}_w_{t:g}"], lim=lim)
            u, v = half(d[f"{s}_u_{t:g}"]), half(d[f"{s}_v_{t:g}"])
            sk = 16
            ax.quiver(xs[::sk // 2], g.x[::sk], u[::sk // 2, ::sk].T, v[::sk // 2, ::sk].T, color="k", scale=25, width=0.012)
            glue_marks(ax, s); ax.set_xlim(-0.15, np.pi + 0.15); ax.set_ylim(-0.15, 2 * np.pi + 0.15)
            ax.set_title(f"t = {t}")
            if i == 0:
                ax.set_ylabel("torus" if s == "Trect" else "Klein bottle", fontsize=13)
    ax = fig.add_subplot(gs[:, 5])
    t = d["t"]
    ax.plot(t, d["Trect_E"], color=COL["Trect"], lw=2.5, label="torus: never slows down")
    ax.plot(t, d["K_E"], color=COL["K"], lw=2.5, label="Klein bottle: flow tears itself up")
    ax.set_xlabel("time"); ax.set_ylabel("kinetic energy per area"); ax.set_ylim(0, 0.55); ax.legend(loc="lower left")
    ax.set_title("Same start: everything moving up")
    fig.suptitle("Fluid moving up through the rectangle. Straight glue: a steady sweep forever. Flipped glue: up meets down at the "
                 "twist, the shear rolls up, and it ends as ONE vortex", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig("figures/fig2-upward-flow-tears-on-the-twist.png", dpi=110); plt.close(fig)


# ---------------------------------------------------------------- Figure 3: matched decaying turbulence
def fig3(seed=0):
    z = {s: np.load(f"{RUNS}/{s}_s{seed}_snaps.npz") for s in COL}
    ts = ["0.0", "5.0", "20.0", "60.0", "400.0"]
    fig, axs = plt.subplots(3, 5, figsize=(14, 12), gridspec_kw=dict(height_ratios=[1, 1, 1]))
    for j, s in enumerate(["T2", "Trect", "K"]):
        for i, t in enumerate(ts):
            w = z[s][f"w_{t}"]
            show(axs[j, i], w, rect=(s != "T2"), title=f"t = {float(t):g}")
            if s != "T2":
                glue_marks(axs[j, i], s)
                axs[j, i].set_xlim(-0.15, np.pi + 0.15); axs[j, i].set_ylim(-0.15, 2 * np.pi + 0.15)
        axs[j, 0].set_ylabel(NAME[s], fontsize=11)
    fig.suptitle("Same random start, same viscosity. Red spins one way, blue the other. Small eddies look alike;\n"
                 "the end states differ: the square torus ends with a pair, the torus with a band, the Klein bottle with ONE vortex",
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig("figures/fig3-matched-turbulence.png", dpi=100); plt.close(fig)


# ---------------------------------------------------------------- final-state classification
def angle(r, s):
    """Mixing angle between the two families of lowest modes: 0 = pure x-shear family, 90 = pure y family."""
    if s == "K":
        A, B = r["ax"][-1], abs(r["siny"][-1])
    elif s == "T2":
        A, B = r["ax"][-1], r["ay"][-1]
    else:
        A, B = r["ax"][-1], r["ay"][-1]
    return float(np.degrees(np.arctan2(B, A)))


def classify(theta):
    return "band" if (theta < 20 or theta > 70) else "vortex"


def fig4(runs, seed=0):
    fig, axs = plt.subplots(1, 4, figsize=(15, 5.4), gridspec_kw=dict(width_ratios=[2, 1, 1, 2]))
    for ax, s in zip([axs[0], axs[1], axs[2]], ["T2", "Trect", "K"]):
        z = np.load(f"{RUNS}/{s}_s{seed}_snaps.npz")
        psi = z["psi_400.0"]; w = z["w_400.0"]
        rect = s != "T2"
        show(ax, w, rect=rect)
        g = Grid(psi.shape[0])
        P = half(psi) if rect else psi
        xs = g.x[: P.shape[0]]
        ax.contour(xs, g.x, P.T, 12, colors="k", linewidths=0.8)
        if rect:
            glue_marks(ax, s); ax.set_xlim(-0.15, np.pi + 0.15); ax.set_ylim(-0.15, 2 * np.pi + 0.15)
        th = angle(runs[s][seed], s)
        ax.set_title({"T2": "square torus: a vortex PAIR", "Trect": "torus: a BAND", "K": "Klein bottle: ONE vortex"}[s]
                     + f"\n(mix angle {th:.0f}°)")
    ax = axs[3]
    z = np.load(f"{RUNS}/K_s{seed}_snaps.npz")
    show(ax, z["w_400.0"], rect=False)
    ax.axvline(np.pi, color="k", ls="--")
    ax.text(np.pi / 2, 6.6, "the Klein bottle", ha="center"); ax.text(1.5 * np.pi, 6.6, "its mirror copy", ha="center")
    ax.set_ylim(0, 7.1)
    ax.set_title("The lone vortex and its mirror image:\non the cover they look like a pair")
    fig.suptitle("End states at t = 400 (vorticity colours, streamlines black)", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig("figures/fig4-end-states.png", dpi=110); plt.close(fig)


def fig5(runs, runs3=None):
    sets = [("smaller starting eddies (k0 = 6)", runs)] + ([("bigger starting eddies (k0 = 3)", runs3)] if runs3 else [])
    fig, AX = plt.subplots(len(sets), 2, figsize=(14, 5.2 * len(sets)), squeeze=False)
    out = {}
    for row, (label, runs) in enumerate(sets):
        out[label] = _fig5_row(AX[row], runs, label)
    fig.tight_layout(); fig.savefig("figures/fig5-end-state-counts.png", dpi=110); plt.close(fig)
    return out


def _fig5_row(axs, runs, label):
    ax = axs[0]
    summary = {}
    for i, s in enumerate(["T2", "Trect", "K"]):
        th = np.array([angle(r, s) for r in runs[s]])
        summary[s] = dict(angles=th.round(1).tolist(), vortex=int(sum(classify(t) == "vortex" for t in th)), n=len(th),
                          shell=[round(r["shell_frac"][-1], 3) for r in runs[s]])
        ax.scatter(th, np.full_like(th, i) + np.random.default_rng(i).uniform(-0.15, 0.15, len(th)), s=70,
                   color=COL[s], edgecolor="k", zorder=3)
    ax.axvspan(20, 70, color="#ffe9b0", zorder=0)
    ax.text(45, 2.45, "vortex states\n(pair on the square torus,\none vortex on the Klein bottle)", ha="center", fontsize=10)
    ax.text(8, 2.45, "band", ha="center", fontsize=10); ax.text(82, 2.45, "band", ha="center", fontsize=10)
    ax.set_yticks([0, 1, 2]); ax.set_yticklabels(["square torus", "torus", "Klein bottle"])
    ax.set_xlim(-5, 95); ax.set_ylim(-0.5, 2.9); ax.set_xlabel("mix angle of the two largest-scale mode families (degrees)")
    ax.set_title(f"How {len(runs['K'])} matched runs per surface end: {label}")
    ax = axs[1]
    labels = ["square torus", "torus", "Klein bottle"]
    vort = [summary[s]["vortex"] for s in ["T2", "Trect", "K"]]
    band = [summary[s]["n"] - summary[s]["vortex"] for s in ["T2", "Trect", "K"]]
    ax.bar(labels, vort, color="#f0b000", label="ends as vortex state")
    ax.bar(labels, band, bottom=vort, color="#9ecae1", label="ends as band")
    ax.set_ylabel("runs"); ax.legend(); ax.set_title("Counts")
    return summary


def fig6(runs):
    fig, axs = plt.subplots(1, 3, figsize=(16, 5))
    for s in ["T2", "Trect", "K"]:
        t = np.array(runs[s][-1]["t"])
        E = np.array([np.interp(t, r["t"], r["E"]) for r in runs[s]])
        Z = np.array([np.interp(t, r["t"], r["Z"]) for r in runs[s]])
        axs[0].plot(t, E.mean(0), color=COL[s], lw=2.5, label=NAME[s])
        axs[0].fill_between(t, E.min(0), E.max(0), color=COL[s], alpha=0.15)
        axs[1].plot(t, Z.mean(0), color=COL[s], lw=2.5, label=NAME[s])
        axs[1].fill_between(t, Z.min(0), Z.max(0), color=COL[s], alpha=0.15)
        for tk, ls in (("5.0", "-"), ("20.0", "--")):
            sp = np.array([r["spectra"][tk] for r in runs[s]]).mean(0)
            k = np.arange(len(sp))
            axs[2].loglog(k[1:], sp[1:], color=COL[s], ls=ls, lw=2, label=f"{NAME[s].split(' (')[0]}, t={float(tk):g}")
    axs[0].set_xscale("log"); axs[0].set_xlabel("time"); axs[0].set_title("Energy (per area)"); axs[0].legend(fontsize=9)
    axs[1].set_xscale("log"); axs[1].set_yscale("log"); axs[1].set_xlabel("time"); axs[1].set_title("Enstrophy (spin squared, per area)")
    axs[2].set_xlabel("wavenumber k (1 = largest scale)"); axs[2].set_title("Energy by scale: the twist only shows at the largest scales")
    axs[2].set_ylim(1e-9, 1); axs[2].legend(fontsize=8)
    fig.tight_layout(); fig.savefig("figures/fig6-energy-enstrophy-spectra.png", dpi=110); plt.close(fig)


def fig7(runs, seed=0):
    fig, axs = plt.subplots(1, 3, figsize=(16, 5))
    for s in ["Trect", "K"]:
        z = np.load(f"{RUNS}/{s}_s{seed}_snaps.npz")
        w = z["w_20.0"].ravel(); w = w / w.std()
        axs[0].hist(w, bins=121, range=(-8, 8), density=True, histtype="step", lw=2.5, color=COL[s], log=True, label=NAME[s])
        for r in runs[s]:
            axs[1].plot(r["t"], r["skew"], color=COL[s], alpha=0.6, lw=1.3)
    axs[0].set_xlabel("vorticity / its spread"); axs[0].set_title("Spin distribution at t = 20 (one run)"); axs[0].legend(fontsize=9)
    axs[1].set_xscale("log"); axs[1].set_xlabel("time"); axs[1].set_title("Skew of the spin (all runs)\nKlein bottle: exactly 0, always")
    axs[1].plot([], [], color=COL["Trect"], label="torus"); axs[1].plot([], [], color=COL["K"], label="Klein bottle"); axs[1].legend()
    # 'total spin' by cut
    ax = axs[2]
    for s in ["Trect", "K"]:
        z = np.load(f"{RUNS}/{s}_s{seed}_snaps.npz")
        w = z["w_20.0"]; N = w.shape[0]
        a = np.arange(N) * 2 * np.pi / N
        vals = [np.sum(w[(np.arange(N // 2) + i) % N, :]) * (2 * np.pi / N) ** 2 for i in range(N)]
        ax.plot(a, vals, color=COL[s], lw=2.5, label=NAME[s])
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xlabel("where the surface is cut open (left edge of the rectangle)")
    ax.set_title("'Total spin' of the fluid: torus 0 wherever you cut;\nKlein bottle depends on the cut (so it is not a real quantity)")
    ax.legend(fontsize=9)
    fig.tight_layout(); fig.savefig("figures/fig7-spin-symmetry-and-cut.png", dpi=110); plt.close(fig)


def fig8():
    """Allowed large-scale flows on each surface."""
    g = Grid(64); X, Y = g.X, g.Y
    cases = [("torus", "uniform sweep up", 0 * X, np.ones_like(X), True),
             ("torus", "uniform sweep right", np.ones_like(X), 0 * X, True),
             ("torus", "band  u = cos y", np.cos(Y), 0 * X, True),
             ("Klein bottle", "uniform sweep right", np.ones_like(X), 0 * X, True),
             ("Klein bottle", "uniform sweep up: NOT allowed", 0 * X, np.ones_like(X), False),
             ("Klein bottle", "band  v = sin x", 0 * X, np.sin(X), True),
             ("Klein bottle", "band  u = cos y", np.cos(Y), 0 * X, True),
             ("Klein bottle", "lone vortex  psi = cos x + sin y", np.cos(Y), np.sin(X), True)]
    fig, axs = plt.subplots(2, 5, figsize=(15, 7.5))
    order = [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (1, 3), (1, 4)]
    for (r, c), (surf, title, u, v, ok) in zip(order, cases):
        ax = axs[r, c]
        sl = slice(0, 32)
        ax.add_patch(Rectangle((0, 0), np.pi, 2 * np.pi, fc="#f7f7f7" if ok else "#ffe0e0", ec="k"))
        ax.quiver(g.x[sl][::3], g.x[::4], u[sl][::3, ::4].T, v[sl][::3, ::4].T, color="k" if ok else "#c00000", scale=18, width=0.012)
        glue_marks(ax, "K" if surf == "Klein bottle" else "Trect")
        ax.set_xlim(-0.3, np.pi + 0.3); ax.set_ylim(-0.3, 2 * np.pi + 0.3); ax.set_aspect("equal"); ax.axis("off")
        ax.set_title(f"{surf}\n{title}", fontsize=10, color="k" if ok else "#c00000")
    axs[0, 3].axis("off"); axs[0, 4].axis("off")
    axs[0, 3].text(0.0, 0.6, "Torus: 2 sweeps that never slow\n(Betti number 2); largest-scale\nmodes: bands along one side only.\n\n"
                   "Klein bottle: 1 sweep (Betti number 1),\nbut 3 largest-scale modes, which mix\ninto bands or one lone vortex.",
                   fontsize=12, va="center")
    fig.suptitle("Which large-scale flows each surface allows (same rectangle; green edges glued straight on the torus, flipped on the Klein bottle)")
    fig.tight_layout(rect=[0, 0, 1, 0.95]); fig.savefig("figures/fig8-allowed-large-flows.png", dpi=110); plt.close(fig)


def fig9():
    """Mode counting: Klein-bottle modes = (torus-cover modes + twisted trace)/2. The twisted trace stays small:
    the twist only changes the count at the largest scales (Weyl's law is set by area alone)."""
    lam_max = 400
    cnt = {}
    for kx in range(-21, 22):
        for ky in range(-21, 22):
            lam = kx * kx + ky * ky
            if lam == 0 or lam > lam_max:
                continue
            cnt.setdefault(lam, []).append((kx, ky))
    lams = sorted(cnt)
    cover, klein, rect, twist = [], [], [], []
    for lam in lams:
        modes = cnt[lam]
        n = len(modes)  # complex exponentials = real dimension of the shell
        # pseudo-scalar glide action on e^{i(kx x + ky y)}: -> -(-1)^kx e^{i(kx x - ky y)}; trace = sum over fixed (ky = 0)
        tr = sum(-((-1) ** kx) for kx, ky in modes if ky == 0)
        trs = sum(((-1) ** kx) for kx, ky in modes)  # straight shift: diagonal, trace (-1)^kx
        cover.append(n); twist.append(tr); klein.append((n + tr) // 2); rect.append((n + trs) // 2)
    cover, klein, rect, twist = map(np.array, (cover, klein, rect, twist))
    fig, axs = plt.subplots(1, 2, figsize=(14, 5))
    ax = axs[0]
    L = np.array(lams)
    ax.step(L, np.cumsum(rect), where="post", color=COL["Trect"], lw=2.5, label="torus (rectangle)")
    ax.step(L, np.cumsum(klein), where="post", color=COL["K"], lw=2, ls="--", label="Klein bottle (same rectangle)")
    ax.plot(L, L * (2 * np.pi**2) / (4 * np.pi), color="k", lw=1, label="Weyl's law: area × λ / 4π")
    ax.set_xlabel("eigenvalue λ (wavenumber squared)"); ax.set_ylabel("number of flow patterns up to λ")
    ax.set_title("Same area, same count of patterns at small scales"); ax.legend()
    ax = axs[1]
    m = L <= 30
    ax.bar(L[m] - 0.2, rect[m], width=0.4, color=COL["Trect"], label="torus")
    ax.bar(L[m] + 0.2, klein[m], width=0.4, color=COL["K"], label="Klein bottle")
    ax.set_xlabel("eigenvalue λ"); ax.set_ylabel("patterns at exactly λ")
    ax.set_title("Shell by shell they differ by a little that never grows;\nthe one that matters is the largest scale, λ = 1: torus 2, Klein bottle 3"); ax.legend()
    fig.tight_layout(); fig.savefig("figures/fig9-mode-counting.png", dpi=110); plt.close(fig)
    return dict(lam=lams[:8], rect=rect[:8].tolist(), klein=klein[:8].tolist(), twist_trace=twist[:8].tolist())


if __name__ == "__main__":
    runs = load_runs()
    fig1(); fig2(); fig3(); fig4(runs)
    runs3 = load_runs("runs_N128_k3") if os.path.isdir("runs_N128_k3") else None
    summary = fig5(runs, runs3)
    fig6(runs); fig7(runs); fig8()
    counts = fig9()
    json.dump(dict(end_states=summary, mode_counts=counts), open("results.json", "w"), indent=1)
    print(json.dumps(dict(end_states=summary, mode_counts=counts), indent=1))
