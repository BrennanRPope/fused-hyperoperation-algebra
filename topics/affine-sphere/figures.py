import numpy as np, json, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa
import os
plt.rcParams.update({"font.size": 11, "axes.titlesize": 12.5, "figure.dpi": 130, "savefig.bbox": "tight"})
EC, PC, HC, GREY, GREEN, PURPLE = "#1f5fbf", "#e0115f", "#e07b00", "#8a8f99", "#2e8b57", "#6a3d9a"
os.makedirs("figures", exist_ok=True)
R2 = np.sqrt(2)


def col(c):
    return PC if abs(c) < 1e-9 else (EC if c > 0 else HC)


def surface(k, d, n=60, rmax=1.6):
    """Chord surface x^2+kxy+dy^2+w^2=y, the sheet through O, as (x, w, y) grids."""
    c = d - k * k / 4
    th = np.linspace(0, 2 * np.pi, n)
    if c > 1e-9:
        a, b, y0 = 1 / (2 * np.sqrt(c)), 1 / (2 * c), 1 / (2 * c)
        ph = np.linspace(0, np.pi, n)
        PH, TH = np.meshgrid(ph, th)
        X, W, Y = a * np.sin(PH) * np.cos(TH), a * np.sin(PH) * np.sin(TH), y0 - b * np.cos(PH)
    else:
        rho = np.linspace(0, rmax, n)
        RH, TH = np.meshgrid(rho, th)
        X, W = RH * np.cos(TH), RH * np.sin(TH)
        if abs(c) < 1e-9:
            Y = X ** 2 + W ** 2
        else:
            ac = -c
            y0 = 1 / (2 * c)
            Y = y0 + np.sqrt((1 / (4 * ac) + X ** 2 + W ** 2) / ac)
    return X - k * Y / 2, W, Y


def centre(k, d):
    c = d - k * k / 4
    return None if abs(c) < 1e-9 else np.array([-k, 2.0]) / (4 * c)


# ---------- Fig 1: one quadric through O for every setting, with its affine normal at O
cases = [(0.6, 0.7, "round side (Δ < 0)\nELLIPSOID"), (1.0, 0.25, "wall (Δ = 0)\nPARABOLOID"),
         (1.0, -0.5, "two-wall side (Δ > 0)\nHYPERBOLOID (2 sheets, this one)")]
fig = plt.figure(figsize=(15, 5.6))
for i, (k, d, lab) in enumerate(cases):
    c = d - k * k / 4
    ax = fig.add_subplot(1, 3, i + 1, projection="3d")
    X, W, Y = surface(k, d, rmax=1.05 if abs(c) < 1e-9 else 1.6)
    ax.plot_surface(X, W, Y, color=col(c), alpha=.28, linewidth=0, rstride=2, cstride=2)
    ax.plot_wireframe(X, W, Y, color=col(c), alpha=.25, linewidth=.4, rstride=6, cstride=6)
    xi = R2 * np.array([-k / 2, 1, 0])
    ax.quiver(0, 0, 0, xi[0], xi[2], xi[1], color="k", lw=2.4, arrow_length_ratio=.15)
    ax.scatter([0], [0], [0], color="k", s=30)
    ax.text(0.08, 0, -0.18, "O (identity)", fontsize=10)
    ax.text(xi[0] * 1.05, 0, xi[1] * 1.08, "affine normal ξ(O)\n= √2·(−k/2, 1, 0)", fontsize=9.5)
    ce = centre(k, d)
    if ce is not None:
        ax.scatter([ce[0]], [0], [ce[1]], color=GREEN, s=40, marker="*")
        ax.plot([0, ce[0]], [0, 0], [0, ce[1]], color=GREEN, ls="--", lw=1.2)
        ax.text(ce[0], 0, ce[1] + .1, "centre", color=GREEN, fontsize=10)
    ax.set_title(f"{lab}\nk = {k}, d = {d}:  H = 2^(3/2)(d − k²/4) = {2 ** 1.5 * c:.2f}", color=col(c), fontsize=11.5)
    ax.set_xlabel("x"); ax.set_ylabel("w"); ax.set_zlabel("y")
    lim = 1.6
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(-.4, 2.8 if c >= 0 else 2.8)
    ax.view_init(elev=14, azim=-62)
fig.suptitle("Figure 1. Every plane-dial setting (k, d) gives one affine sphere through the identity: "
             "x² + kxy + dy² + w² = y", fontsize=13.5, y=1.0)
fig.savefig("figures/fig1-one-affine-sphere-per-setting.png"); plt.close(fig)


# ---------- Fig 2: affine normals all meet at the centre (the w = 0 slice is the chord conic)
def Pc(sv, k, d):
    q = 1 + k * sv + d * sv * sv
    return np.array([sv / q, sv * sv / q])


fig, axs = plt.subplots(1, 3, figsize=(15, 5.2))
for ax, (k, d, lab) in zip(axs, cases):
    c = d - k * k / 4
    lim = 2.4
    xs = np.linspace(-lim, lim, 500); XX, YY = np.meshgrid(xs, np.linspace(-1, 3.8, 500))
    ax.contour(XX, YY, XX ** 2 + k * XX * YY + d * YY ** 2 - YY, [0], colors=[col(c)], linewidths=2.5)
    ce = centre(k, d)
    svals = np.linspace(-2.5, 2.5, 11) if c > 0 else (np.linspace(-6, 6, 13) if abs(c) < 1e-9 else np.concatenate([np.linspace(-6, 6, 25), [-12, 12, 30, -30]]))
    for j, sv in enumerate(svals):
        q = 1 + k * sv + d * sv * sv
        if q <= 0.05:
            continue
        p = Pc(sv, k, d)
        if np.abs(p).max() > 3.6:
            continue
        # Euclidean normal: gradient of x^2+kxy+dy^2-y
        g = np.array([2 * p[0] + k * p[1], k * p[0] + 2 * d * p[1] - 1]); g /= np.linalg.norm(g)
        ax.plot([p[0] - 3 * g[0], p[0] + 3 * g[0]], [p[1] - 3 * g[1], p[1] + 3 * g[1]], color=GREY, lw=.6, alpha=.6)
        if ce is not None:
            dirn = ce - p
        else:
            dirn = np.array([-k / 2, 1.0])
        dirn = dirn / np.linalg.norm(dirn)
        if ce is not None and c > 0:
            ax.plot([p[0], ce[0]], [p[1], ce[1]], color="k", lw=.9)
        else:
            ax.plot([p[0] - 6 * dirn[0], p[0] + 6 * dirn[0]], [p[1] - 6 * dirn[1], p[1] + 6 * dirn[1]], color="k", lw=.9)
        ax.plot(*p, "o", color=col(c), ms=4)
    if ce is not None:
        ax.plot(*ce, "*", color=GREEN, ms=16, zorder=6)
        ax.annotate("all affine normals\nmeet here (centre)", ce, xytext=(ce[0] + .5, ce[1] + (1.6 if c > 0 else .5)), color=GREEN, fontsize=10,
                    arrowprops=dict(arrowstyle="->", color=GREEN), bbox=dict(facecolor="white", edgecolor="none", alpha=.9))
    else:
        ax.text(-2.25, 3.25, "all affine normals parallel\n(centre at infinity)", color=GREEN, fontsize=10,
                bbox=dict(facecolor="white", edgecolor="none", alpha=.9))
    ax.plot(0, 0, "ko", ms=6); ax.text(.08, -.3, "O", fontsize=11)
    ax.set_xlim(-lim, lim); ax.set_ylim(-1, 3.8); ax.set_aspect("equal")
    ax.set_title(f"{lab.splitlines()[0]}: k = {k}, d = {d}", color=col(c))
axs[0].plot([], [], color="k", lw=.9, label="affine normals (black)")
axs[0].plot([], [], color=GREY, lw=.8, label="ordinary normals (grey) do not meet")
axs[0].legend(loc="lower left", fontsize=9)
fig.suptitle("Figure 2. The w = 0 slice is the chord conic x² + kxy + dy² = y. Its affine normals meet in one point: "
             "that is what 'affine sphere' means", fontsize=13, y=1.02)
fig.savefig("figures/fig2-affine-normals-meet.png"); plt.close(fig)

# ---------- Fig 3 and 4 from measured.json
M = json.load(open("measured.json"))
k = np.array([m["k"] for m in M]); d = np.array([m["d"] for m in M]); H = np.array([m["H"] for m in M])
xi0 = np.array([m["xi0"] for m in M])
c = d - k * k / 4
fig, axs = plt.subplots(1, 2, figsize=(13, 5.2))
sc = axs[0].scatter(d, H, c=k, cmap="coolwarm", s=55, edgecolor="k", vmin=-2, vmax=2)
axs[0].set_xlabel("curvature dial d (a setting)"); axs[0].set_ylabel("measured affine mean curvature H")
axs[0].set_title("H against d alone: scattered (k matters too)")
axs[1].scatter(c, H, c=k, cmap="coolwarm", s=55, edgecolor="k", vmin=-2, vmax=2)
cc = np.linspace(c.min() - .1, c.max() + .1, 10)
axs[1].plot(cc, 2 ** 1.5 * cc, color="k", lw=1.2, label="H = 2^(3/2)·(d − k²/4)")
axs[1].set_xlabel("curvature reading d − k²/4 = −Δ/4 (computed from k and d)")
axs[1].set_title("H against the curvature reading: one line")
for ax in axs:
    ax.axhline(0, color=GREY, lw=.6); ax.axvline(0, color=GREY, lw=.6)
axs[1].text(.03, .93, "ellipsoids (H > 0)", transform=axs[1].transAxes, color=EC, weight="bold")
axs[1].text(.22, .03, "two-sheeted hyperboloids (H < 0)", transform=axs[1].transAxes, color=HC, weight="bold")
axs[1].legend(loc="lower right", bbox_to_anchor=(1, .12))
cb = fig.colorbar(sc, ax=axs, shrink=.85); cb.set_label("homothety k")
fig.suptitle(f"Figure 3. Measured on {len(M)} random settings: the shape's mean curvature reads d − k²/4 and nothing else",
             fontsize=13, y=1.0)
fig.savefig("figures/fig3-mean-curvature-reads-curvature.png"); plt.close(fig)

fig, axs = plt.subplots(1, 3, figsize=(17, 5))
fig.subplots_adjust(wspace=.42)
lean = -xi0[:, 0] / xi0[:, 1]
sc = axs[0].scatter(k, lean, c=d, cmap="viridis", s=55, edgecolor="k")
kk = np.linspace(-2.1, 2.1, 10); axs[0].plot(kk, kk / 2, color="k", lw=1.1, label="lean = k/2")
axs[0].set_xlabel("homothety k (a setting)"); axs[0].set_ylabel("lean of ξ(O): −ξx/ξy")
axs[0].set_title("lean of the normal at O reads k"); axs[0].legend(loc="upper left")
fig.colorbar(sc, ax=axs[0], shrink=.8).set_label("curvature dial d")
axs[1].scatter(d, np.linalg.norm(xi0, axis=1) / np.sqrt(1 + lean ** 2), c=k, cmap="coolwarm", s=55, edgecolor="k", vmin=-2, vmax=2)
axs[1].axhline(R2, color="k", lw=1.1); axs[1].set_ylim(1.2, 1.6)
axs[1].text(-.9, R2 + .02, "ξy(O) = √2 for every setting", fontsize=10)
axs[1].set_xlabel("curvature dial d"); axs[1].set_ylabel("height ξy(O)")
axs[1].set_title("the normal at O carries no d")
# shear panel
ax = axs[2]
kk0, dd0 = 1.2, 0.6
cc0 = dd0 - kk0 ** 2 / 4
xs = np.linspace(-2.2, 2.2, 400); XX, YY = np.meshgrid(xs, np.linspace(-.6, 3.6, 400))
ax.contour(XX, YY, XX ** 2 + kk0 * XX * YY + dd0 * YY ** 2 - YY, [0], colors=[EC], linewidths=2.2)
ax.contour(XX, YY, XX ** 2 + cc0 * YY ** 2 - YY, [0], colors=[EC], linewidths=2.2, linestyles="--")
ax.annotate("", (-kk0 / 2 * R2, R2), (0, 0), arrowprops=dict(arrowstyle="-|>", color="k", lw=2))
ax.annotate("", (0, R2), (0, 0), arrowprops=dict(arrowstyle="-|>", color=GREY, lw=2))
ax.annotate("", (0.1, 2.35), (-kk0 / 2 * 2.35 + .1, 2.35), arrowprops=dict(arrowstyle="->", color=PURPLE, lw=1.4))
ax.text(-1.95, 3.0, "shear x → x + (k/2)y\n(area-preserving; on the curve it is\nthe Cayley change of reading)",
        color=PURPLE, fontsize=9.5, bbox=dict(facecolor="white", edgecolor="none", alpha=.9))
ax.text(-2.1, .9, "k = 1.2, d = 0.6\n(solid)", color=EC, fontsize=9.5, bbox=dict(facecolor="white", edgecolor="none", alpha=.9))
ax.text(.75, 1.0, "k = 0, same\nd − k²/4 = 0.24\n(dashed)", color=EC, fontsize=9.5, bbox=dict(facecolor="white", edgecolor="none", alpha=.9))
ax.set_aspect("equal"); ax.set_xlim(-2.2, 2.2); ax.set_ylim(-.6, 3.6)
ax.set_title("a shear stands the normal up, removing k")
fig.suptitle("Figure 4. The homothety k is read off where the normal at the identity leans; "
             "a shear removes it, so k is placement, not shape", fontsize=13, y=1.02)
fig.savefig("figures/fig4-normal-lean-reads-k.png"); plt.close(fig)

# ---------- Fig 5: where a third (cubic-form) dial can live
fig, ax = plt.subplots(figsize=(12, 7))
kk = np.linspace(-3, 3, 400)
ax.fill_between(kk, kk ** 2 / 4, 3.6, color=EC, alpha=.10)
ax.fill_between(kk, -2.4, kk ** 2 / 4, color=HC, alpha=.13)
ax.plot(kk, kk ** 2 / 4, color=PC, lw=3)
ax.axhline(0, color="k", lw=.5); ax.axvline(0, color="k", lw=.5)
ax.set_xlim(-3, 3); ax.set_ylim(-2.4, 3.6)
ax.set_xlabel("homothety k"); ax.set_ylabel("curvature dial d")
ax.text(-2.9, 3.3, "round side: a complete elliptic affine sphere is an ELLIPSOID\n→ cubic form C = 0, no room for a third dial",
        color=EC, fontsize=11, weight="bold", va="top")
ax.annotate("wall d = k²/4: a complete improper affine sphere is a PARABOLOID\n→ C = 0 (Jörgens, Calabi, Pogorelov)", (2.3, 1.3225),
            xytext=(.15, 2.2), color=PC, fontsize=10.5, weight="bold", arrowprops=dict(arrowstyle="->", color=PC))
ax.text(-2.9, -0.15, "two-wall side: complete hyperbolic affine spheres ↔ convex cones\n"
        "(Cheng–Yau). C can be nonzero, one sphere per cone shape:\nthe only place a cubic-form (shape) dial has room.",
        color=HC, fontsize=10, weight="bold", va="top")
ins = ax.inset_axes([.70, .02, .29, .38], projection="3d")
s_ = np.linspace(-1.3, 1.3, 40); S_, T_ = np.meshgrid(s_, s_)
ins.plot_surface(np.exp(S_), np.exp(T_), np.exp(-S_ - T_), color=HC, alpha=.5, linewidth=0)
ins.set_xlim(0, 4); ins.set_ylim(0, 4); ins.set_zlim(0, 4)
ins.set_title("xyz = 1 (Țițeica): triangle cone,\n|K|² = 0.877 ≠ 0 (check S8b)", fontsize=9)
ins.set_xticks([]); ins.set_yticks([]); ins.set_zticks([])
ins2 = ax.inset_axes([.42, .01, .2, .21], projection="3d")
X, W, Y = surface(1.0, -0.5, n=40, rmax=1.2)
ins2.plot_surface(X, W, Y, color=HC, alpha=.4, linewidth=0)
ins2.set_title("our shape here: hyperboloid,\nround cone, C = 0", fontsize=9)
ins2.set_xticks([]); ins2.set_yticks([]); ins2.set_zticks([])
ax.set_title("Figure 5. Our sphere is a quadric for every setting (C = 0, check S5). "
             "Known theorems leave room for a nonzero cubic form only on the two-wall side", fontsize=12.5)
fig.savefig("figures/fig5-where-a-cubic-form-could-live.png"); plt.close(fig)

# ---------- Fig 6: the matrix-orbit alternative
fig = plt.figure(figsize=(13, 5.4))
lab6 = [(1.0, "det X₀ = d − k²/4 > 0 (round side)\n→ TWO-sheeted hyperboloid:\nhyperbolic type", EC),
        (0.0, "= 0 (wall)\n→ cone (not smooth)", PC),
        (-1.0, "< 0 (two-wall side)\n→ ONE-sheeted hyperboloid:\nindefinite, not convex", HC)]
for i, (lv, lab, cl) in enumerate(lab6):
    ax = fig.add_subplot(1, 3, i + 1, projection="3d")
    p_ = np.linspace(-1.8, 1.8, 50); th = np.linspace(0, 2 * np.pi, 60)
    if lv > 0:   # two sheets: r = +-sqrt(p^2+q^2+lv)
        Pg, TH = np.meshgrid(np.linspace(0, 2.4, 50), th)
        A, B = Pg * np.cos(TH), Pg * np.sin(TH)
        Rp = np.sqrt(A ** 2 + B ** 2 + lv)
        ax.plot_surface(A, B, Rp, color=cl, alpha=.4, linewidth=0)
        ax.plot_surface(A, B, -Rp, color=cl, alpha=.4, linewidth=0)
    else:        # one connected piece (or cone): p^2+q^2 = r^2 - lv
        Rg, TH = np.meshgrid(np.linspace(-2.5, 2.5, 80), th)
        rad = np.sqrt(Rg ** 2 - lv)
        ax.plot_surface(rad * np.cos(TH), rad * np.sin(TH), Rg, color=cl, alpha=.4, linewidth=0)
    ax.set_title(lab, color=cl, fontsize=10.5); ax.set_zlim(-2.6, 2.6); ax.view_init(elev=8, azim=-60)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
    ax.set_xlabel("p"); ax.set_ylabel("q"); ax.set_zlabel("r")
fig.suptitle("Figure 6. The other candidate: the generator matrix X = [[k, 1], [−d, 0]]. Its trace is k and its traceless part X₀ lies on det X₀ = d − k²/4.\n"
             "Same split of k from the curvature, but the shapes have the wrong types, so this is not the affine sphere", fontsize=12, y=1.04)
fig.savefig("figures/fig6-matrix-orbits-wrong-types.png"); plt.close(fig)
print("ok")
