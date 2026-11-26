"""Pictures for the shape-dials page (DRAFT). Run: python3 figures.py"""
import sys, os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "..", "simplex"))
from facedial import *
from shapes_catalog import REGULAR, NONREGULAR
from wythoff import CHAINS, wythoff_polytope, orbit, roots_for, chamber_corners
import shape_dial as SD
FIG = os.path.join(HERE, "figures")
os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({"font.size": 10, "axes.titlesize": 10})
PCOL = {2: "#1f77b4", 3: "#2a9d8f", 4: "#57a773", 6: "#e9c46a", 10: "#f4a261", 25: "#e76f51", np.inf: "#222222"}
PLIST = [2, 3, 4, 6, 10, 25, np.inf]


def plab(p):
    return "p = ∞ (the polytope)" if np.isinf(p) else ("p = 2 (round)" if p == 2 else f"p = {p}")


def rad(f, U, p):
    return f.radial_polytope(U) if np.isinf(p) else f.radial(U, p)


# ----------------------------------------------------------- 3D surface helper
NU, NV = 140, 70
uu, vv = np.meshgrid(np.linspace(0, 2 * np.pi, NU), np.linspace(1e-3, np.pi - 1e-3, NV))
DIRS = np.stack([np.sin(vv) * np.cos(uu), np.sin(vv) * np.sin(uu), np.cos(vv)], -1).reshape(-1, 3)
LIGHT = np.array([0.5, -0.6, 0.8]); LIGHT /= np.linalg.norm(LIGHT)


def surf(ax, R, color):
    P = (R[:, None] * DIRS).reshape(NV, NU, 3)
    X, Y, Z = P[..., 0], P[..., 1], P[..., 2]
    gu = np.gradient(P, axis=1); gv = np.gradient(P, axis=0)
    nrm = np.cross(gu, gv); nrm /= np.linalg.norm(nrm, axis=-1, keepdims=True) + 1e-12
    sh = np.clip(np.abs(nrm @ LIGHT), 0.15, 1)
    base = np.array(matplotlib.colors.to_rgb(color))
    fc = np.clip(base[None, None, :] * (0.35 + 0.75 * sh[..., None]), 0, 1)
    ax.plot_surface(X, Y, Z, facecolors=fc, linewidth=0, antialiased=False, shade=False, rcount=NV, ccount=NU)
    m = R.max() * 0.62
    ax.set_xlim(-m, m); ax.set_ylim(-m, m); ax.set_zlim(-m, m)
    ax.set_box_aspect((1, 1, 1)); ax.set_axis_off(); ax.view_init(22, 35)


def norm_r(R):
    return R / R.max()


# ----------------------------------------------------------- fig1: polygons
def fig1():
    names = ["3-gon", "4-gon", "5-gon", "6-gon", "8-gon"]
    fig, axs = plt.subplots(1, 5, figsize=(15, 3.6))
    th, U = circle_dirs(1500)
    for ax, nm in zip(axs, names):
        f = FaceDial(REGULAR[nm])
        for p in PLIST:
            r = rad(f, U, p)
            ax.plot(*(r[:, None] * U).T, color=PCOL[p], lw=1.8 if p in (2, np.inf) else 1.1, label=plab(p))
        ax.set_aspect("equal"); ax.axis("off")
        ax.set_title(f"{nm}\noffset a = {f.a:.3f}" + (" (symmetric)" if f.a == 0 else ""))
    axs[0].legend(loc="upper left", bbox_to_anchor=(0, -0.02), ncol=7, frameon=False, fontsize=9)
    fig.suptitle("Figure 1. One dial p, aimed at each polygon's faces: round at p = 2, the polygon at p = ∞", y=1.02)
    fig.savefig(os.path.join(FIG, "fig1-polygons-round-to-sharp.png"), dpi=130, bbox_inches="tight"); plt.close(fig)


def grid3d(rows, ps, fname, title, colors, getter):
    fig = plt.figure(figsize=(2.6 * len(ps) + 1.6, 2.5 * len(rows)))
    for i, nm in enumerate(rows):
        for j, p in enumerate(ps):
            ax = fig.add_subplot(len(rows), len(ps), i * len(ps) + j + 1, projection="3d")
            surf(ax, getter(nm, p), colors[i % len(colors)])
            if i == 0:
                ax.set_title(plab(p), fontsize=11)
            if j == 0:
                ax.text2D(-0.12, 0.5, nm, transform=ax.transAxes, rotation=90, va="center", ha="center", fontsize=11)
    fig.suptitle(title, y=1.0, fontsize=12)
    fig.subplots_adjust(wspace=0.0, hspace=0.02)
    fig.savefig(os.path.join(FIG, fname), dpi=110, bbox_inches="tight"); plt.close(fig)


CACHE = {}


def dial_for(nm):
    if nm not in CACHE:
        CACHE[nm] = FaceDial({**REGULAR, **NONREGULAR}[nm], nm)
    return CACHE[nm]


def fig2():
    rows = ["tetrahedron", "cube", "octahedron", "dodecahedron", "icosahedron"]
    grid3d(rows, [2, 4, 10, np.inf], "fig2-platonic-round-to-sharp.png",
           "Figure 2. The five Platonic solids, each grown out of the same ball by turning p",
           ["#e76f51", "#2a9d8f", "#457b9d", "#e9c46a", "#8d6cab"], lambda nm, p: rad(dial_for(nm), DIRS, p))


# ----------------------------------------------------------- fig3: 4D regular, central slices
def fig3():
    names = ["5-cell", "tesseract", "16-cell", "24-cell", "120-cell", "600-cell"]
    rng = np.random.default_rng(11)
    Q, _ = np.linalg.qr(rng.standard_normal((4, 2)))
    th, U = circle_dirs(1500)
    D = U @ Q.T
    fig, axs = plt.subplots(1, 6, figsize=(17, 3.4))
    for ax, nm in zip(axs, names):
        f = dial_for(nm)
        for p in PLIST:
            r = rad(f, D, p)
            ax.plot(*(r[:, None] * U).T, color=PCOL[p], lw=1.8 if p in (2, np.inf) else 1.0)
        ax.set_aspect("equal"); ax.axis("off")
        ax.set_title(f"{nm}\n{f.m} faces, a = {f.a:g}")
    handles = [plt.Line2D([], [], color=PCOL[p], lw=2) for p in PLIST]
    fig.legend(handles, [plab(p) for p in PLIST], loc="lower center", ncol=7, frameon=False, bbox_to_anchor=(0.5, -0.06))
    fig.suptitle("Figure 3. The six regular 4D shapes, cut by the same flat plane through the centre, p from round to sharp", y=1.03)
    fig.savefig(os.path.join(FIG, "fig3-4d-regular-slices.png"), dpi=130, bbox_inches="tight"); plt.close(fig)


# ----------------------------------------------------------- fig4: non-regular 3D + elliptope
def ell_radials(p, a=0.5):
    out = np.empty(len(DIRS))
    for k, u in enumerate(DIRS):
        lo, hi = 0.0, 4.0
        for _ in range(45):
            mid = (lo + hi) / 2
            x, y, z = mid * u
            lam = np.linalg.eigvalsh(np.array([[1, x, y], [x, 1, z], [y, z, 1.0]]))
            t = np.abs(a + 1 - lam)
            val = t.max() if np.isinf(p) else np.mean(t ** p) ** (1 / p)
            lo, hi = (mid, hi) if val <= a + 1 else (lo, mid)
        out[k] = lo
    return out


def fig4():
    rows = ["permutohedron (3D)", "associahedron (3D, Loday)", "zonotope (3D, 5 random zones)", "cyclic polytope C(7,3)",
            "elliptope"]
    def get(nm, p):
        if nm == "elliptope":
            return ell_radials(p)
        return rad(dial_for(nm), DIRS, p)
    grid3d(rows, [2, 4, 10, np.inf], "fig4-nonregular-round-to-sharp.png",
           "Figure 4. Shapes we use that are not regular. Elliptope: the faces are its matrix's eigenvalues (curved sides)",
           ["#2a9d8f", "#e76f51", "#457b9d", "#8d6cab", "#c9a227"], get)


# ----------------------------------------------------------- fig5: Wythoff depth (shape-choice dial) x p
WCACHE = {}


def wy(nm, t):
    k = (nm, round(t, 4))
    if k not in WCACHE:
        WCACHE[k] = FaceDial(wythoff_polytope(nm, t))
    return WCACHE[k]


WNAMES = {("B3", 0): "cube", ("B3", 0.5): "truncated cube", ("B3", 1): "cuboctahedron", ("B3", 1.5): "truncated octahedron",
          ("B3", 2): "octahedron", ("A3", 0): "tetrahedron", ("A3", 0.5): "truncated tetrahedron", ("A3", 1): "octahedron",
          ("A3", 1.5): "truncated tetrahedron\n(dual side)", ("A3", 2): "tetrahedron (dual)", ("H3", 0): "dodecahedron",
          ("H3", 0.5): "truncated dodecahedron", ("H3", 1): "icosidodecahedron", ("H3", 1.5): "truncated icosahedron",
          ("H3", 2): "icosahedron"}


def fig5():
    fams, ts = ["A3", "B3", "H3"], [0, 0.5, 1, 1.5, 2]
    cols = {"A3": "#e76f51", "B3": "#2a9d8f", "H3": "#8d6cab"}
    fig = plt.figure(figsize=(14.5, 9))
    for i, fam in enumerate(fams):
        for j, t in enumerate(ts):
            ax = fig.add_subplot(3, 5, i * 5 + j + 1, projection="3d")
            f = wy(fam, t)
            D = DIRS @ np.eye(3)
            surf(ax, f.radial_polytope(D), cols[fam])
            ax.set_title(f"depth {t}: {WNAMES[(fam, t)]}", fontsize=9.5)
            if j == 0:
                ax.text2D(-0.1, 0.5, {"A3": "simplex family\n(A3 mirrors)", "B3": "cube family\n(B3 mirrors)",
                                       "H3": "dodecahedron family\n(H3 mirrors)"}[fam], transform=ax.transAxes,
                          rotation=90, va="center", ha="center", fontsize=10.5)
    fig.suptitle("Figure 5. The shape-choice dial: truncation depth (Wythoff's kaleidoscope). Whole depth = regular or "
                 "rectified; half depth = the truncations between. Shown at p = ∞", y=0.99, fontsize=11.5)
    fig.subplots_adjust(wspace=0.0, hspace=0.12)
    fig.savefig(os.path.join(FIG, "fig5-truncation-depth-dial.png"), dpi=110, bbox_inches="tight"); plt.close(fig)


def fig6():
    ts, ps = [0, 0.5, 1, 1.5, 2], [2, 3, 6, np.inf]
    fig = plt.figure(figsize=(13.5, 10.5))
    for i, p in enumerate(ps[::-1]):
        for j, t in enumerate(ts):
            ax = fig.add_subplot(4, 5, i * 5 + j + 1, projection="3d")
            f = wy("B3", t)
            surf(ax, rad(f, DIRS, p), "#2a9d8f" if not np.isinf(p) else "#1d6f66")
            if i == 0:
                ax.set_title(f"depth {t}\n{WNAMES[('B3', t)]}", fontsize=10)
            if j == 0:
                ax.text2D(-0.1, 0.5, plab(p), transform=ax.transAxes, rotation=90, va="center", ha="center", fontsize=10.5)
    fig.suptitle("Figure 6. Two dials on the cube family: depth picks the polytope (left to right), p sets how sharp "
                 "(top to bottom).\nEvery column meets the same ball at p = 2", y=0.995, fontsize=11.5)
    fig.subplots_adjust(wspace=0.0, hspace=0.06)
    fig.savefig(os.path.join(FIG, "fig6-two-dial-map-cube-family.png"), dpi=110, bbox_inches="tight"); plt.close(fig)


def fig7():
    ts = [0, 0.15, 0.3, 0.5, 0.7, 0.85, 1]
    fig, axs = plt.subplots(3, 7, figsize=(14, 6.4))
    th, U = circle_dirs(1200)
    for i, n in enumerate([3, 4, 5]):
        for j, t in enumerate(ts):
            ax = axs[i, j]
            f = FaceDial(wythoff_polytope("I2(n)", t, n))
            for p in [2, 4, np.inf]:
                Uo = U
                r = rad(f, Uo, p)
                ax.plot(*(r[:, None] * U).T, color=PCOL[p] if p != 4 else "#e9c46a", lw=1.6 if np.isinf(p) else 1)
            ax.set_aspect("equal"); ax.axis("off")
            nv = len(f.Vw)
            lab = {0: f"{n}-gon", 0.5: f"regular {2 * n}-gon", 1: f"{n}-gon (turned)"}.get(t, f"{nv} corners")
            ax.set_title(f"depth {t}: {lab}", fontsize=9)
        axs[i, 0].text(-0.25, 0.5, f"n = {n}", transform=axs[i, 0].transAxes, fontsize=11, va="center")
    handles = [plt.Line2D([], [], color=c, lw=2) for c in [PCOL[2], "#e9c46a", PCOL[np.inf]]]
    fig.legend(handles, ["p = 2", "p = 4", "p = ∞"], loc="lower center", ncol=3, frameon=False)
    fig.suptitle("Figure 7. In 2D: depth is continuous (n-gon → regular 2n-gon → n-gon); the side count n is a whole number", y=1.0)
    fig.savefig(os.path.join(FIG, "fig7-polygons-depth-and-side-count.png"), dpi=125, bbox_inches="tight"); plt.close(fig)


# ----------------------------------------------------------- fig8: simplex test
def fig8():
    R = json.load(open(os.path.join(HERE, "results.json")))
    tri = FaceDial(regular_polygon(3))
    th, U = circle_dirs(4000)
    fig, axs = plt.subplots(1, 3, figsize=(15.5, 4.8))
    ax = axs[0]
    rfp = tri.radial(U, 6); rfp /= rfp.max()
    best = R["fp_vs_temperature"]["6"][0]
    rT = lse_radial(tri, U, best); rT /= rT.max()
    ax.plot(*(rfp[:, None] * U).T, color="#1f77b4", lw=2, label="face-power, p = 6 (a = 1/2)")
    ax.plot(*(rT[:, None] * U).T, color="#e76f51", lw=1.4, ls="--", label=f"temperature route, β = {best:.2f}")
    rp = tri.radial_polytope(U); rp /= rp.max()
    ax.plot(*(rp[:, None] * U).T, color="#999", lw=1, label="triangle (p = ∞)")
    ax.set_aspect("equal"); ax.axis("off"); ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.2), frameon=False, fontsize=9)
    ax.set_title("A. Face-power and temperature at matched sharpness:\nnearly the same shape (radius off by at most 1.1%)")
    ax = axs[1]
    kfp = R["khat_facepower"]
    x = [v[0] for v in kfp.values()]; y = [v[1] for v in kfp.values()]
    ax.plot(x, y, "o-", color="#1f77b4", label="face-power (p = 2 … 200)")
    O = np.zeros((1, 2))
    xs, ys = [], []
    for b in [0.01, 0.1, 0.3, 1, 2, 3, 6, 10, 30, 100]:
        r = lse_radial(tri, U, b)
        xs.append(r.min() / r.max()); ys.append(SD.khat(SD.Body(r[:, None] * U), O)[0])
    ax.plot(xs, ys, "s--", color="#e76f51", label="temperature (β = 0.01 … 100)")
    for name, mk, bodies in [("blend (disk + triangle)", "^", [SD.minkowski(t) for t in [0.25, 0.5, 0.75, 0.9]]),
                             ("trim (triangle ∩ disk)", "v", [SD.cut_triangle(s) for s in [0.55, 0.65, 0.75, 0.9]])]:
        xx, yy = [], []
        for bo in bodies:
            rr = np.linalg.norm(bo.boundary, axis=1)
            xx.append(rr.min() / rr.max()); yy.append(SD.khat(bo, O)[0])
        ax.plot(xx, yy, mk + ":", color="#777", label=name)
    ax.set_xlabel("in-to-out radius (1 = round, 0.5 = triangle)"); ax.set_ylabel("centre turn K̂ (−1 = round, 0 = flat)")
    ax.invert_xaxis(); ax.legend(fontsize=8.5, frameon=False); ax.grid(alpha=0.3)
    ax.set_title("B. Centre turn along four round-to-triangle routes")
    ax = axs[2]
    lim = R["limit"]
    a = np.array(list(map(float, lim.keys()))); e = np.array(list(lim.values()))
    ax.loglog(a, e, "o-", color="#1f77b4", label="face-power with offset a, p = 3a")
    ax.loglog(a, e[0] * a[0] / a, "k:", label="slope −1 (error ∝ 1/a)")
    ax.set_xlabel("offset a (the triangle's own a = 1/2)"); ax.set_ylabel("largest radius gap to temperature β = 3")
    ax.legend(frameon=False, fontsize=9); ax.grid(alpha=0.3, which="both")
    ax.set_title("C. The temperature route is face-power's far edge:\nlet a → ∞ with p = β·a")
    fig.suptitle("Figure 8. The simplex test: the round-to-triangle routes against face-power", y=1.03)
    fig.savefig(os.path.join(FIG, "fig8-simplex-test.png"), dpi=125, bbox_inches="tight"); plt.close(fig)


def fig9():
    R = json.load(open(os.path.join(HERE, "results.json")))["sharpness_gap"]
    PS = np.array([2, 3, 4, 8, 16, 64, 256])
    fig, ax = plt.subplots(figsize=(8.5, 5.2))
    cmap = plt.get_cmap("viridis")
    items = list(R.items())
    for k, (nm, g) in enumerate(items):
        g = np.array(g)
        ax.loglog(PS[1:], g[1:], color=cmap(k / len(items)), lw=1, alpha=0.85)
    for nm in ["3-gon", "cube", "600-cell", "cyclic polytope C(8,4)"]:
        g = np.array(R[nm]); ax.loglog(PS[1:], g[1:], lw=2.2, label=nm)
    ax.loglog(PS[1:], 2 / PS[1:], "k:", label="slope −1 (gap ∝ 1/p)")
    ax.set_xlabel("p"); ax.set_ylabel("how far K_p still bulges past the polytope (largest radius ratio − 1)")
    ax.set_title("Figure 9. All 34 shapes sharpen at the same 1/p rate (bound: (m^{1/p} − 1)(a + 1))")
    ax.legend(frameon=False); ax.grid(alpha=0.3, which="both")
    fig.savefig(os.path.join(FIG, "fig9-sharpening-rate.png"), dpi=130, bbox_inches="tight"); plt.close(fig)


if __name__ == "__main__":
    which = sys.argv[1:] or ["1", "2", "3", "4", "5", "6", "7", "8", "9"]
    for w in which:
        globals()["fig" + w]()
        print("fig", w, "done", flush=True)
