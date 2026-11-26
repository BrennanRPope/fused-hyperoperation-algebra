"""Pictures for the radial-plane note (DRAFT). Run after checks.py (reads results.json).
Writes figures/g1, g2, g3. Colour key as in the curved-space note: orange = known, blue = ours, grey = context."""
import os, json, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.special import lambertw

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures"); os.makedirs(FIG, exist_ok=True)
RES = json.load(open(os.path.join(HERE, "results.json")))
BLUE, ORANGE, GREY, INK, INK2, GRID = "#257de0", "#f47230", "#9a9a9a", "#191919", "#666666", "#e6e6e6"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
                     "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 2,
                     "legend.frameon": False, "axes.titleweight": "bold", "axes.titlesize": 11})


def save(fig, name):
    fig.tight_layout(rect=[0, 0, 1, 0.9])
    fig.savefig(os.path.join(FIG, name), dpi=150, bbox_inches="tight"); plt.close(fig); print("wrote", name)


rstar = lambda r, M: r + 2 * M * np.log(r / (2 * M) - 1)
r_of = lambda rs, M: 2 * M * (1 + lambertw(np.exp(rs / (2 * M) - 1)).real)


def path(M, r0, E, tau, n=200):
    f = lambda r: 1 - 2 * M / r
    rd0 = -math.sqrt(E**2 - f(r0))
    rhs = lambda s, y: [y[1], -M / y[0]**2, E / f(y[0])]
    out = {}
    for sgn in (-1, 1):
        s = np.linspace(0, sgn * tau / 2, n)
        out[sgn] = solve_ivp(rhs, [0, sgn * tau / 2], [r0, rd0, 0.0], t_eval=s, rtol=1e-11, atol=1e-12).y
    return out


# ---------- g1: the spherical diamond and the cancellation ----------
M, r0, E, tau = 1.0, 6.0, 1.0, 4.0
p = path(M, r0, E, tau)
rp, tp = p[-1][0, -1], p[-1][2, -1]
rq, tq = p[1][0, -1], p[1][2, -1]
up, vp = tp - rstar(rp, M), tp + rstar(rp, M)
uq, vq = tq - rstar(rq, M), tq + rstar(rq, M)
fig, ax = plt.subplots(1, 2, figsize=(11, 4.8), gridspec_kw={"width_ratios": [1, 1.15]})
a = ax[0]
for (u0, u1, v0, v1) in [(up, uq, vp, vp), (up, uq, vq, vq), (up, up, vp, vq), (uq, uq, vp, vq)]:
    uu, vv = np.linspace(u0, u1, 200), np.linspace(v0, v1, 200)
    a.plot(r_of((vv - uu) / 2, M), (uu + vv) / 2, color=INK2, lw=1.4)
U, V = np.meshgrid(np.linspace(up, uq, 9)[1:-1], np.linspace(vp, vq, 9)[1:-1])
R, T = r_of((V - U) / 2, M), (U + V) / 2
a.scatter(R, T, s=2.2 * R**2, facecolor="none", edgecolor=BLUE, lw=1, zorder=3)
path_r = np.concatenate([p[-1][0][::-1], p[1][0]]); path_t = np.concatenate([p[-1][2][::-1], p[1][2]])
a.plot(path_r, path_t, color=ORANGE, lw=2.2, label="river parcel falling in (E = 1)")
a.plot([rp, rq], [tp, tq], "o", color=ORANGE, ms=6)
a.text(rp + 0.15, tp - 0.15, "p", color=INK, fontsize=11); a.text(rq + 0.15, tq + 0.05, "q", color=INK, fontsize=11)
a.axvline(2 * M, color=GREY, ls="--", lw=1.2); a.text(2 * M + 0.08, tq + 0.05, "horizon r = 2M", color=GREY, fontsize=9)
a.set_xlabel("radius r (each point is a whole sphere)"); a.set_ylabel("Schwarzschild time t")
a.set_title("A spherical diamond in the radial plane")
a.text(0.06, 0.2, "circle size = sphere area 4πr²:\nsmaller toward the hole", transform=a.transAxes, va="top",
       fontsize=9, color=BLUE)
a.legend(loc="upper right", fontsize=9); a.set_xlim(1.6, 9.5)

b = ax[1]
mr = 0.2                                   # M / r for the black-hole case; river parcel rdot^2 = 2M/r
cases = [("hovering, empty space", 1.0, 0.0), ("river parcel, black hole\n(M/r = 0.2)", 1 + 2 * (2 * mr) - 2 * mr, -2 * mr)]
x = np.arange(len(cases)); w = 0.26
for i, (lab, sph, dial) in enumerate(cases):
    b.bar(i - w, sph, w, color=BLUE, label="sphere sizes: 1 + 2ṙ² − 2M/r" if i == 0 else None)
    b.bar(i, dial, w, color=ORANGE, label="radial-plane curvature: −2M/r" if i == 0 else None)
    b.bar(i + w, sph + dial, w, color=INK2, label="total = 2E² − 1" if i == 0 else None)
    for dx, val in [(-w, sph), (0, dial), (w, sph + dial)]:
        b.text(i + dx, val + (0.04 if val >= 0 else -0.12), f"{val:+.1f}" if val else "0", ha="center", fontsize=9)
b.axhline(0, color=INK2, lw=0.8)
b.set_xticks(x); b.set_xticklabels([c[0] for c in cases])
b.set_ylabel("first correction to the count  (× τ²/24r²)"); b.set_ylim(-0.6, 1.75)
b.set_title("The mass shows twice and cancels")
b.legend(loc="upper left", fontsize=9)
fig.suptitle("g1  Counting spheres in a diamond: the mass bends both parts, and the two bends cancel (ours)",
             fontsize=12, fontweight="bold", y=0.97)
save(fig, "g1_spherical_diamond.png")

# ---------- g2: measured coefficients against the predictions ----------
rows = RES["leading_rows"]
fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
a = ax[0]
Es = np.linspace(0.75, 1.4, 200)
a.plot(Es, 2 * Es**2 - 1, color=BLUE, label="prediction 2E² − 1 (no M in it)")
for rw in rows:
    flat = rw["M"] == 0
    a.plot(rw["E"], rw["4D coeff x 24 r0^2"], "s" if flat else "o", ms=8, mfc="white" if flat else BLUE,
           mec=INK if flat else BLUE, zorder=3)
a.plot([], [], "s", mfc="white", mec=INK, label="empty space (M = 0)")
a.plot([], [], "o", color=BLUE, label="black hole, M/r0 from 0.04 to 0.17")
a.set_xlabel("free-fall energy E"); a.set_ylabel("4D count coefficient × 24 r0²")
a.set_title("Sphere-weighted count: only E shows"); a.legend(fontsize=9, loc="upper left")
b = ax[1]
Ks = np.array([-rw["-K"] for rw in rows]); c2 = np.array([rw["2D coeff x 24"] for rw in rows])
kk = np.linspace(0, Ks.max() * 1.1, 50)
b.plot(kk * 1e3, -kk * 1e3, color=ORANGE, label="prediction −K, K = 2M/r0³ (known: BD)")
b.plot(Ks * 1e3, c2 * 1e3, "o", color=BLUE, ms=8, label="measured, 8 cases (any E)")
b.set_xlabel("K = 2M/r0³  (× 10⁻³)"); b.set_ylabel("2D count coefficient × 24  (× 10⁻³)")
b.set_title("Plain radial-plane count: reads the local curvature"); b.legend(fontsize=9)
fig.suptitle("g2  What a small diamond's count reads, measured against the predictions (checks D1, D2)",
             fontsize=12, fontweight="bold", y=0.97)
save(fig, "g2_coefficients.png")

# ---------- g3: the mass appears at the next order ----------
nxt = RES["next_order_E1"]; fit = RES["next_order_fit"]
xs = np.array([float(k) for k in nxt]); y2 = np.array([nxt[k][0] for k in nxt]); y4 = np.array([nxt[k][1] for k in nxt])
xx = np.linspace(0, 0.21, 100)
fig, ax = plt.subplots(figsize=(7.5, 4.4))
ax.plot(xx, (fit["2D"][0] * xx + fit["2D"][1] * xx**2) * 1e3, color=ORANGE, lw=1.4)
ax.plot(xs, y2 * 1e3, "o", color=ORANGE, ms=8, label="plain radial-plane count")
ax.plot(xx, (fit["4D"][0] * xx + fit["4D"][1] * xx**2) * 1e3, color=BLUE, lw=1.4)
ax.plot(xs, y4 * 1e3, "o", color=BLUE, ms=8, label="sphere-weighted count")
ax.axhline(0, color=INK2, lw=0.8)
ax.set_xlabel("M / r0  (river parcel, E = 1)"); ax.set_ylabel("next-order coefficient of (τ/r0)⁴  (× 10⁻³)")
ax.legend(fontsize=9, loc="lower left")
ax.text(0.006, -0.12, "zero in empty space", color=INK2, fontsize=9)
fig.suptitle("g3  One order further out, the mass shows up directly (check D6, fitted)", fontsize=12,
             fontweight="bold", y=0.97)
save(fig, "g3_next_order.png")
