"""DRAFT figures for hypersphere/orbits/README.md. Run: python3 figures_orbits.py (numpy, scipy, matplotlib)."""
import json
import math
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import integrate

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "figures")
os.makedirs(OUT, exist_ok=True)
KNOWN, OURS, GREY, RED = "#d9822b", "#2b6cb0", "#888888", "#c53030"
plt.rcParams.update({"font.size": 11, "axes.titlesize": 12})
res = json.load(open(os.path.join(HERE, "results.json")))


def save(fig, name):
    fig.savefig(os.path.join(OUT, name), dpi=140, bbox_inches="tight")
    plt.close(fig)
    print("wrote", name)


def run_sr(c, v0, T):
    def rhs(t, y):
        x, yy, px, py = y
        r = math.hypot(x, yy)
        g = math.sqrt(1 + (px * px + py * py) / c ** 2)
        return [px / g, py / g, -x / r ** 3, -yy / r ** 3]
    return integrate.solve_ivp(rhs, (0, T), [1, 0, 0, v0], rtol=1e-11, atol=1e-12, max_step=0.005)


# ------------------------------------------------------------ O1 two spheres
c = 3.0
sol = run_sr(c, 1.1, 60)
x, y, px, py = sol.y
r = np.hypot(x, y)
H = c * c * np.sqrt(1 + (px ** 2 + py ** 2) / c ** 2) - 1 / r
p0 = math.sqrt(-2 * (H[0] - c * c))
p2 = px ** 2 + py ** 2
W = np.c_[(p2 - p0 ** 2) / (p2 + p0 ** 2), 2 * p0 * px / (p2 + p0 ** 2), 2 * p0 * py / (p2 + p0 ** 2)]
fig = plt.figure(figsize=(16, 5.4))
ax = fig.add_subplot(1, 3, 1)
ax.plot(x, y, color=OURS, lw=0.8)
ax.plot(0, 0, "k*", ms=12)
ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
ax.set_title("A. A relativistic orbit (c = 3 in these units):\nthe ellipse turns a little every lap")
ax = fig.add_subplot(1, 3, 2, projection="3d")
u, v = np.mgrid[0:2 * np.pi:40j, 0:np.pi:20j]
ax.plot_wireframe(np.cos(u) * np.sin(v), np.sin(u) * np.sin(v), np.cos(v), color=GREY, alpha=0.1, lw=0.5)
gx = (x * px + y * py)
peri = np.where((gx[:-1] < 0) & (gx[1:] >= 0))[0]
cols = plt.cm.Blues(np.linspace(0.4, 1, len(peri)))
for (i0, i1), col in zip(zip(peri[:-1], peri[1:]), cols):
    ax.plot(W[i0:i1, 1], W[i0:i1, 2], W[i0:i1, 0], color=col, lw=1.6)
ax.quiver(0, 0, 0, 0, 0, 1.3, color=RED, lw=2, arrow_length_ratio=0.1)
ax.text(0, 0, 1.4, "energy axis", color=RED)
ax.set_axis_off(); ax.set_box_aspect((1, 1, 1))
ax.set_title("B. Fock's momentum sphere:\neach lap is a great circle, tilting lap by lap")
ax = fig.add_subplot(1, 3, 3)
L = x[0] * py[0] - y[0] * px[0]
R = 1 / L
rho = math.atanh(R / c)
t = np.linspace(0, 2 * np.pi, 300)
ax.plot(np.cos(t), np.sin(t), color=GREY, ls=":")
rk = math.tanh(rho / 2)        # Poincaré-disc radius of a rapidity circle rho
ax.fill(rk * np.cos(t), rk * np.sin(t), color=OURS, alpha=0.35)
ax.text(-0.95, -1.25, f"disc of speed-radius k/L = {R / c:.2f} c\narea × (1/c²) = {2 * math.pi * (math.cosh(rho) - 1):.4f}"
        f" = the turn per lap ({res['sr_c3.0']['adv']:.4f})", fontsize=10)
ax.text(-0.35, 0.9, "edge = light speed", color=GREY, fontsize=9)
ax.set_aspect("equal"); ax.axis("off"); ax.set_xlim(-1.1, 1.1); ax.set_ylim(-1.45, 1.1)
ax.set_title("C. Relativity's velocity space (imaginary\nradius c, drawn as a disc): the turn is this area")
fig.suptitle("Figure O1. Two spheres: the orbit lives on Fock's sphere; how much it turns per lap is an area on "
             "relativity's velocity sphere", y=1.02)
save(fig, "figO1-two-spheres.png")

# ------------------------------------------------------------ O2 advance vs area, many L
fig, axs = plt.subplots(1, 2, figsize=(13, 4.6))
Ls, advs = [], []
for v0 in np.linspace(0.75, 1.3, 8):
    s = run_sr(2.5, v0, 80)
    xx, yy, ppx, ppy = s.y
    g = np.sqrt(1 + (ppx ** 2 + ppy ** 2) / 6.25)
    rr = np.hypot(xx, yy)
    rd = (xx * ppx / g + yy * ppy / g) / rr
    idx = np.where((rd[:-1] < 0) & (rd[1:] >= 0))[0]
    ang = np.unwrap([math.atan2(yy[i], xx[i]) for i in idx])
    Ls.append(xx[0] * ppy[0] - yy[0] * ppx[0]); advs.append(np.diff(ang).mean())
Ls = np.array(Ls)
LL = np.linspace(Ls.min() * 0.95, Ls.max() * 1.05, 200)
axs[0].plot(LL, 2 * np.pi * (1 / np.sqrt(1 - (1 / (2.5 * LL)) ** 2) - 1), color=OURS,
            label="area of the velocity disc, radius k/L")
axs[0].plot(Ls, advs, "o", color=KNOWN, label="simulated orbits (c = 2.5)")
axs[0].set_xlabel("angular momentum L"); axs[0].set_ylabel("perihelion advance per lap (rad)"); axs[0].legend()
axs[0].set_title("A. Eight orbits, eight eccentricities: the turn per lap\nis the disc's area, and only L sets it")
beta = np.linspace(0.01, 0.95, 100)
gam = 1 / np.sqrt(1 - beta ** 2)
axs[1].plot(beta, 2 * np.pi * (gam - 1), color=OURS, lw=4, alpha=0.5, label="Thomas turn per lap of a circular orbit")
axs[1].plot(beta, 2 * np.pi * (1 / np.sqrt(1 - beta ** 2) - 1), "k--", label="Sommerfeld advance per lap, same L")
axs[1].set_xlabel("orbital speed v/c"); axs[1].set_ylabel("radians per lap"); axs[1].legend(fontsize=9)
axs[1].set_title("B. They are the same curve: an orbit precesses exactly as\nmuch as a gyroscope on the circular orbit turns (R1)")
save(fig, "figO2-advance-is-the-disc-area.png")

# ------------------------------------------------------------ O3 coefficient
fig, axs = plt.subplots(1, 2, figsize=(13, 4.6))
names = ["scalar gravity\n(Nordström)", "special relativity\n(1/r force, Coulomb)", "general relativity\n(Schwarzschild)"]
co = res["coef"]
vals = [co["scalar (Nordström)"], co["SR (1/r force)"], co["GR (Schwarzschild)"]]
axs[0].bar(range(3), [-1, 1, 6], color=[GREY, OURS, KNOWN], alpha=0.6, label="first-order theory")
axs[0].plot(range(3), vals, "ko", label="simulated")
for i, v in enumerate(vals):
    axs[0].text(i + 0.12, v, f"{v:+.3f}", ha="left", va="center")
axs[0].axhline(0, color="black", lw=0.6)
axs[0].set_xticks(range(3)); axs[0].set_xticklabels(names)
axs[0].set_ylabel("advance ÷ (hodograph area / c²)"); axs[0].legend()
axs[0].set_ylim(-1.8, 7.2)
axs[0].set_title("A. One area, three coefficients: −1, +1, +6")
cs = np.array([5, 7, 10, 15, 25, 40, 80])
def run_gr(c, h=1.0, e=0.3):
    def rhs(phi, y): return [y[1], 1 / h ** 2 + 3 * y[0] ** 2 / c ** 2 - y[0]]
    s = integrate.solve_ivp(rhs, (0, 30 * math.pi), [(1 + e) / h ** 2, 0], rtol=1e-12, atol=1e-14, max_step=0.002)
    u, du = s.y
    idx = np.where((du[:-1] > 0) & (du[1:] <= 0))[0]
    ph = [s.t[i] + (-du[i] / (du[i + 1] - du[i])) * (s.t[i + 1] - s.t[i]) for i in idx]
    return np.diff(ph).mean() - 2 * math.pi
co_gr = [run_gr(c) / (math.pi / c ** 2) for c in cs]
axs[1].plot(1 / cs ** 2, co_gr, "o-", color=KNOWN)
axs[1].axhline(6, color="black", ls="--")
axs[1].set_xlabel("field strength (GM/hc)²"); axs[1].set_ylabel("GR coefficient")
axs[1].set_title("B. GR's 6 holds only in weak fields; strong fields\nbend away, so GR is not one exact area law")
save(fig, "figO3-one-area-three-coefficients.png")

# ------------------------------------------------------------ O4 fine structure
alpha = 7.2973525693e-3
mc2 = 510998.95
def E_d(n, k):
    d = k - math.sqrt(k * k - alpha ** 2)
    return mc2 / math.sqrt(1 + (alpha / (n - d)) ** 2) - mc2
fig, axs = plt.subplots(1, 2, figsize=(14, 5))
names_k = {(1, 1): "1S½", (2, 1): "2S½, 2P½", (2, 2): "2P3/2", (3, 1): "3S½, 3P½", (3, 2): "3P3/2, 3D3/2", (3, 3): "3D5/2"}
for n in (1, 2, 3):
    x0 = (n - 1) * 1.3
    axs[0].plot([x0, x0 + 0.9], [0, 0], color=GREY, lw=3)
    axs[0].text(x0 + 0.45, 4, f"n = {n} Bohr level: {n * n} state" + ("s" if n > 1 else ""), ha="center", fontsize=9, color=GREY)
    for k in range(1, n + 1):
        sh = (E_d(n, k) - (-alpha ** 2 * mc2 / 2 / n ** 2)) * 1e6
        axs[0].plot([x0, x0 + 0.9], [sh, sh], color=OURS, lw=2)
        axs[0].text(x0 + 0.92, sh, f"κ={k}: {names_k[(n, k)]}", fontsize=8, va="center")
axs[0].set_xlim(-0.1, 3.9); axs[0].set_ylim(-200, 14); axs[0].set_xticks([])
axs[0].set_ylabel("shift from the Bohr level (μeV)")
axs[0].set_title("A. The 3-sphere makes each Bohr level n²-fold (grey);\nthe area per κ splits it (blue: Dirac's fine structure)")
ks = np.arange(1, 6)
A = 2 * np.pi * (1 / np.sqrt(1 - (alpha / ks) ** 2) - 1)
axs[1].bar(ks, A * 1e5, color=OURS)
for kk, a in zip(ks, A):
    axs[1].text(kk, a * 1e5 + 0.4, f"{a * 1e5:.2f}", ha="center", fontsize=9)
axs[1].set_xlabel("κ = j + ½"); axs[1].set_ylabel("disc area (× 10⁻⁵ rad)")
axs[1].set_title("B. One number per κ: the area of the velocity disc of\nradius αc/κ. 2P3/2 − 2P1/2 = 10.950 GHz (measured 10.969)")
save(fig, "figO4-fine-structure-from-the-area.png")

# ------------------------------------------------------------ O5 Mercury
fig, axs = plt.subplots(1, 2, figsize=(13, 4.8))
GM, a, e = 1.32712440018e20, 5.7909050e10, 0.205630
h = math.sqrt(GM * a * (1 - e * e))
Rh = GM / h
t = np.linspace(0, 2 * np.pi, 300)
pr = a * (1 - e * e) / (1 + e * np.cos(t))
axs[0].plot(pr * np.cos(t) / 1e9, pr * np.sin(t) / 1e9, color=OURS)
axs[0].plot(0, 0, "*", color=KNOWN, ms=14)
axs[0].set_aspect("equal"); axs[0].set_xlabel("million km"); axs[0].set_title("A. Mercury's orbit")
axs[1].plot(Rh * np.cos(t) / 1e3, (Rh * np.sin(t) + e * Rh) / 1e3, color=OURS)
axs[1].fill(Rh * np.cos(t) / 1e3, (Rh * np.sin(t) + e * Rh) / 1e3, color=OURS, alpha=0.2)
axs[1].plot(0, 0, "k+")
axs[1].set_aspect("equal"); axs[1].set_xlabel("km/s")
axs[1].set_title(f"B. Its velocities: a circle of radius {Rh / 1e3:.1f} km/s.\n"
                 f"6 × area/c² per orbit = {res['mercury']['arcsec_per_century']:.2f}″ per century")
save(fig, "figO5-mercury-hodograph.png")
