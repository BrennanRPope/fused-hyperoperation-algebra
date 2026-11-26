"""Pictures for the new-law pass (run after new_laws.py)."""
import json, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mpmath as mp

N = json.load(open("new_laws.json"))
BLUE, ORANGE, GREY, RED, GREEN = "#2b6cb0", "#dd6b20", "#888888", "#c53030", "#2f855a"
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})
def save(fig, name):
    fig.savefig(f"figures/{name}.png", dpi=130, bbox_inches="tight"); plt.close(fig); print("wrote", name)
def pball(p, m=400):
    t = np.linspace(0, 2*np.pi, m); D = np.column_stack([np.cos(t), np.sin(t)])
    g = np.max(np.abs(D), 1) if np.isinf(p) else np.sum(np.abs(D)**p, 1)**(1/p)
    return D/g[:, None]
def ngon(n, rot=np.pi/2):
    a = rot + 2*np.pi*np.arange(n + 1)/n; return np.column_stack([np.cos(a), np.sin(a)])

# ---------------------------------------------------------------- fig13: Bell floor and saddle
fig = plt.figure(figsize=(17, 4.9)); gs = fig.add_gridspec(1, 3, width_ratios=[1.25, 1, 1.15], wspace=0.28)
ax0 = fig.add_subplot(gs[0])
h = N["bell_floor"]["hist"]; edges = np.linspace(2.7, 4.05, len(h) + 1)
ax0.bar(edges[:-1], h, width=np.diff(edges), align="edge", color=BLUE, alpha=.8, edgecolor="white")
ax0.axvline(2*np.sqrt(2), color=ORANGE, lw=2.5); ax0.text(2*np.sqrt(2) + .03, max(h)*.92, "2√2", color=ORANGE, ha="left", fontsize=12)
ax0.axvspan(2.7, 2*np.sqrt(2), color=GREY, alpha=.15)
ax0.annotate("grey: no symmetric shape lands here", (2.76, max(h)*.4), (3.0, max(h)*.62), color="#555555", fontsize=10, arrowprops=dict(arrowstyle="->", color="#555555"))
ax0.set_xlabel("best Bell score (CHSH)"); ax0.set_ylabel("shapes"); ax0.set_xlim(2.7, 4.05)
ax0.set_title(f"{N['bell_floor']['n_shapes']} random symmetric shapes (polygons and smooth)\nlowest {N['bell_floor']['min']:.4f}: Tsirelson's value is a floor")
ax1 = fig.add_subplot(gs[1])
ps = np.concatenate([np.linspace(1, 2, 60), np.linspace(2, 12, 120)])
ax1.plot(np.log2(ps), [2**(1 + max(1/p, 1 - 1/p)) for p in ps], color=BLUE, lw=2.5, label="roundness dial: 2^(1 + max(1/p, 1 − 1/p))")
for p, mk in [(1.25, "o"), (1.5, "o"), (3, "o"), (5, "o")]:
    ax1.plot(np.log2(p), N["bell_3d"][str(float(p))], "s", color=GREEN, ms=8, mfc="none", mew=2)
ax1.plot([], [], "s", color=GREEN, mfc="none", mew=2, label="same shapes in 3D (polished)")
ax1.axhline(2*np.sqrt(2), color=ORANGE, ls="--", lw=1.2)
ax1.set_xticks(np.log2([1, 1.5, 2, 3, 4, 8])); ax1.set_xticklabels(["1\ncross-polytope", "1.5", "2\nsphere", "3", "4", "8 → cube"])
ax1.set_ylabel("CHSH"); ax1.legend(fontsize=9, loc="upper center"); ax1.set_title("Smooth → angular, staying symmetric:\nthe score only goes up from 2√2")
ax2 = fig.add_subplot(gs[2]); ax2.set_aspect("equal"); ax2.axis("off")
# saddle sketch: centre disk, arrows to octagon/square (up) and to pentagon/triangle (down)
def draw(P, cx, cy, sc, c):
    ax2.fill(P[:, 0]*sc + cx, P[:, 1]*sc + cy, color=c, alpha=.25); ax2.plot(P[:, 0]*sc + cx, P[:, 1]*sc + cy, color=c, lw=1.6)
draw(pball(2), 0, 0, .55, ORANGE); ax2.text(0.62, -0.05, "disk\n2.828", ha="left", va="center", fontsize=10)
draw(ngon(8, np.pi/8), -1.7, 1.35, .45, BLUE); ax2.text(-1.7, 0.72, "octagon 2.828", ha="center", fontsize=9)
draw(ngon(4, np.pi/4), 1.7, 1.35, .45, BLUE); ax2.text(1.7, 0.72, "square 4", ha="center", fontsize=9)
draw(ngon(5), -1.7, -1.55, .45, RED); ax2.text(-1.7, -2.2, "pentagon 2.683*", ha="center", fontsize=9)
draw(ngon(3), 1.7, -1.55, .45, RED); ax2.text(1.7, -2.2, "triangle 2*", ha="center", fontsize=9)
for (x, y, c) in [(-1.2, .95, BLUE), (1.2, .95, BLUE), (-1.2, -1.05, RED), (1.2, -1.05, RED)]:
    ax2.annotate("", (x, y), (0.45*np.sign(x), 0.45*np.sign(y)), arrowprops=dict(arrowstyle="->", color=c, lw=1.6))
ax2.text(0, 1.95, "symmetric and angular: up (or flat)", ha="center", color=BLUE, fontsize=10)
ax2.text(0, -2.75, "toward the triangle (shape dial): down\n*Janotta et al. 2011", ha="center", color=RED, fontsize=9)
ax2.set_xlim(-2.6, 2.6); ax2.set_ylim(-3.1, 2.3)
ax2.set_title("Where the two dials cross, 2√2 is a saddle")
fig.suptitle("N1 · Bell's bound is a floor on the roundness dial and a ceiling toward the triangle", fontsize=14, y=1.04)
save(fig, "fig13-new-bell-floor")

# ---------------------------------------------------------------- fig14: Fermat for every rational p
fig, ax = plt.subplots(1, 2, figsize=(16, 4.4), gridspec_kw={"width_ratios": [1.7, 1]})
from fractions import Fraction
fr = sorted({Fraction(m, n) for n in range(1, 9) for m in range(n, 8*n + 1) if math.gcd(m, n) == 1 and m/n <= 6})
for f in fr:
    rich = f.numerator <= 2
    ax[0].plot(float(f), 0 if rich else 0, "o", color=BLUE if rich else RED, ms=11 if rich else 3.5, alpha=1 if rich else .55, zorder=3 if rich else 2)
ax[0].plot(6.25, 0, "o", color=BLUE, ms=11); ax[0].text(6.25, .16, "∞", ha="center", fontsize=13, color=BLUE)
for x, t in [(1, "p = 1\ntaxicab"), (2, "p = 2\nPythagoras")]:
    ax[0].text(x, .16, t, ha="center", color=BLUE, fontsize=10)
for f, dx in [(Fraction(3, 2), .12), (Fraction(4, 3), -.12), (Fraction(5, 2), 0), (Fraction(3), 0), (Fraction(4), 0)]:
    ax[0].annotate(f"{f}", (float(f), 0), (float(f) + dx, -.22), ha="center", fontsize=9, color=RED, arrowprops=dict(arrowstyle="-", color=RED, lw=.6))
ax[0].set_ylim(-.4, .4); ax[0].set_yticks([]); ax[0].set_xlim(.8, 6.5); ax[0].spines["left"].set_visible(False)
ax[0].set_xlabel("roundness p (every fraction with denominator up to 8 is drawn)")
ax[0].set_title("Blue: infinitely many rational points.  Red: only the 4 trivial ones.\np = m/n behaves like Fermat's exponent m (the numerator)")
for p, c, lab in [(1.5, RED, "p = 3/2 ↔ u³ + v³ = 1 (x = u²)"), (4/3, GREEN, "p = 4/3 ↔ u⁴ + v⁴ = 1 (x = u³)")]:
    B = pball(p); ax[1].plot(B[:, 0], B[:, 1], color=c, lw=2, label=lab)
B = pball(2); ax[1].plot(B[:, 0], B[:, 1], color=ORANGE, lw=1, alpha=.5)
ax[1].plot([1, 0, -1, 0], [0, 1, 0, -1], "o", color="k", ms=6)
ax[1].set_aspect("equal"); ax[1].axis("off"); ax[1].legend(fontsize=9, loc="lower center", bbox_to_anchor=(.5, -.28))
ax[1].set_title("In-between curves keep only\nthe 4 corner points")
fig.suptitle("N2 · Fermat for every rational roundness: only taxicab, round and square circles have rational points to spare", fontsize=14, y=1.05)
save(fig, "fig14-new-fermat-rational-p")

# ---------------------------------------------------------------- fig15: Lee-Yang rule a <= p
fig, ax = plt.subplots(1, 3, figsize=(18, 5), gridspec_kw={"width_ratios": [1.15, 1, 1]})
pp = np.linspace(1, 6, 400)
ax[0].fill_between(pp, 0, pp, color=BLUE, alpha=.18); ax[0].fill_between(pp, pp, 8, color=RED, alpha=.12)
ax[0].plot(pp, pp, color="k", lw=2)
for pe in [2, 4, 6]: ax[0].plot([pe, pe], [0, 8], color=BLUE, lw=3)
far = N["lee_yang_far"]
for k_, v in far.items():
    p = float(k_.split(",")[0][2:]); a = float(k_.split("a=")[1])
    if v.get("window") is None: ax[0].plot(p, a, "x", color=GREY, ms=7); continue
    on = v["changes"] > 10
    ax[0].plot(p, a, "o", color=BLUE if on else RED, ms=7, mec="k", mew=.6)
ax[0].text(4.6, 1.2, "zeros on the line", color=BLUE, fontsize=11); ax[0].text(1.2, 6.5, "zeros leave the line", color=RED, fontsize=11)
ax[0].text(2.05, 7.3, "even p:\nalways on", color=BLUE, fontsize=9)
ax[0].text(5.1, 5.55, "a = p", fontsize=11, rotation=33)
ax[0].plot([], [], "x", color=GREY, label="crossover too far out to test")
ax[0].set_xlabel("roundness p"); ax[0].set_ylabel("a = (n − 1)/p"); ax[0].set_ylim(0, 8); ax[0].legend(fontsize=8, loc="upper right")
ax[0].set_title("The rule: Lee–Yang holds iff a ≤ p\n(or p is an even whole number)")
# panel 2: G(s) far out just below and above the line, p = 1.5
mp.mp.dps = 80
def Gs(p, a, ss):
    c = [(-1)**k/mp.factorial(2*k)*mp.beta(mp.mpf(2*k + 1)/p, a + 1)/p for k in range(170)]
    out = []
    for s in ss:
        s2 = mp.mpf(s)**2; acc = mp.mpf(0)
        for ck in reversed(c): acc = acc*s2 + ck
        out.append(float(acc*mp.mpf(s)**(a + 1)))
    out = np.array(out); return out/np.max(np.abs(out))
ss = np.linspace(20, 80, 500)
for a, c, lab in [(1.0, BLUE, "a = 1.0 < p: keeps crossing zero"), (2.5, RED, "a = 2.5 > p: stops crossing")]:
    ax[1].plot(ss, Gs(1.5, a, ss), color=c, lw=1.6, label=lab)
ax[1].axhline(0, color=GREY, lw=.8); ax[1].legend(fontsize=9, loc="upper right")
ax[1].set_xlabel("s (imaginary field)"); ax[1].set_ylabel("G(s) · s^(a+1), scaled to fit"); ax[1].set_ylim(-1.2, 1.6); ax[1].set_title("p = 1.5, exact sums: the rim's wave\nagainst the centre's steady push")
# panel 3: in spin terms, largest dimension n for which a uniform p-ball spin keeps Lee-Yang
pq = np.linspace(1, 3.2, 300); nmax = np.floor(1 + pq**2)
ax[2].step(pq, nmax, where="post", color=BLUE, lw=2.2)
for p, lab, xy in [(1, "octahedron (p = 1): n ≤ 2,\nso the 3D octahedron fails", (1.05, 0.6)), (1.5, "p = 1.5: n ≤ 3", (1.0, 5.2)), (3, "p = 3: n ≤ 10", (2.45, 11.6))]:
    ax[2].plot(p, math.floor(1 + p*p), "o", color=ORANGE, ms=8); ax[2].annotate(lab, (p, math.floor(1 + p*p)), xy, fontsize=9, arrowprops=dict(arrowstyle="-", color=GREY, lw=.6))
ax[2].axvline(2, color=GREY, ls=":"); ax[2].text(1.6, 8.6, "p = 2 itself:\nevery n (Bessel)", fontsize=9, color="#555555")
ax[2].set_xlabel("roundness p"); ax[2].set_ylabel("largest spin dimension n"); ax[2].set_ylim(0, 13)
ax[2].set_title("A spin uniform in an n-dimensional p-ball:\nLee–Yang iff n ≤ 1 + p²")
fig.suptitle("N3 · Where Lee–Yang breaks on the roundness dial: one line, a = p", fontsize=14, y=1.04)
save(fig, "fig15-new-lee-yang-rule")
