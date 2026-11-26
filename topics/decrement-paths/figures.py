"""Pictures for decrement-paths/README.md. Self-contained (numpy + matplotlib)."""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "figures"); os.makedirs(OUT, exist_ok=True)
C1, C2, C3, C4, GR = "#257de0", "#f47230", "#00b894", "#f2a600", "#868e96"
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
gx, gw = np.polynomial.legendre.leggauss(200); TT = 0.5*(gx + 1); WW = 0.5*gw

def save(fig, name):
    fig.tight_layout(); fig.savefig(os.path.join(OUT, name), dpi=140); plt.close(fig); print("wrote", name)

def split(P, dP, grad):
    return (grad(P(TT))*dP(TT)*WW).sum(axis=1)
union_grad = lambda p: np.prod(1 - p, axis=0)/(1 - p)

# ---------------------------------------------------------------- Fig 1: paths through the square
q = np.array([0.30, 0.50]); tot = 1 - np.prod(1 - q); mu = -np.log(1 - q)
paths = {
    "uniform timing (UDD): straight in the cube": (lambda s: q[:, None]*s, lambda s: q[:, None]*np.ones_like(s), C1, "-"),
    "constant force: straight in log coordinates": (lambda s: 1 - np.exp(-mu[:, None]*s), lambda s: mu[:, None]*np.exp(-mu[:, None]*s), C2, "-"),
    "cause 1 early, cause 2 late": (lambda s: q[:, None]*np.vstack([1-(1-s)**3, s**3]), lambda s: q[:, None]*np.vstack([3*(1-s)**2, 3*s**2]), C3, "--"),
    "cause 2 early, cause 1 late": (lambda s: q[:, None]*np.vstack([s**3, 1-(1-s)**3]), lambda s: q[:, None]*np.vstack([3*s**2, 3*(1-s)**2]), C4, "--"),
}
fig, axs = plt.subplots(1, 2, figsize=(11.5, 5.2), gridspec_kw={"width_ratios": [1.1, 1]})
ax = axs[0]
A, B = np.meshgrid(np.linspace(0, 0.32, 200), np.linspace(0, 0.52, 200))
cs = ax.contour(A, B, 1 - (1-A)*(1-B), levels=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6], colors=GR, linewidths=0.6)
ax.clabel(cs, fmt="total %.1f", fontsize=7)
s = np.linspace(0, 1, 200); rows = []
for lab, (P, dP, c, ls) in paths.items():
    pts = P(s); ax.plot(pts[0], pts[1], color=c, ls=ls, lw=2, label=lab)
    rows.append((lab, split(P, dP, union_grad), c))
ax.plot([0], [0], "o", color="k", ms=6); ax.plot([q[0]], [q[1]], "o", color="k", ms=6)
ax.annotate("start of year\n(nothing has struck)", (0, 0), (0.06, 0.035), fontsize=8, arrowprops=dict(arrowstyle="->", color="k", lw=0.8))
ax.annotate("end of year\n(q₁, q₂) = (0.30, 0.50)", q, (0.15, 0.45), fontsize=8, ha="center")
ax.set_xlabel("chance cause 1 has struck by time t"); ax.set_ylabel("chance cause 2 has struck by time t")
ax.set_title("(a) A timing assumption is a path through the square.\nGrey lines: the total chance the policy has ended.", fontsize=9)
ax.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=2, frameon=False)
ax = axs[1]
y = np.arange(len(rows))[::-1]
for yi, (lab, sp_, c) in zip(y, rows):
    ax.barh(yi, sp_[0], color=C1, height=0.55, edgecolor="white", lw=2)
    ax.barh(yi, sp_[1], left=sp_[0], color=C2, height=0.55, edgecolor="white", lw=2)
    ax.text(sp_[0]/2, yi, f"{sp_[0]:.3f}", ha="center", va="center", color="white", fontsize=8)
    ax.text(sp_[0] + sp_[1]/2, yi, f"{sp_[1]:.3f}", ha="center", va="center", color="white", fontsize=8)
ax.axvline(tot, color="k", lw=1); ax.text(tot, len(rows) - 0.45, f" total {tot:.2f}\n on every path", fontsize=8, va="top")
ax.set_yticks(y); ax.set_yticklabels([r[0].replace(": ", ":\n") for r in rows], fontsize=7.5)
ax.set_xlim(0, 0.8); ax.set_xlabel("dependent rate: policies ending by each cause")
ax.bar(0, 0, color=C1, label="cause 1"); ax.bar(0, 0, color=C2, label="cause 2")
ax.legend(fontsize=7.5, loc="lower right", frameon=False)
ax.set_title("(b) Same end point, same total. The split depends on the path.", fontsize=9)
save(fig, "fig1-paths.png")

# ---------------------------------------------------------------- Fig 2: the equal corner split
q3 = np.array([0.30, 0.50, 0.15]); names = ["cause 1", "cause 2", "cause 3"]
fig, axs = plt.subplots(1, 2, figsize=(11.5, 4.3), gridspec_kw={"width_ratios": [1, 1.2]})
ax = axs[0]; ax.axis("off")
# a drawn cube with its corners named
o = np.array([0, 0]); ex, ey, ez = np.array([1, 0]), np.array([0.45, 0.35]), np.array([0, 1])
V = {(a, b, c): o + a*ex + b*ey + c*ez for a in (0, 1) for b in (0, 1) for c in (0, 1)}
for u in V:
    for i in range(3):
        if u[i] == 0:
            w = list(u); w[i] = 1; w = tuple(w)
            ax.plot(*zip(V[u], V[w]), color=GR, lw=1)
lab = {(0,0,0): "∅\n0", (1,0,0): "{1}\nq₁", (0,1,0): "{2}\nq₂", (0,0,1): "{3}\nq₃",
       (1,1,0): "{1,2}\n−q₁q₂", (1,0,1): "{1,3}\n−q₁q₃", (0,1,1): "{2,3}\n−q₂q₃", (1,1,1): "{1,2,3}\n+q₁q₂q₃"}
for u, txt in lab.items():
    sz = sum(u); col = [GR, C1, C2, C3][sz]
    ax.plot(*V[u], "o", color=col, ms=9)
    ax.annotate(txt, V[u], (V[u][0] + (0.07 if u[0] else -0.07), V[u][1] + 0.06), fontsize=8,
                ha="left" if u[0] else "right")
ax.set_xlim(-0.5, 1.9); ax.set_ylim(-0.25, 1.6); ax.set_aspect("equal")
ax.set_title("(a) Each corner of the cube is a group of causes that would collide.\n"
             "Its worth (the Harsanyi dividend) is ±∏ q over the group.", fontsize=9)
ax = axs[1]
parts = []
for j in range(3):
    o_ = [i for i in range(3) if i != j]
    own = q3[j]; pair = -q3[j]*q3[o_[0]]/2 - q3[j]*q3[o_[1]]/2; trip = q3.prod()/3
    parts.append((own, pair, trip))
x = np.arange(3)
for j, (own, pair, trip) in enumerate(parts):
    ax.bar(j, own, color=C1, width=0.55, edgecolor="white", lw=2, label="own corner, whole" if j == 0 else None)
    ax.bar(j, pair, bottom=own, color="white", hatch="///", width=0.55, edgecolor=C2, lw=1.5, label="pair corners, halved: taken away" if j == 0 else None)
    ax.bar(j, trip, bottom=own + pair, color=C3, width=0.55, edgecolor="white", lw=2, label="triple corner, in thirds" if j == 0 else None)
    tot_j = own + pair + trip
    ax.plot([j - 0.33, j + 0.33], [tot_j, tot_j], color="k", lw=1.5)
    ax.text(j + 0.36, tot_j, f"{tot_j:.4f}\n= UDD\n= Shapley", fontsize=7.5, va="center")
ax.set_xticks(x); ax.set_xticklabels(names); ax.set_xlim(-0.5, 2.9)
ax.set_ylabel("share of the total"); ax.legend(fontsize=7.5, frameon=False, loc="upper right")
ax.set_title("(b) Uniform timing hands each corner to its members in equal shares.\n"
             "q = (0.30, 0.50, 0.15). Exact for 2 to 5 causes (checks A).", fontsize=9)
save(fig, "fig2-corner-split.png")

# ---------------------------------------------------------------- Fig 3: every line dial
xs = np.array([0.30, 0.50]); ks = np.linspace(-1.9, 1.9, 77)
sh, pr = [], []
for kv in ks:
    totk = ((1 + kv*xs[0])*(1 + kv*xs[1]) - 1)/kv if abs(kv) > 1e-9 else xs.sum()
    sh.append(xs[0] + kv*xs[0]*xs[1]/2)
    L = np.log1p(kv*xs)/kv if abs(kv) > 1e-9 else xs
    pr.append(L[0]/L.sum()*totk)
fig, ax = plt.subplots(figsize=(7.5, 4.2))
ax.plot(ks, sh, color=C1, lw=2, label="straight in the cube = Shapley (UDD at k = −1)")
ax.plot(ks, pr, color=C2, lw=2, label="straight in log coordinates = proportional (constant force at k = −1)")
ax.axvline(-1, color=GR, lw=0.8, ls=":"); ax.text(-0.98, 0.17, "k = −1: actuarial\n(union rule)", fontsize=7.5)
ax.axvline(0, color=GR, lw=0.8, ls=":"); ax.text(0.02, 0.17, "k = 0: plain sum,\nrules agree", fontsize=7.5)
ax.set_xlabel("dial k  (game v(S) = (∏(1 + k x) − 1)/k)"); ax.set_ylabel("share of player 1  (x₁ = 0.30, x₂ = 0.50)")
ax.legend(fontsize=7.5, frameon=False, loc="lower right")
ax.set_title("Both textbook rules exist on every line dial: one diagonal rule, read in two coordinates", fontsize=9)
save(fig, "fig3-every-dial.png")

# ---------------------------------------------------------------- Fig 4: dependent causes
def plane_grad(p, kv, dv):
    a, b = p; den = 1 - dv*a*b
    return np.array([(1 + kv*b)/den + dv*b*(a + b + kv*a*b)/den**2, (1 + kv*a)/den + dv*a*(a + b + kv*a*b)/den**2])
ths = np.linspace(0, 0.95, 40); gap = []; first = []
for th in ths:
    kv, dv = -(1 + th), th
    st = split(lambda s: xs[:, None]*s, lambda s: xs[:, None]*np.ones_like(s), lambda p: plane_grad(p, kv, dv))
    F = (xs.sum() + kv*xs.prod())/(1 - dv*xs.prod())
    gap.append((F + xs[0] - xs[1])/2 - st[0]); first.append(-dv*xs.prod()*(xs[0] - xs[1])/6)
fig, ax = plt.subplots(figsize=(7.5, 4.2))
ax.plot(ths, gap, color=C1, lw=2, label="exact gap (path integral)")
ax.plot(ths, first, color=C2, lw=2, ls="--", label="first-order formula −θ·x·y·(x − y)/6")
ax.errorbar([0.5], [0.210811 - 0.20849], yerr=[3*np.sqrt(0.20849/2e6)], fmt="o", color="k", ms=6, capsize=3,
            label="simulation, 2 million lives (θ = 0.5, ±3 s.e.)")
ax.axhline(0, color=GR, lw=0.8)
ax.set_xlabel("dependence θ between the two causes (Ali–Mikhail–Haq copula; θ = 0 independent)")
ax.set_ylabel("Shapley share − uniform-timing share, cause 1")
ax.legend(fontsize=7.5, frameon=False, loc="upper left")
ax.set_title("With dependent causes the total is a plane law, and uniform timing\nstops giving the Shapley split "
             "(q₁ = 0.30, q₂ = 0.50)", fontsize=9)
save(fig, "fig4-dependent.png")
