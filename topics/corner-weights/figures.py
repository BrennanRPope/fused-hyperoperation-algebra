"""Figures for the corner-weights page. Run: python3 figures.py (numpy, matplotlib). Reads results.json."""
import json, itertools, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
R = json.load(open("results.json"))
plt.rcParams.update({"font.size": 10, "axes.titlesize": 11})
BLUE, ORANGE, GREEN, RED, PURPLE, GRAY = "#2a6fdb", "#e8833a", "#2e9e5b", "#d64545", "#8a5cc2", "#777777"

# ---------- fig1: the two choices ----------
fig, ax = plt.subplots(figsize=(11, 6.2)); ax.set_xlim(0, 3.3); ax.set_ylim(0, 2.45); ax.axis("off")
cols = ["Diagonal matrix\n(each wave alone)", "Cauchy matrix  1/(pᵢ+qⱼ)\n(one shared channel)", "Generic matrix\n(two or more channels)"]
rows = ["Evaluate at a\nnumber k", "Evaluate at a 2×2\nmatrix (trace k, det d)"]
cells = [["OUR OPERATION ∘ₖ\nsize weight k^(m−1)\nno pair factors (all A = 1)",
          "SIZE × PAIR\nk^(m−1) · ∏ pairs Aᵢⱼ\nKdV, KP line solitons\n(k only slides the waves)",
          "pair factors + 3-body terms\nnot a soliton solution\n(any 3-body term breaks KdV)"],
         ["THE PLANE LAW\nLucas weights Uₘ(k, d)\n(Fibonacci, tan, tanh ...)",
          "LUCAS × PAIR\n(0, 1): sine-Gordon kinks\n(0, −1): sinh-Gordon\nd = 0: back to KdV",
          "Lucas × pairs × 3-body\n(not explored)"]]
colors = [[BLUE, GREEN, GRAY], [PURPLE, GREEN, GRAY]]
for j, c in enumerate(cols):
    ax.text(0.75 + 0.85 * j + 0.42, 2.32, c, ha="center", va="center", fontsize=10.5, fontweight="bold")
for i, r in enumerate(rows):
    ax.text(0.35, 1.62 - 0.9 * i, r, ha="center", va="center", fontsize=10.5, fontweight="bold")
    for j in range(3):
        x0, y0 = 0.75 + 0.85 * j + 0.02, 1.2 - 0.9 * i
        ax.add_patch(FancyBboxPatch((x0, y0), 0.8, 0.82, boxstyle="round,pad=0.01", fc=colors[i][j], alpha=0.13, ec=colors[i][j], lw=1.5))
        ax.text(x0 + 0.4, y0 + 0.41, cells[i][j], ha="center", va="center", fontsize=9.3)
ax.annotate("", xy=(3.27, 2.15), xytext=(0.77, 2.15), arrowprops=dict(arrowstyle="->", color="k"))
ax.text(2.0, 2.19, "which MATRIX: its principal minors are the cube's corners (pair weights live here)", ha="center", fontsize=9, style="italic")
ax.annotate("", xy=(0.72, 0.25), xytext=(0.72, 2.05), arrowprops=dict(arrowstyle="->", color="k"))
ax.text(0.69, 0.18, "which DIAL t in det(I + tM)\n(size and Lucas weights live here)", ha="right", fontsize=9, style="italic")
ax.set_title("Fig 1. Every corner weighting so far is det(I + tM): two independent choices", fontsize=12, pad=4)
plt.savefig("fig1-two-choices.png", dpi=150, bbox_inches="tight"); plt.close()

# ---------- fig2: cube corners ----------
def cube(ax, labels, title, edge_labels=None):
    P = {(0,0,0):(0,0),(1,0,0):(2,0),(0,1,0):(0.9,0.75),(0,0,1):(0,2),(1,1,0):(2.9,0.75),(1,0,1):(2,2),(0,1,1):(0.9,2.75),(1,1,1):(2.9,2.75)}
    for a, b in itertools.combinations(P, 2):
        if sum(abs(x - y) for x, y in zip(a, b)) == 1:
            ax.plot(*zip(P[a], P[b]), color="#bbbbbb", lw=1.2, zorder=1)
    for v, (x, y) in P.items():
        m = sum(v); col = [GRAY, BLUE, GREEN, RED][m]
        ax.scatter([x], [y], s=60, color=col, zorder=3)
        ax.text(x + 0.08, y + 0.12, labels[v], fontsize=8.6, color=col, zorder=4)
    ax.set_xlim(-0.4, 4.6); ax.set_ylim(-0.5, 3.3); ax.axis("off"); ax.set_title(title, fontsize=10)
fig, axs = plt.subplots(1, 4, figsize=(15, 4.2))
nm = lambda v: "".join(s for s, b in zip("abc", v) if b) or "∅"
cube(axs[0], {v: (nm(v) if sum(v) else "1") + ("" if sum(v) < 2 else ("·k" if sum(v) == 2 else "·k²")) for v in itertools.product((0,1), repeat=3)},
     "Ours: weight by SIZE\nrank m → k^(m−1)")
pairs = {(1,1,0): "ab·A_ab", (1,0,1): "ac·A_ac", (0,1,1): "bc·A_bc", (1,1,1): "abc·A_ab·A_ac·A_bc"}
cube(axs[1], {v: pairs.get(v, nm(v) if sum(v) else "1") for v in itertools.product((0,1), repeat=3)},
     "Solitons: weight by PAIRS\ncorner = product over its pairs")
both = {(1,1,0): "ab·k·A_ab", (1,0,1): "ac·k·A_ac", (0,1,1): "bc·k·A_bc", (1,1,1): "abc·k²·A_abA_acA_bc"}
cube(axs[2], {v: both.get(v, nm(v) if sum(v) else "1") for v in itertools.product((0,1), repeat=3)},
     "BOTH = Cauchy matrix at dial k\n(k only slides the waves)")
luc = {0: "1", 1: "{n}", 2: "{n}: U₂=k | −d·U₁", 3: "abc: U₃=k²−d | −d·U₂"}
cube(axs[3], {v: luc[sum(v)].format(n=nm(v)) for v in itertools.product((0,1), repeat=3)},
     "Lucas: weight by SIZE on the plane\ntop weight | bottom weight")
for a in axs: pass
fig.suptitle("Fig 2. The 3-cube's corners under each weighting (grey = empty corner, blue = one wave, green = pairs, red = triple)", y=1.02)
plt.savefig("fig2-cube-corners.png", dpi=150, bbox_inches="tight"); plt.close()

# ---------- fig3: matrices and the three-body term ----------
p = np.array([0.6, 1.1, 1.9, 2.8]); x = np.array([0.9, 0.8, 1.1, 0.7])
Mdiag = np.diag(x)
Mcau = np.array([[2*np.sqrt(p[i]*p[j]*x[i]*x[j])/(p[i]+p[j]) for j in range(4)] for i in range(4)])
rng = np.random.default_rng(2)
r2, s2 = rng.uniform(0.5, 1.5, 4), rng.uniform(0.5, 1.5, 4)
Mgen = Mcau + 0.5 * np.outer(r2, s2) / (p[:, None] + p[None, :])
fig, axs = plt.subplots(1, 4, figsize=(16, 3.6), gridspec_kw={"width_ratios": [1, 1, 1, 1.4]}); fig.subplots_adjust(wspace=0.45)
for a, Mx, t in zip(axs[:3], [Mdiag, Mcau, Mgen], ["Diagonal: our operation\n(no coupling, all A = 1)", "Cauchy (KdV, symmetrised)\n2√(pᵢpⱼxᵢxⱼ)/(pᵢ+pⱼ)", "Cauchy + a second channel\n(rank-2 source)"]):
    im = a.imshow(Mx, cmap="Blues", vmin=0, vmax=1.2); a.set_title(t, fontsize=9.5)
    a.set_xticks(range(4)); a.set_yticks(range(4)); a.set_xticklabels(["wave 1","2","3","4"], fontsize=8); a.set_yticklabels(["wave 1","2","3","4"], fontsize=8)
    for i in range(4):
        for j in range(4): a.text(j, i, f"{Mx[i,j]:.2f}", ha="center", va="center", fontsize=7.5, color="k")
sr = R["source_rank"]
e = [v[0] for v in sr][1:]; tb = [abs(v[1]) for v in sr][1:]
axs[3].loglog(e, tb, "o-", color=RED); axs[3].axhline(abs(sr[0][1]) + 1e-16, color=GREEN, ls="--")
axs[3].text(0.012, 3e-14, "one channel (Cauchy): 3-body term = 0 (rounding)", color=GREEN, fontsize=8.5)
axs[3].set_xlabel("strength of the second channel"); axs[3].set_ylabel("|three-body term| of the triple corner")
axs[3].set_title("Pairs only ⇔ one shared channel", fontsize=10); axs[3].set_ylim(1e-15, 2)
fig.suptitle("Fig 3. The matrix decides the pair weights: off-diagonal coupling makes them, a second channel adds three-body terms", y=1.03)
plt.savefig("fig3-matrices-three-body.png", dpi=150, bbox_inches="tight"); plt.close()

# ---------- fig4: dressed waves and A = tanh^2 ----------
fig, axs = plt.subplots(1, 2, figsize=(12, 4.2))
th = np.linspace(0, 6, 300)
axs[0].plot(th, np.tanh(th/2)**2, color=GREEN, lw=2, label="A = tanh²(Δθ/2) = ((pᵢ−pⱼ)/(pᵢ+pⱼ))²")
axs[0].plot(th, 1/np.cosh(th/2)**2, color=ORANGE, lw=2, label="coupling 1 − A = sech²(Δθ/2)")
axs[0].axhline(1, color=BLUE, ls=":"); axs[0].text(4.2, 1.03, "A = 1: our operation", color=BLUE)
axs[0].text(0.15, 0.06, "A = 0 at equal speeds:\nexclusion", color=GREEN, fontsize=9)
axs[0].set_xlabel("rapidity difference Δθ = |log pᵢ − log pⱼ|"); axs[0].set_ylabel("value"); axs[0].legend(loc="center right", fontsize=8.5)
axs[0].set_title("KdV pair weight: velocity subtraction, squared")
sp_ = R["spread"]; rr = [s[0] for s in sp_]
axs[1].semilogx(rr, [s[1] for s in sp_], "o-", color=ORANGE, label="largest coupling 1 − A")
axs[1].semilogx(rr, [s[2] for s in sp_], "s-", color=PURPLE, label="largest gap |dressed eigenvalue − bare wave|")
axs[1].set_xlabel("speed ratio between neighbouring waves"); axs[1].legend(fontsize=8.5)
axs[1].set_title("Spread the speeds and KdV turns into our operation")
fig.suptitle("Fig 4. The pair weight is a distance between waves; far apart, the dressing vanishes", y=1.02)
plt.savefig("fig4-pair-weight-dressing.png", dpi=150, bbox_inches="tight"); plt.close()

# ---------- fig5: dial plane ----------
fig, ax = plt.subplots(figsize=(8.8, 6.6))
kk = np.linspace(-3.2, 3.2, 400)
ax.fill_between(kk, kk**2/4, 3.2, color=PURPLE, alpha=0.08); ax.fill_between(kk, -2.6, kk**2/4, color=ORANGE, alpha=0.08)
ax.plot(kk, kk**2/4, color="k", lw=1.3); ax.axhline(0, color=BLUE, lw=2.2)
ax.text(-3.0, 2.75, "ELLIPTIC  d > k²/4\nweights oscillate (Chebyshev)\nvalue = tan(sum of arctans)", color=PURPLE, fontsize=9)
ax.text(1.45, -2.45, "HYPERBOLIC  d < k²/4\nweights grow; value = ratio of\ntwo determinants det(I+αM)/det(I+βM)", color=ORANGE, fontsize=9)
ax.text(1.9, 0.12, "d = 0: our dial k\nweights k^(m−1)", color=BLUE, fontsize=9)
ax.text(-3.1, 1.95, "parabola k² = 4d:\nUₘ = m(k/2)^(m−1)\n(a derivative / a count)", fontsize=8.5)
pts = {(0,0):"plain addition",(1,0):"our ∘₁",(-1,0):"incl.–excl.",(0,1):"tan of a sum\n→ SINE-GORDON",(0,-1):"tanh / velocity\n→ SINH-GORDON",
       (1,-1):"Fibonacci",(2,-1):"Pell",(1,-2):"Jacobsthal",(3,2):"Mersenne",(2,1):"Uₘ = m",(1,1):"period 6"}
for (a, b), lab in pts.items():
    big = (a, b) in [(0,1),(0,-1)]
    ax.scatter([a], [b], s=70 if big else 35, color=GREEN if big else "k", zorder=3)
    dx_, ha_ = (-0.08, "right") if big else (0.07, "left")
    ax.text(a + dx_, b + 0.07, lab, fontsize=8.5 if not big else 9, color=GREEN if big else "k", fontweight="bold" if big else None, ha=ha_)
ax.set_xlim(-3.2, 3.6); ax.set_ylim(-2.6, 3.2); ax.set_xlabel("k  (trace of the 2×2 dial)"); ax.set_ylabel("d  (determinant of the 2×2 dial)")
ax.set_title("Fig 5. The Lucas plane: where you evaluate det(I + tM).\nOn the soliton (Cauchy) matrix, green points give known integrable equations")
plt.savefig("fig5-lucas-dial-plane.png", dpi=150, bbox_inches="tight"); plt.close()

# ---------- fig6: KdV tests ----------
fig, axs = plt.subplots(1, 2, figsize=(12, 4.2))
it = R["interp_s"]
axs[0].semilogy([a for a, b in it], [max(b, 1e-175) for a, b in it], "o-", color=RED)
axs[0].annotate("s = 1: Cauchy matrix\nexact KdV solution\n(residual 10⁻¹⁷¹ = zero)", (1, 1e-171), (0.45, 1e-120), arrowprops=dict(arrowstyle="->"), fontsize=9)
axs[0].set_xlabel("blend s: diagonal (0, ours) → Cauchy (1, KdV)"); axs[0].set_ylabel("KdV residual (relative)")
axs[0].set_title("Only the Cauchy end solves KdV")
pK = np.array([0.7, 1.3]); X = np.linspace(-14, 14, 1400)
def u_of(kgauge):
    eta = pK[:, None] * X[None, :]
    A12 = ((pK[0]-pK[1])/(pK[0]+pK[1]))**2
    tau = 1 + kgauge*np.exp(eta[0]) + kgauge*np.exp(eta[1]) + kgauge**2*A12*np.exp(eta[0]+eta[1])
    L = np.log(tau); return 2*np.gradient(np.gradient(L, X), X)
for kg, c in [(0.05, BLUE), (1, GREEN), (20, ORANGE)]:
    axs[1].plot(X, u_of(kg), color=c, label=f"k = {kg}")
axs[1].set_xlabel("x"); axs[1].set_ylabel("u = 2 (log τ)''"); axs[1].legend(); axs[1].set_xlim(-12, 12)
axs[1].set_title("Size dial k: each wave slides by −log(k)/pᵢ\n(same two solitons, met at a different moment)")
fig.suptitle("Fig 6. KdV: pairs must come from the Cauchy matrix; the size dial is a gauge", y=1.02)
plt.savefig("fig6-kdv-tests.png", dpi=150, bbox_inches="tight"); plt.close()

# ---------- fig7: KP field and BZ ----------
fig, axs = plt.subplots(1, 2, figsize=(12.5, 4.6), gridspec_kw={"width_ratios": [1.4, 1]})
pk = np.array([-0.35, 0.45, 1.3]); qk = np.array([0.95, 1.05, 1.15])
xs = np.linspace(-20, 20, 480); ys = np.linspace(-20, 20, 480); XX, YY = np.meshgrid(xs, ys)
def tauKP(XX, YY, t=0.0):
    E = [np.exp((pk[i]+qk[i])*XX + (pk[i]**2-qk[i]**2)*YY + (pk[i]**3+qk[i]**3)*t) for i in range(3)]
    tot = np.zeros_like(XX)
    for r in range(4):
        for S in itertools.combinations(range(3), r):
            w = math.prod((pk[i]-pk[j])*(qk[i]-qk[j])/((pk[i]+qk[j])*(pk[j]+qk[i])) for i, j in itertools.combinations(S, 2))
            term = np.ones_like(XX)
            for i in S: term = term*E[i]
            tot = tot + w*term
    return tot
L = np.log(tauKP(XX, YY)); dx = xs[1]-xs[0]
U = 2*np.gradient(np.gradient(L, dx, axis=1), dx, axis=1)
im = axs[0].imshow(U, extent=[xs[0], xs[-1], ys[0], ys[-1]], origin="lower", cmap="viridis")
plt.colorbar(im, ax=axs[0], label="u = 2 ∂²ₓ log det(I + M)")
axs[0].set_xlabel("x"); axs[0].set_ylabel("y"); axs[0].set_title("KP: three line solitons in three directions\nfrom one Cauchy matrix (regular: every A > 0)")
bz = R["bz"]
axs[1].bar(["one channel\n(static, Weyl)", "two channels\n(rotating, Kerr-type)"], [abs(bz["threebody_one_channel"]) + 1e-16, abs(bz["threebody_two_channel"])], color=[GREEN, RED])
axs[1].set_yscale("log"); axs[1].set_ylabel("|three-body term| for three solitons")
axs[1].set_title("Belinski–Zakharov gravitational solitons:\nsame rule (Cauchy 1/(ρ²+μₖμₗ))")
fig.suptitle("Fig 7. The same matrix picture in 2-D water waves (KP) and in gravity (Belinski–Zakharov)", y=1.02)
plt.savefig("fig7-kp-gravity.png", dpi=150, bbox_inches="tight"); plt.close()

# ---------- fig8: k against A ----------
fig, ax = plt.subplots(figsize=(8.6, 6))
ax.set_xlim(-0.1, 3.1); ax.set_ylim(-0.15, 2.7)
ax.axhspan(0, 1, color=GREEN, alpha=0.08); ax.axhline(1, color=BLUE, lw=2.2); ax.axhline(0, color=GRAY, lw=2.2)
kk = np.linspace(0.05, 3.1, 200)
for c in (0.25, 0.5, 1.0, 2.0):
    ax.plot(kk, c/kk, color=GRAY, ls=":", lw=1); ax.text(3.0, c/3.0 + 0.03, f"kA = {c}", fontsize=7.5, color=GRAY, ha="right")
for y in (0.5, 1.6, 2.3):
    ax.annotate("", xy=(2.9, y), xytext=(0.2, y), arrowprops=dict(arrowstyle="<->", color=PURPLE, alpha=0.6))
ax.text(1.0, 2.36, "k: gauge for solitons (slides each wave by −log k / p)", color=PURPLE, fontsize=9)
ax.text(0.05, 1.04, "A = 1: our operation ∘ₖ (no pairs)", color=BLUE, fontsize=9.5, fontweight="bold")
ax.text(0.05, 0.04, "A = 0: hard core = plain addition", color=GRAY, fontsize=9.5, fontweight="bold")
ax.text(1.6, 0.45, "KdV and most KP pairs:\n0 ≤ A < 1 (repulsive)", color=GREEN, fontsize=9.5)
ax.axhspan(1, 2.7, xmin=0, xmax=1, color=ORANGE, alpha=0.06)
ax.text(1.6, 1.35, "three KP waves sharing ONE A:\nfound only for A > 1\n(numerical; none below 1)", color=ORANGE, fontsize=9.5)
ax.set_xlabel("size dial k"); ax.set_ylabel("pair dial A (one shared pair weight)")
ax.set_title("Fig 8. The two dials, k against A. Dotted lines: same two-input law a + b + kA·ab;\nthey split at three inputs (k²A³ against k²A²)")
plt.savefig("fig8-k-against-A.png", dpi=150, bbox_inches="tight"); plt.close()
print("8 figures written")

# ---------- fig9: the face picture as vertex operators ----------
fig = plt.figure(figsize=(14, 4.6))
ax = fig.add_subplot(1, 3, 1); ax.set_xlim(0, 10); ax.set_ylim(0, 6); ax.axis("off")
ax.set_title("(a) each wave = a MULTIPLY part × a SHIFT part", fontsize=10.5)
def box(x, y, w, txt, c):
    ax.add_patch(FancyBboxPatch((x, y), w, 0.9, boxstyle="round,pad=0.05", fc=c, alpha=0.18, ec=c, lw=1.5))
    ax.text(x + w / 2, y + 0.45, txt, ha="center", va="center", fontsize=9.5)
ax.text(0.2, 5.3, "wave 1 = X(p)", fontsize=10, fontweight="bold"); ax.text(5.3, 5.3, "wave 2 = X(p′)", fontsize=10, fontweight="bold")
box(0.2, 4.1, 2.2, "multiply\ne^(2p·t)", BLUE); box(2.6, 4.1, 2.2, "shift\nt → t − 1/p", ORANGE)
box(5.3, 4.1, 2.2, "multiply\ne^(2p′·t)", BLUE); box(7.7, 4.1, 2.2, "shift\nt → t − 1/p′", ORANGE)
ax.annotate("", xy=(6.4, 3.9), xytext=(3.7, 3.9), arrowprops=dict(arrowstyle="<->", color=RED, lw=2, connectionstyle="arc3,rad=0.4"))
ax.text(5.05, 2.55, "slide wave 1's shift past\nwave 2's multiply: one factor", ha="center", fontsize=9.5, color=RED)
ax.text(5.0, 1.35, "cost = exp(−4 artanh(p′/p)) = ((p − p′)/(p + p′))²  =  Aᵢⱼ", ha="center", fontsize=10.5, fontweight="bold")
ax.text(5.0, 0.45, "every PAIR pays once; a free field has nothing else to pay\n(t stands for all the times t₁, t₃, t₅, ...; shown schematically)", ha="center", fontsize=8.5, style="italic", color=GRAY)
ax = fig.add_subplot(1, 3, 2)
x = np.linspace(0.0, 0.995, 400)
ax.plot(x, ((1 - x) / (1 + x)) ** 2, color="k", lw=2.5, label="exact pair weight ((p−p′)/(p+p′))²")
for J, c in ((1, ORANGE), (3, PURPLE), (9, BLUE)):
    part = sum(4 * x ** j / j for j in range(1, J + 1, 2))
    ax.plot(x, np.exp(-part), color=c, lw=1.3, ls="--", label=f"keeping modes up to j = {J}")
g = np.array(R["vo_exclusion"]); ax.plot(1 - g[:, 0], g[:, 1], "o", color=RED, ms=7, label="computed from the operators")
ax.annotate("equal speeds: weight 0\n(two identical waves can't both be on)", xy=(0.99, 0.002), xytext=(0.30, 0.32), fontsize=9, arrowprops=dict(arrowstyle="->", color=RED))
ax.set_xlabel("speed ratio p′/p"); ax.set_ylabel("pair weight Aᵢⱼ"); ax.set_ylim(-0.02, 1.02)
ax.set_title("(b) the pair weight is built up mode by mode", fontsize=10.5); ax.legend(fontsize=8, loc="upper right")
ax = fig.add_subplot(1, 3, 3)
sw = np.array(R["vo_beta_sweep"])
ax.plot(sw[:, 0], sw[:, 1], "o-", color=RED, lw=2, label="interacting field (squared term in the shift)")
ax.axhline(0, color=GREEN, lw=2.5, label="free field: 3-body term exactly 0")
ax.set_xlabel("strength β of the squared term"); ax.set_ylabel("3-body term at the all-three corner")
ax.set_title("(c) free field: pairs only. Add an interaction: 3-body", fontsize=10.5); ax.legend(fontsize=8.5, loc="upper left")
ax.text(0.40, 0.05, "Gravity (static Weyl solutions) is the\nfree-field case: the potential adds, and\nthe other metric function is built from\npairs only (3-body part 0).", fontsize=8.5, color=GRAY, transform=ax.transAxes, ha="left", va="bottom")
fig.suptitle("Fig 9. Face picture: waves as vertex operators of a free field. The pair weights are the cost of putting them in order (DRAFT)", fontsize=12)
plt.tight_layout(); plt.savefig("fig9-vertex-operators.png", dpi=150, bbox_inches="tight"); plt.close()
