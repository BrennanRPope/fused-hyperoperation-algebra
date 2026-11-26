"""DRAFT · Pictures for tunnelling/README.md. Run: python3 figures.py  (writes figures/*.png)"""
import math, cmath, os
import numpy as np
import mpmath as mp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

os.makedirs("figures", exist_ok=True)
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
BLUE, GREEN, ORANGE, RED, GREY, PURPLE = "#2b6cb0", "#2f855a", "#dd6b20", "#c53030", "#718096", "#6b46c1"

def delta_M(g, kw=1.0):
    b = g/(2*kw); return np.array([[1 - 1j*b, -1j*b], [1j*b, 1 + 1j*b]])
def prop(ph): return np.array([[np.exp(1j*ph), 0], [0, np.exp(-1j*ph)]])
def odds(M): return abs(M[0, 1])**2
def length(M): return math.acosh(abs(M[0, 0])**2 + abs(M[0, 1])**2)
def mixer(s, pa=0.0, pb=0.0):
    a = math.sqrt(1 - s)*cmath.exp(1j*pa); b = math.sqrt(s)*cmath.exp(1j*pb)
    return np.array([[a, b], [-b.conjugate(), a.conjugate()]])
def circ(k, a, b): return a + b + k*a*b

# ---------------------------------------------------------------- Fig 1: one law, three curvatures
fig, ax = plt.subplots(1, 3, figsize=(12, 4.2))
# (a) sphere: orthographic view of a spherical triangle
a0 = ax[0]; a0.set_aspect("equal"); a0.axis("off")
a0.add_patch(Circle((0, 0), 1, fill=False, color=GREY))
def sph(lat, lon):
    x = math.cos(lat)*math.sin(lon); y = math.sin(lat); z = math.cos(lat)*math.cos(lon)
    return np.array([x, y, z])
def view(p):  # rotate a little for a 3-D look
    tilt = 0.35; c, s = math.cos(tilt), math.sin(tilt)
    return np.array([p[0], c*p[1] - s*p[2]]), c*p[2] + s*p[1]
def arc(p, q, n=60):
    out = []
    for i in range(n + 1):
        w = i/n; v = (1 - w)*p + w*q; out.append(view(v/np.linalg.norm(v))[0])
    return np.array(out)
P0, P1, P2 = sph(0.9, -0.2), sph(0.05, -0.6), sph(0.15, 0.55)
for p, q, col, lab in [(P0, P1, BLUE, "θ₁"), (P0, P2, BLUE, "θ₂"), (P1, P2, "k", "θ")]:
    A = arc(p, q); a0.plot(A[:, 0], A[:, 1], color=col, lw=2.2 if col != "k" else 1.6)
    m = A[len(A)//2]; a0.text(m[0]*1.08, m[1]*1.08 - 0.03, lab, color=col, ha="center", fontsize=11)
a0.text(view(P0)[0][0], view(P0)[0][1] + 0.07, "φ", ha="center", color=RED, fontsize=12)
a0.set_title("(a) sphere, K = +1\n mean law ∘₋₂ on s = sin²(θ/2)", fontsize=10)
a0.text(0, -1.25, "two-mode mixers with a random phase,\nbit flips (piling-up), freely rotating polymers",
        ha="center", fontsize=8.5, color=GREY)
a0.set_xlim(-1.3, 1.3); a0.set_ylim(-1.5, 1.25)
# (b) flat
a1 = ax[1]; a1.set_aspect("equal"); a1.axis("off")
O = np.array([0, 0]); A_ = np.array([1.1, 0.25]); B_ = np.array([0.35, 0.95])
a1.plot(*zip(O, A_), color=BLUE, lw=2.2); a1.plot(*zip(O, B_), color=BLUE, lw=2.2); a1.plot(*zip(A_, B_), color="k", lw=1.6)
a1.text(0.6, 0.02, "D₁", color=BLUE); a1.text(0.05, 0.55, "D₂", color=BLUE); a1.text(0.8, 0.67, "D", color="k")
a1.text(0.12, 0.12, "φ", color=RED, fontsize=12)
a1.set_title("(b) flat, K = 0\n mean law ∘₀ (plain +) on t = D²/4", fontsize=10)
a1.text(0.55, -0.45, "random walks: variances add", ha="center", fontsize=8.5, color=GREY)
a1.set_xlim(-0.3, 1.4); a1.set_ylim(-0.6, 1.2)
# (c) hyperbolic: Poincare disk, geodesic triangle 0, u1, u1 (+) u2
a2 = ax[2]; a2.set_aspect("equal"); a2.axis("off")
a2.add_patch(Circle((0, 0), 1, fill=False, color=GREY))
def mobius_add(u, v): return (u + v)/(1 + np.conj(u)*v)
def geod(p, q, n=80):  # geodesic in the disk: map p to 0, straight line, map back
    out = []
    w = (q - p)/(1 - np.conj(p)*q)
    for i in range(n + 1):
        z = w*i/n; out.append((z + p)/(1 + np.conj(p)*z))
    return np.array(out)
u1 = 0.7 + 0.0j; u2 = 0.75*cmath.exp(1j*1.75)
V = mobius_add(u1, u2)
for p, q, col, lab in [(0, u1, BLUE, "D₁"), (u1, V, BLUE, "D₂"), (0, V, "k", "D")]:
    G = geod(p, q); a2.plot(G.real, G.imag, color=col, lw=2.2 if col != "k" else 1.6)
    m = G[len(G)//2]; a2.text(m.real + 0.04, m.imag - 0.09 if lab == "D₁" else m.imag + 0.03, lab, color=col)
a2.text(u1.real + 0.03, u1.imag - 0.12, "φ", color=RED, fontsize=12)
a2.plot([0, u1.real, V.real], [0, u1.imag, V.imag], "o", color="k", ms=3)
a2.set_title("(c) hyperbolic, K = −1\n mean law ∘₂ on t = sinh²(D/2) = R/T", fontsize=10)
a2.text(0, -1.25, "barriers: the disk of reflection amplitudes", ha="center", fontsize=8.5, color=GREY)
a2.set_xlim(-1.3, 1.3); a2.set_ylim(-1.5, 1.25)
fig.suptitle("Fig 1. Average a triangle over its angle φ: the squared half-chord composes by ∘ₖ with k = −2K", fontsize=11)
fig.tight_layout(); fig.savefig("figures/fig1-one-law-three-curvatures.png", dpi=150); plt.close(fig)

# ---------------------------------------------------------------- Fig 2: three averages = three arithmetics
fig, ax = plt.subplots(1, 3, figsize=(13, 4))
# (a) Jensen: two paths
r = np.linspace(0.01, 3, 400); phi = np.linspace(0, 2*np.pi, 4000, endpoint=False)
mean_l = [np.log(np.abs(1 + x*np.exp(1j*phi))**2).mean() for x in r]
mean_p = [np.mean(np.abs(1 + x*np.exp(1j*phi))**2) for x in r]
ax[0].plot(r, np.log(mean_p), color=GREEN, lw=2, label="ln ⟨|A + B e^(iφ)|²⟩ = ln(|A|² + |B|²)\n(probabilities add: classical)")
ax[0].plot(r, mean_l, color=BLUE, lw=2.5, label="⟨ln |A + B e^(iφ)|²⟩")
ax[0].plot(r, np.log(np.maximum(1, r)**2), "--", color="k", lw=1, label="ln max(|A|², |B|²)  (tropical)")
ax[0].set_xlabel("|B| / |A|   (A = 1)"); ax[0].set_ylabel("log weight")
ax[0].set_title("(a) two paths, random phase: the log-average\n is the max exactly (Jensen)", fontsize=10)
ax[0].legend(fontsize=7.5, frameon=False, loc="upper left")
# (b) sphere: which path wins, with the wall
S1, S2 = np.meshgrid(np.linspace(0, 1, 300), np.linspace(0, 1, 300))
win = np.where((1 - S1)*(1 - S2) >= S1*S2, 0, 1)
ax[1].contourf(S1, S2, win, levels=[-0.5, 0.5, 1.5], colors=["#bee3f8", "#feebc8"])
ax[1].plot([0, 1], [1, 0], color=RED, lw=2.5)
ax[1].text(0.22, 0.25, "stay–stay wins\n⟨ln(1−s)⟩ adds:\n∘₋₁ on s", ha="center", fontsize=9)
ax[1].text(0.77, 0.73, "switch–switch wins\n⟨ln(1−s)⟩ = ln s₁s₂", ha="center", fontsize=9)
ax[1].text(0.62, 0.44, "wall s₁ + s₂ = 1", color=RED, rotation=-45, fontsize=9)
ax[1].set_xlabel("s₁ (switch probability of mixer 1)"); ax[1].set_ylabel("s₂")
ax[1].set_title("(b) two mixers: the typical law ∘₋₁\n holds up to a wall", fontsize=10)
ax[1].set_aspect("equal")
# (c) barriers: one pair over the phase with the three laws
M1, M2 = delta_M(1.3), delta_M(0.8); t1, t2 = odds(M1), odds(M2)
ph = np.linspace(0, np.pi, 600)
rho = np.array([odds(M2 @ prop(p) @ M1) for p in ph])
ax[2].plot(ph, rho, color="k", lw=1.6, label="the pair's odds t at fixed phase")
ax[2].axhline(circ(2, t1, t2), color=GREEN, lw=2, label=f"mean ⟨t⟩ = t₁ ∘₂ t₂")
ax[2].axhline(circ(1, t1, t2), color=BLUE, lw=2, label="typical: t₁ ∘₁ t₂ (log-average)")
ax[2].axhline(t1 + t2, color=ORANGE, lw=2, label="1/⟨T⟩ − 1 = t₁ + t₂ (∘₀)")
ax[2].set_xlabel("loop phase φ (half a round trip)"); ax[2].set_ylabel("odds of reflection t = R/T")
ax[2].set_title("(c) two barriers: no wall, the direct\n path always wins (cosh·cosh > sinh·sinh)", fontsize=10)
ax[2].legend(fontsize=7.5, frameon=False, loc="upper right")
fig.suptitle("Fig 2. One interference, three arithmetics: amplitudes at fixed phase, + on average, max on log-average", fontsize=11)
fig.tight_layout(); fig.savefig("figures/fig2-three-arithmetics.png", dpi=150); plt.close(fig)

# ---------------------------------------------------------------- Fig 3: edges and the twist
fig, ax = plt.subplots(1, 2, figsize=(11, 5))
D1, D2 = length(M1), length(M2)
ph = np.linspace(0.001, np.pi - 0.001, 600)
amp = np.array([math.tanh(length(M2 @ prop(p) @ M1)/2) for p in ph])
r1, r2 = math.tanh(D1/2), math.tanh(D2/2)
ax[0].plot(ph, amp, color="k", lw=1.8, label="|r| of the pair at fixed phase")
ax[0].axhline((r1 + r2)/(1 + r1*r2), color=PURPLE, lw=2, ls="--", label="(r₁ + r₂)/(1 + r₁r₂): law (k=0, d=−1)")
ax[0].axhline(abs(r1 - r2)/(1 - r1*r2), color=PURPLE, lw=1.2, ls=":", label="(r₁ − r₂)/(1 − r₁r₂): the same law in reverse")
ax[0].set_xlabel("loop phase φ"); ax[0].set_ylabel("reflection amplitude |r| = tanh(D/2)")
ax[0].set_title("(a) the band's edges are relativity's velocity addition\n on reflection amplitudes", fontsize=10)
ax[0].legend(fontsize=8, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=1)
def hyp_area(a, b, c):
    def ang(x, y, z): return math.acos(max(-1, min(1, (math.cosh(x)*math.cosh(y) - math.cosh(z))/(math.sinh(x)*math.sinh(y)))))
    return math.pi - (ang(b, c, a) + ang(a, c, b) + ang(a, b, c))
def wrap(x): return (x + math.pi) % (2*math.pi) - math.pi
exc, half = [], []
for p in ph:
    M = M2 @ prop(p) @ M1
    exc.append(wrap(np.angle(M[0, 0]) - np.angle(M1[0, 0]) - np.angle(M2[0, 0]) - p))
    half.append(hyp_area(D1, D2, length(M))/2)
exc, half = np.array(exc), np.array(half)
ax[1].plot(ph, exc, color="k", lw=1.8, label="extra transmission phase")
ax[1].plot(ph, half, color=RED, lw=1, ls="--", label="+ ½ area of the barrier triangle")
ax[1].plot(ph, -half, color=RED, lw=1, ls=":", label="− ½ area")
ax[1].axhline(0, color=GREY, lw=0.6)
ax[1].set_xlabel("loop phase φ"); ax[1].set_ylabel("radians")
ax[1].set_title("(b) the twist: the phase beyond the parts is half the\n triangle's area (holonomy ω, Wigner's angle)", fontsize=10)
ax[1].legend(fontsize=8, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=3)
fig.suptitle("Fig 3. Fixed phase: the edges read the curvature dial d = K; the twist reads the holonomy", fontsize=11)
fig.tight_layout(); fig.savefig("figures/fig3-edges-and-twist.png", dpi=150); plt.close(fig)

# ---------------------------------------------------------------- Fig 4: Legendre ladder and Legendre's K
mp.mp.dps = 15
fig, ax = plt.subplots(1, 3, figsize=(13, 4))
nus = np.linspace(-1.5, 2, 300)
ax[0].plot(nus, nus*(nus + 1), color="k", lw=1.5)
for nu, lab, col in [(1, "ν = 1: mean odds, ∘₂\nrate 2", GREEN), (-0.5, "ν = −½: mean T\nrate ¼ (bottom of\nthe spectrum)", ORANGE)]:
    ax[0].plot(nu, nu*(nu + 1), "o", color=col, ms=8); ax[0].annotate(lab, (nu, nu*(nu + 1)), xytext=((nu + 0.15) if nu == 1 else -1.45, (nu*(nu + 1) + 1.2) if nu == 1 else 3.6), fontsize=8.5, color=col)
ax[0].plot(0, 0, "o", color=BLUE, ms=8); ax[0].annotate("slope at ν = 0: typical, ∘₁\nrate 1", (0, 0), xytext=(-1.45, 1.9), fontsize=8.5, color=BLUE)
ax[0].plot([-0.4, 0.4], [-0.4, 0.4], color=BLUE, lw=1.5)
ax[0].set_xlabel("Legendre index ν"); ax[0].set_ylabel("ν(ν + 1): growth per unit of typical length")
ax[0].set_title("(a) the three standard averages are three\n points of one Legendre family", fontsize=10)
tt = np.linspace(0.0, 6, 200)
PK = [float(mp.legenp(-0.5, 0, 1 + 2*x, type=3)) for x in tt]
Kf = [2/math.pi*math.sqrt(1/(1 + x))*float(mp.ellipk(x/(1 + x))) for x in tt[::10]]
ax[1].plot(tt, PK, color=ORANGE, lw=2, label="P₋½(1 + 2t)")
ax[1].plot(tt[::10], Kf, "o", color="k", ms=4, label="(2/π)·√T·K(|r|)")
ax[1].plot(tt, 1 - tt/4, "--", color=GREY, lw=1, label="1 − t/4 (weak barrier)")
ax[1].set_ylim(0, 1.05)
ax[1].set_xlabel("one barrier's odds t"); ax[1].set_ylabel("per-barrier factor of ⟨T⟩")
ax[1].set_title("(b) mean transmission falls per barrier by\n Legendre's elliptic integral K of |r|", fontsize=10)
ax[1].legend(fontsize=8, frameon=False)
MB = delta_M(1.0); tB = odds(MB); zB = 1 + 2*tB
lams = np.linspace(0, 8, 801)
PcB = np.array([float(mp.re(mp.legenp(-0.5 + 1j*l, 0, zB, type=3))) for l in lams])
w = lams*np.tanh(np.pi*lams)*2*np.pi/np.cosh(np.pi*lams)
Ns = np.arange(1, 61)
mT = [np.trapezoid(w*PcB**N, lams) for N in Ns]
PB = float(mp.legenp(-0.5, 0, zB, type=3))
ax[2].semilogy(Ns, mT, color=ORANGE, lw=2, label="⟨T_N⟩ from the Legendre spectrum")
ax[2].semilogy(Ns, [(1/(1 + tB))**N for N in Ns], color=BLUE, lw=1.5, label="typical T: (1 + t₁)^−N")
ax[2].semilogy(Ns, [1/(1 + ((1 + 2*tB)**N - 1)/2) for N in Ns], color=GREEN, lw=1.5, label="1/(1 + mean odds)")
rng = np.random.default_rng(3)
for N in (2, 4, 8, 16):
    trials = 40000; M = np.tile(MB, (trials, 1, 1))
    for j in range(1, N):
        p = rng.uniform(0, np.pi, trials); Pm = np.zeros((trials, 2, 2), complex); Pm[:, 0, 0] = np.exp(1j*p); Pm[:, 1, 1] = np.exp(-1j*p)
        M = MB @ Pm @ M
    ax[2].semilogy(N, (1/np.abs(M[:, 1, 1])**2).mean(), "o", color="k", ms=4)
ax[2].set_xlabel("number of barriers N (t₁ = ¼, random spacings)"); ax[2].set_ylabel("transmission")
ax[2].set_title("(c) dots: simulated ⟨T⟩; three averages,\n three different decay rates", fontsize=10)
ax[2].legend(fontsize=8, frameon=False)
fig.suptitle("Fig 4. The Legendre ladder: mean odds, typical value and mean transmission", fontsize=11)
fig.tight_layout(); fig.savefig("figures/fig4-legendre-ladder.png", dpi=150); plt.close(fig)

# ---------------------------------------------------------------- Fig 5: rewinds
fig, ax = plt.subplots(1, 3, figsize=(13, 4))
zc = 3.0
lam = np.linspace(0, 12, 600)
ch = np.array([float(mp.re(mp.legenp(-0.5 + 1j*l, 0, zc, type=3))) for l in lam])
ax[0].plot(lam, ch, color="k", lw=1.6); ax[0].axhline(0, color=GREY, lw=0.6)
zs = [lam[i] for i in range(len(lam) - 1) if ch[i]*ch[i + 1] < 0]
ax[0].plot(zs, [0]*len(zs), "o", color=RED, ms=6, label="erased components")
ax[0].set_xlabel("spectral parameter λ"); ax[0].set_ylabel("barrier character, t₁ = 1")
ax[0].set_title("(a) one barrier multiplies each component\n of the spread by this; zeros are erased", fontsize=10)
ax[0].legend(fontsize=8, frameon=False)
X, Y = np.meshgrid(np.linspace(0.05, 8, 160), np.linspace(-1.5, 1.5, 80))
Z = np.zeros_like(X)
for i in range(X.shape[0]):
    for j in range(X.shape[1]):
        Z[i, j] = abs(complex(mp.legenp(-0.5 + 1j*complex(X[i, j], Y[i, j]), 0, zc, type=3)))
ax[1].contourf(X, Y, np.log10(Z + 1e-6), levels=30, cmap="Blues_r")
ax[1].plot([z for z in zs if z < 8], [0]*len([z for z in zs if z < 8]), "o", color=RED, ms=5)
ax[1].set_xlabel("Re λ"); ax[1].set_ylabel("Im λ")
ax[1].set_title("(b) log |character| over complex λ:\n every zero sits on the real line (Pólya)", fontsize=10)
s1 = np.linspace(0, 1, 400)
ax[2].plot(s1, 1 - 2*s1, color=GREEN, lw=2, label="mixers: bias factor 1 − 2s₁ (rewind divides by it)")
ax[2].plot(s1, 1 + 2*s1, color=PURPLE, lw=2, label="barriers: factor 1 + 2t₁, never zero")
ax[2].axhline(0, color=GREY, lw=0.6); ax[2].axvline(0.5, color=RED, lw=1, ls="--")
ax[2].text(0.52, 1.8, "wall: a 50/50 mixer\ncannot be rewound", color=RED, fontsize=8.5)
ax[2].set_xlabel("s₁ or t₁"); ax[2].set_ylabel("factor on the mean")
ax[2].set_title("(c) single numbers rewind exactly, except\n through a wall", fontsize=10)
ax[2].legend(fontsize=7.5, frameon=False, loc="lower left")
fig.suptitle("Fig 5. Rewinds: the mean and the typical value rewind; the full spread loses a comb of components", fontsize=11)
fig.tight_layout(); fig.savefig("figures/fig5-rewinds.png", dpi=150); plt.close(fig)

# ---------------------------------------------------------------- Fig 6: the dial line
fig, ax = plt.subplots(figsize=(11, 3.4))
ax.axhline(0, color="k", lw=1.2)
marks = [(-2, "−2", "sphere mean\n(bit flips, piling-up,\nmixers, polymers)", GREEN),
         (-1, "−1", "sphere typical (up to the wall);\nbarriers' blocking chance R,\ntypical (k = 1 read on R)", BLUE),
         (0, "0", "flat (variances add);\npair mean T; the band's\nedges (with d = K)", ORANGE),
         (1, "1", "barriers typical\n(always exact)", BLUE),
         (2, "2", "barriers mean\n(Anderson–Thouless–\nAbrahams–Fisher)", GREEN)]
for x, lab, txt, col in marks:
    ax.plot(x, 0, "o", color=col, ms=11); ax.text(x, 0.12, "k = " + lab, ha="center", fontsize=10, weight="bold")
    ax.text(x, -0.18, txt, ha="center", va="top", fontsize=8.5, color=col)
ax.annotate("", xy=(2.6, 0.45), xytext=(-2.6, 0.45), arrowprops=dict(arrowstyle="<->", color=GREY))
ax.text(-2.6, 0.52, "positive curvature (sphere)", fontsize=8.5, color=GREY)
ax.text(2.6, 0.52, "negative curvature (hyperbolic)", fontsize=8.5, color=GREY, ha="right")
ax.text(0, 0.52, "mean: k = −2K   ·   typical: k = −K", fontsize=9, ha="center")
ax.set_xlim(-2.9, 2.9); ax.set_ylim(-0.95, 0.7); ax.axis("off")
ax.set_title("Fig 6. Where tunnelling and its twin put the homothety k", fontsize=11)
fig.tight_layout(); fig.savefig("figures/fig6-dial-line.png", dpi=150); plt.close(fig)
print("figures written")
