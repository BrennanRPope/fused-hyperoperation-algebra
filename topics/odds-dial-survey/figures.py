"""Figures for the odds-dial survey. Run: python3 figures.py  -> figures/*.png  (about 15 s)
Colour key: blue = our reading, orange = already known, grey = endpoints/background."""
import json, math, itertools
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 14, "axes.titlesize": 15, "axes.labelsize": 14, "legend.fontsize": 12,
                     "axes.spines.top": False, "axes.spines.right": False})
B, O, G = "#1f5fbf", "#e07b00", "#8a8a8a"
DPI = 130
def save(fig, name):
    fig.tight_layout(); fig.savefig("figures/" + name, dpi=DPI); plt.close(fig); print(name)
circ = lambda a, b: a + b + a*b
r2, phi = math.sqrt(2) - 1, (1 + math.sqrt(5))/2
def subsets(n): return itertools.product([0, 1], repeat=n)

# ---- fig0 hero: the dial and its proved stops ----
fig, ax = plt.subplots(figsize=(12, 5.6))
t = np.logspace(-3, 3, 600)
ax.plot(t, t/(1+t), color=G, lw=2.5, label="chance t/(1 + t)  (a reading)")
ax.plot(t, -(t/(1+t)*np.log2(t/(1+t)) + 1/(1+t)*np.log2(1/(1+t))), color=B, lw=2.5, label="bits per element")
ax.set_xscale("log"); ax.set_xlim(1e-3, 1e3); ax.set_ylim(-0.02, 1.32)
for x, c in [(1e-3, G), (1/phi, G), (r2, B), (1.0, O), (1e3, G)]:
    ax.axvline(x, color=c, ls="--", lw=1.5)
ax.text(1.2e-3, 1.08, "t → 0\nelement dropped\n(empty set)", ha="left", fontsize=11.5, color=G)
ax.text(9e2, 1.08, "t → ∞\nelement forced\n(full set)", ha="right", fontsize=11.5, color=G)
ax.text(0.39, 1.08, "√2 − 1\nunion half", ha="right", fontsize=12, color=B)
ax.text(1.08, 1.08, "t = 1\nplain power set,\nself-dual", ha="left", fontsize=12, color=O)
ax.text(0.64, 0.66, "1/φ", ha="left", fontsize=11.5, color=G)
ax.set_xlabel("odds t of one element  (log scale; the mirror t ↔ 1/t is the complement)")
ax.set_ylabel("chance, bits")
ax.annotate("small odds: Fermi = Boltzmann = Bose,\nbinomial → Poisson, (1+1/n)ⁿ → e", (3e-3, 0.03), (1.5e-3, 0.45), fontsize=11.5, color=O,
            arrowprops=dict(arrowstyle="->", color=O))
ax.text(30, 0.33, "t = 1: chance ½, Bayes neutral,\nfermion level half full, boson wall,\nRT mutual information switches on,\nLee–Yang zeros pinch here", fontsize=11.5, color=O)
ax.legend(loc="center left", bbox_to_anchor=(0.0, 0.62), frameon=False)
ax.set_title("The odds dial t: one number per element, with its proved stops")
save(fig, "fig0-odds-dial.png")

# ---- sorting board ----
SHORT = {3: "powers (1+t)^n", 9: "velocity addition", 10: "Bayes odds", 17: "q-counting", 18: "partitions",
         20: "graphs: reliability", 22: "logistic, softmax", 23: "belief propagation", 24: "Shannon entropy",
         25: "repetition dial", 26: "thermodynamics", 27: "Lee-Yang circle", 31: "selection, Elo", 32: "point processes",
         40: "e = lim (1+1/n)^n", 44: "binding, pH", 58: "holographic gluing", 61: "primes, squarefree",
         2: "corner expansion", 15: "power-set dials", 21: "union rule (odds_or)"}
L = json.load(open("laws.json"))["rows"]
order = ["slides", "self", "no"]
colors = {"slides": B, "self": O, "no": "#d0d0d0"}
fig, ax = plt.subplots(figsize=(12, 7.5)); ax.axis("off")
y0 = 0.97
for v in order:
    rows = [r for r in L if r["verdict"] == v]
    ax.text(0.0, y0, f"{v}: {len(rows)}", fontsize=16, weight="bold", color=colors[v] if v != "no" else G, transform=ax.transAxes)
    y0 -= 0.05
    if v == "no":
        per = 15
        for i, r in enumerate(rows):
            x = (i % per)/per; y = y0 - (i//per)*0.06
            ax.add_patch(plt.Rectangle((x, y-0.045), 1/per-0.006, 0.05, color=colors[v], transform=ax.transAxes))
            ax.text(x+0.5/per, y-0.02, str(r["id"]), ha="center", va="center", fontsize=11, transform=ax.transAxes)
        y0 -= 0.06*((len(rows)-1)//per + 1)
    else:
        per = 3
        for i, r in enumerate(rows):
            x = (i % per)/per; y = y0 - (i//per)*0.055
            ax.add_patch(plt.Rectangle((x, y-0.045), 1/per-0.006, 0.05, color=colors[v], alpha=0.85, transform=ax.transAxes))
            lab = f"{r['id']} {r['candidate']}: " + SHORT.get(r["id"], r["law"][:24])
            ax.text(x+0.005, y-0.02, lab, va="center", fontsize=10.5, color="white", transform=ax.transAxes)
        y0 -= 0.055*((len(rows)-1)//per + 1) + 0.02
ax.text(0, 0.0, "No group is 'plausible', 'partly' or 'disguise'. Numbers in grey are catalog ids (dial-catalog/catalog.py).", fontsize=12, color=G, transform=ax.transAxes)
ax.set_title("All 65 law groups, sorted by what turning the odds t does", fontsize=16)
save(fig, "fig-sorting-board.png")

# ---- C1 ----
fig, axs = plt.subplots(1, 2, figsize=(12, 4.8))
t = np.logspace(-2, 2, 400)
axs[0].plot(t, np.log2(1+t), color=O, lw=2.5)
for x, a_, lab in [(r2, 0.5, "√2 − 1 → half (a = ½)"), (1, 1, "1 → whole power set (a = 1)"), (3, 2, "3 → P twice (a = 2)")]:
    axs[0].plot(x, a_, "o", color=B, ms=9); axs[0].annotate(lab, (x, a_), (x*0.07, a_+0.6), fontsize=12, color=B, arrowprops=dict(arrowstyle="->", color=B))
axs[0].set_xscale("log"); axs[0].set_xlabel("uniform odds t"); axs[0].set_ylabel("exponent a = log₂(1 + t)")
axs[0].set_title("Uniform t is the fractional exponent")
n = 6; T = np.logspace(-2, 3, 300)
axs[1].plot(T, (1/(1+T))**n, color=G, lw=2.5, label="P(empty set)")
axs[1].plot(T, (T/(1+T))**n, color=G, lw=2.5, ls="--", label="P(everything)")
axs[1].plot(T, [math.comb(n, 3)*(x/(1+x))**3/(1+x)**3 for x in T], color=B, lw=2.5, label="P(exactly half), peak at t = 1")
axs[1].set_xscale("log"); axs[1].set_xlabel("uniform odds t (n = 6)"); axs[1].set_ylabel("chance"); axs[1].legend(frameon=False, loc="center left")
axs[1].set_title("Empty ↔ plain ↔ full")
save(fig, "fig1-power-set-stops.png")

# ---- C2 ----
rng = np.random.default_rng(3)
fig, axs = plt.subplots(1, 2, figsize=(12, 4.8))
for rep in range(5):
    Ls = np.exp(rng.normal(0.15, 0.5, 40)); lo = np.log(0.2/0.8) + np.concatenate([[0], np.cumsum(np.log(Ls))])
    axs[0].plot(range(41), 1/(1+np.exp(-lo)), color=B, alpha=0.8)
axs[0].set_xlabel("looks at the evidence"); axs[0].set_ylabel("posterior chance"); axs[0].set_title("Bayes: each look multiplies the odds")
O_ = np.logspace(-2, 2, 300)
axs[1].plot(np.log(O_), (O_-1)/(O_+1), color=O, lw=2.5)
axs[1].set_xlabel("log odds  =  rapidity"); axs[1].set_ylabel("velocity (O − 1)/(O + 1)  (units of c)")
axs[1].set_title("The same odds as a boost")
axs[1].text(-4.4, 0.55, "multiply odds\n= add rapidities\n= Einstein addition", color=O, fontsize=12)
save(fig, "fig2-bayes-boost.png")

# ---- C3 ----
fig, axs = plt.subplots(1, 2, figsize=(12, 4.8))
t = np.logspace(-3, 3, 400)
axs[0].plot(t, t/(1+t), color=O, lw=2.5, label="chance = t/(1+t): logistic in log t")
axs[0].plot(t, t**2/(1+t**2), color=B, lw=2, ls="--", label="two cooperative sites (ω large): steeper")
axs[0].axvline(1, color=G, ls=":"); axs[0].set_xscale("log")
axs[0].set_xlabel("odds t  (e.g. [ligand]/K, or 10^(pH − pKa))"); axs[0].set_ylabel("fraction filled"); axs[0].legend(frameon=False, loc="upper left", fontsize=11)
axs[0].set_title("One reading, many fields")
t1 = np.logspace(-2, 2, 300); others = np.array([1.0, 0.5])
axs[1].plot(t1, np.full_like(t1, others[0]/(1+others[0])), color=O, lw=2.5, label="cube (independent yes/no)")
axs[1].plot(t1, others[0]/(1+t1+others.sum()), color=B, lw=2.5, label="simplex (one of many, with 'none')")
axs[1].set_xscale("log"); axs[1].set_xlabel("odds t₁ of element 1, turned alone"); axs[1].set_ylabel("chance of element 2")
axs[1].legend(frameon=False); axs[1].set_title("Turning one element alone")
save(fig, "fig3-logistic-cube-simplex.png")

# ---- C4 ----
fig, axs = plt.subplots(1, 2, figsize=(12, 4.8))
t = np.logspace(-3, 0.5, 400)
axs[0].plot(t, t/(1+t), color=O, lw=2.5, label="Fermi (k = 1)")
axs[0].plot(t, t, color=G, lw=2.5, label="Boltzmann (k = 0)")
tb = t[t < 0.98]; axs[0].plot(tb, tb/(1-tb), color=O, lw=2.5, ls="--", label="Bose (k = −1): wall at t = 1")
axs[0].set_xscale("log"); axs[0].set_yscale("log"); axs[0].set_ylim(1e-3, 50)
axs[0].set_xlabel("odds t = fugacity × Boltzmann factor"); axs[0].set_ylabel("mean copies n = t/(1 + k t)")
axs[0].legend(frameon=False); axs[0].set_title("Small odds: the three agree")
from scipy.stats import binom, poisson
Ns = np.array([5, 10, 20, 50, 100, 200, 500, 1000, 3000]); kk = np.arange(60)
tv = [0.5*np.abs(binom.pmf(kk, N, (2/N)/(1+2/N)) - poisson.pmf(kk, 2)).sum() for N in Ns]
axs[1].loglog(Ns, tv, "o-", color=B, lw=2); axs[1].loglog(Ns, 1/Ns, color=G, ls=":", label="1/N")
axs[1].set_xlabel("N elements, each at odds 2/N"); axs[1].set_ylabel("distance to Poisson(2)"); axs[1].legend(frameon=False)
axs[1].set_title("Binomial → Poisson as t → 0")
save(fig, "fig4-classical-limit.png")

# ---- C5 ----
res = json.load(open("results.json"))["C5"]
def gas_poly(n_sites, edges, k):
    c = np.zeros(n_sites+1)
    for S in subsets(n_sites):
        c[sum(S)] += k**sum(S[i] and S[j] for i, j in edges)
    return c
L4 = 4; tor = []
for i in range(L4):
    for j in range(L4):
        s0 = i*L4+j; tor += [(s0, i*L4+(j+1) % L4), (s0, ((i+1) % L4)*L4+j)]
fig, axs = plt.subplots(1, 3, figsize=(14, 4.8))
th = np.linspace(0, 2*np.pi, 300)
axs[0].plot(np.cos(th), np.sin(th), color=O, lw=1.5)
z = np.array(res["torus"]); axs[0].plot(z[:, 0], z[:, 1], "o", color=B, ms=6, label="zeros")
axs[0].plot([1], [0], "*", color=O, ms=16)
axs[0].set_aspect("equal"); axs[0].set_xlabel("Re normalised odds"); axs[0].set_ylabel("Im"); axs[0].legend(frameon=False, loc="center", fontsize=11)
axs[0].set_title("4×4 torus, k = 2: zeros on |t·k²| = 1", fontsize=13)
zz = np.logspace(-1.2, 1.2, 300)
for kk_, c in [(1.5, G), (3.0, B), (6.0, O)]:
    cz = gas_poly(16, tor, kk_) * (kk_**2)**(-np.arange(17))
    rho = [sum(i*ci*x**i for i, ci in enumerate(cz))/sum(ci*x**i for i, ci in enumerate(cz))/16 for x in zz]
    axs[1].plot(zz, rho, color=c, lw=2.3, label=f"k = {kk_}")
axs[1].set_xscale("log"); axs[1].axvline(1, color=G, ls=":")
axs[1].set_xlabel("normalised real odds"); axs[1].set_ylabel("filled fraction"); axs[1].legend(frameon=False)
axs[1].set_title("Steepest at odds 1")
hc = np.roots(gas_poly(12, [(i, (i+1) % 12) for i in range(12)], 0.0)[::-1])
axs[2].plot(hc.real, hc.imag, "o", color=O, ms=8); axs[2].axvline(-0.25, color=B, ls="--")
axs[2].text(-0.24, 0.2, "−1/4", color=B); axs[2].set_xscale("symlog", linthresh=0.1)
axs[2].set_xlabel("odds t (hard core, k = 0, ring of 12)"); axs[2].set_ylabel("Im"); axs[2].set_ylim(-0.5, 0.5)
axs[2].set_title("Repulsion: zeros on t < −1/4")
save(fig, "fig5-lee-yang.png")

# ---- C6 ----
from sympy import primerange
import mpmath as mp
P = np.array(list(primerange(2, 20000)), dtype=float)
fig, axs = plt.subplots(1, 2, figsize=(12, 4.8))
for s, c in [(1.5, G), (2.0, B), (3.0, O)]:
    part = np.cumprod(1 + P**-s); axs[0].semilogx(P, part, color=c, lw=2.3, label=f"s = {s}")
    axs[0].axhline(float(mp.zeta(s)/mp.zeta(2*s)), color=c, ls=":")
axs[0].set_xlabel("primes included up to"); axs[0].set_ylabel("∏ (1 + p^−s)"); axs[0].legend(frameon=False)
axs[0].set_title("Odds t_p = p^−s per prime")
ss = np.linspace(1.05, 4, 200)
axs[1].plot(ss, [float(mp.zeta(s)/mp.zeta(2*s)) for s in ss], color=O, lw=2.5)
axs[1].axvline(1, color=G, ls="--"); axs[1].text(1.1, 9, "s = 1: Σ t_p diverges\n(the pole)", color=G, fontsize=12)
axs[1].plot([2], [15/math.pi**2], "o", color=B, ms=9); axs[1].text(2.1, 2.2, "15/π² at s = 2", color=B)
axs[1].set_ylim(0, 12); axs[1].set_xlabel("s  (all odds shrink as s grows)"); axs[1].set_ylabel("ζ(s)/ζ(2s)")
axs[1].set_title("Squarefree numbers along the dial")
save(fig, "fig6-primes.png")

# ---- C7 and N2 ----
lp = lambda X: np.where(X > 1, np.log(np.maximum(X, 1e-300)), 0.0)
fig, axs = plt.subplots(1, 2, figsize=(12, 5))
X = np.logspace(-1.5, 1.5, 400)
axs[0].plot(X, lp(X), color=O, lw=2.5); axs[0].axvline(1, color=G, ls=":")
axs[0].set_xscale("log"); axs[0].set_xlabel("odds X of the connected RT phase"); axs[0].set_ylabel("I(A:B) / (c/3)")
axs[0].set_title("Mutual information switches on at X = 1")
g = np.linspace(0, 1.4, 400); X1, X2 = np.meshgrid(g, g)
I3 = lp(X1) + lp(X2) - lp(circ(X1, X2))
axs[1].contourf(X1, X2, I3 < -1e-12, levels=[-0.5, 0.5, 1.5], colors=["white", "#cfe0f7"])
xs = np.linspace(0, 1, 200); axs[1].plot(xs, (1-xs)/(1+xs), color=B, lw=2.5, label="partner curve (1 − X₁)/(1 + X₁)")
axs[1].plot([r2], [r2], "o", color=B, ms=10); axs[1].annotate("√2 − 1", (r2, r2), (0.05, 0.12), color=B, fontsize=13, arrowprops=dict(arrowstyle="->", color=B))
axs[1].axvline(1, color=G, ls=":"); axs[1].axhline(1, color=G, ls=":")
axs[1].text(0.45, 0.75, "I₃ < 0", color=B, fontsize=14); axs[1].text(0.05, 0.45, "I₃ = 0", color=G, fontsize=14)
axs[1].set_xlabel("odds X₁ (A with B₁)"); axs[1].set_ylabel("odds X₂ (A with B₂)"); axs[1].legend(frameon=False, loc="upper right", fontsize=11)
axs[1].set_title("Where three-party sharing starts")
save(fig, "fig7-holographic-odds.png")

def X_of(A, Bb):
    (a1, a2), (b1, b2) = A, Bb; return (a2-a1)*(b2-b1)/((b1-a2)*(b2-a1))
gs = np.linspace(0.05, 3, 400); rows = []
for g_ in gs:
    x1, x2 = X_of((-1, 0), (g_, g_+1)), X_of((-1, 0), (g_+1, g_+2)); x12 = circ(x1, x2)
    rows.append((x1, x2, x12, float(lp(np.array(x1)) + lp(np.array(x2)) - lp(np.array(x12)))))
rows = np.array(rows)
fig, axs = plt.subplots(1, 2, figsize=(12, 4.8))
axs[0].semilogy(gs, rows[:, 0], color=G, lw=2, label="X₁ (A, B₁)"); axs[0].semilogy(gs, rows[:, 1], color=G, lw=2, ls="--", label="X₂ (A, B₂)")
axs[0].semilogy(gs, rows[:, 2], color=B, lw=2.5, label="X₁ ∘₁ X₂ (A, B₁B₂)"); axs[0].axhline(1, color=O, ls=":")
axs[0].set_xlabel("gap between A = [−1, 0] and B₁ (B₁, B₂ of length 1)"); axs[0].set_ylabel("odds"); axs[0].legend(frameon=False, fontsize=11)
axs[0].set_title("Real intervals: the glued odds")
axs[1].plot(gs, rows[:, 3], color=B, lw=2.5); axs[1].axhline(0, color=G, lw=1)
zone = (rows[:, 0] < 1) & (rows[:, 2] > 1)
axs[1].axvspan(gs[zone].min(), gs[zone].max(), color="#cfe0f7")
axs[1].text(gs[zone].max()+0.05, -0.25, "shaded: A shares nothing with B₁ or B₂\nalone, yet I₃ < 0 (the glued odds\npassed 1)", color=B, fontsize=11.5)
axs[1].set_xlabel("gap"); axs[1].set_ylabel("I₃ / (c/3)"); axs[1].set_title("I₃ < 0 once the glued odds passes 1")
save(fig, "fig12-new-union-half.png")

# ---- C8 ----
fig, axs = plt.subplots(1, 2, figsize=(12, 4.8))
t = np.logspace(-2, 2, 300)
for kk_, c, ls in [(-1, O, "-"), (0, G, "-"), (0.5, B, "--"), (1, O, "--")]:
    mu = 0.3; mu_t = t*mu/(1 + kk_*mu*(1-t))
    ok = mu_t > 0
    axs[0].loglog(t[ok], mu_t[ok], color=c, ls=ls, lw=2.3, label={-1: "fermion mode (k = −1)", 0: "Poisson (k = 0)", 0.5: "k = 1/2", 1: "boson mode (k = 1)"}[kk_])
axs[0].set_xlabel("tilt t (multiplies the mode's odds)"); axs[0].set_ylabel("mean count μ_t (μ = 0.3)"); axs[0].legend(frameon=False, fontsize=11)
axs[0].set_title("Tilting stays in the family")
mu = np.array([0.1, 0.4, 0.7, 0.9])
for m_ in mu:
    axs[1].semilogx(t, t*m_/(1 - m_ + t*m_), color=B, lw=2)
axs[1].set_xlabel("tilt t"); axs[1].set_ylabel("DPP kernel eigenvalue"); axs[1].set_title("L-ensemble: eigenvalues of L are odds")
axs[1].text(0.012, 0.8, "μ_t/(1 − μ_t) = t · μ/(1 − μ)", color=O, fontsize=13)
save(fig, "fig8-k-copies-tilt.png")

# ---- C9 ----
fig, axs = plt.subplots(1, 2, figsize=(12, 4.8))
nq = 12
for q, c in [(1.0, G), (0.8, B), (0.5, O)]:
    poly = np.array([1.0])
    for i in range(nq): poly = np.convolve(poly, [1, q**i])
    axs[0].semilogy(range(nq+1), poly, "o-", color=c, lw=2, label=f"q = {q}")
axs[0].set_xlabel("subset size j"); axs[0].set_ylabel("weight of size-j subsets"); axs[0].legend(frameon=False)
axs[0].set_title("Odds t_i = qⁱ: binomial → q-binomial")
m = 30; c = [0]*(m+1); c[0] = 1
for part in range(1, m+1):
    for s0 in range(m, part-1, -1): c[s0] += c[s0-part]
axs[1].plot(range(m+1), c, "o-", color=O, lw=2)
axs[1].set_xlabel("m"); axs[1].set_ylabel("partitions of m into distinct parts"); axs[1].set_title("n → ∞ at odds qⁱ: Euler's product")
save(fig, "fig9-geometric-odds.png")

# ---- C10 ----
fig, axs = plt.subplots(1, 2, figsize=(12, 4.8), gridspec_kw={"width_ratios": [1, 1.3]})
pos = {0: (0, 0.5), 1: (0.5, 1), 2: (0.5, 0), 3: (1, 0.5)}
E = [(0, 1), (0, 2), (1, 3), (2, 3), (1, 2)]
for i, (u, v) in enumerate(E):
    axs[0].plot(*zip(pos[u], pos[v]), color=B if i == 4 else G, lw=5 if i == 4 else 3)
for k_, (x, y) in pos.items(): axs[0].plot(x, y, "o", color="k", ms=10)
axs[0].text(0.55, 0.45, "bridge bond:\nturn its odds alone", color=B, fontsize=12); axs[0].text(-0.08, 0.6, "s", fontsize=14); axs[0].text(1.03, 0.6, "t", fontsize=14)
axs[0].axis("off"); axs[0].set_title("Wheatstone bridge")
def rel(pp, pe):
    par_ = lambda x, y: x + y - x*y
    return (1-pe)*par_(pp*pp, pp*pp) + pe*par_(pp, pp)**2
tb = np.logspace(-2, 2, 300); pe = tb/(1+tb)
for pp, c in [(0.3, B), (0.6, O)]:
    axs[1].semilogx(tb, rel(pp, pe), color=c, lw=2.5, label=f"other bonds p = {pp}")
axs[1].set_xlabel("odds of the bridge bond"); axs[1].set_ylabel("chance s connects to t")
axs[1].text(0.012, 0.36, "t → 0: bridge deleted", color=G); axs[1].text(4, 0.36, "t → ∞: contracted", color=G)
axs[1].legend(frameon=False); axs[1].set_title("Deletion ↔ contraction")
save(fig, "fig10-deletion-contraction.png")

# ---- N1: the t = 1 hub ----
fig, ax = plt.subplots(figsize=(12, 5))
t = np.logspace(-2, 2, 400)
ax.plot(np.log(t), t/(1+t), color=G, lw=2.5)
ax.plot(-np.log(t), t/(1+t), color=G, lw=1, ls=":")
ax.axvline(0, color=O, lw=2)
laws = ["complement t ↔ 1/t fixes only t = 1 (Lean: compl_fixed_dial)", "Bayes: likelihood ratio 1 changes nothing",
        "logistic chance ½; entropy at its maximum", "fermion level half full (level and hole swap t, 1/t)",
        "RT mutual information switches on (X = 1)", "Lee–Yang: zeros pinch the real axis only here",
        "(the boson wall is also at t = 1, for another reason)"]
for i, s in enumerate(laws):
    ax.text(0.25, 0.95 - 0.12*i, ("• " if i < 6 else "") + s, fontsize=12.5, color=O if i < 6 else G, va="top")
ax.set_xlabel("log odds  (the mirror log t → −log t is the complement)"); ax.set_ylabel("chance")
ax.set_title("t = 1: one stop, five proved laws, one reason (the mirror t ↔ 1/t)")
save(fig, "fig11-new-self-dual-stop.png")

# ---- N3 golden ----
fig, axs = plt.subplots(1, 2, figsize=(12, 4.6))
p = np.linspace(0, 1, 400)
for nn, c in [(2, B), (3, G), (5, G)]:
    axs[0].plot(p, (1-p)**nn, color=c, lw=2.3, label=f"all {nn} copies miss")
axs[0].plot(p, p, color=O, lw=2.3, label="one copy hits")
axs[0].plot([1/phi**2], [1/phi**2], "o", color=B, ms=10); axs[0].text(0.4, 0.5, "p = 1/φ² (t = 1/φ)", color=B)
axs[0].set_xlabel("chance p of one copy"); axs[0].legend(frameon=False, fontsize=11); axs[0].set_title("The golden point is (1 − p)² = p")
fib = [gas_poly(n_, [(i, i+1) for i in range(n_-1)], 0.0).sum() for n_ in range(2, 15)]
axs[1].semilogy(range(2, 15), fib, "o-", color=O, lw=2); axs[1].semilogy(range(2, 15), [phi**(n_+2)/math.sqrt(5) for n_ in range(2, 15)], color=G, ls=":", label="φⁿ⁺²/√5")
axs[1].set_xlabel("path length n"); axs[1].set_ylabel("hard-core count at t = 1"); axs[1].legend(frameon=False)
axs[1].set_title("φ as a growth rate, at t = 1")
save(fig, "fig13-new-golden.png")
