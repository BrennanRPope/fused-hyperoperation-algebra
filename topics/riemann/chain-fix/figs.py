"""Chain-fix figures. DRAFT."""
import numpy as np, json, sys, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import mpmath as mp
sys.path.insert(0, "../knauf-field"); from knauf import chain
def walsh(f):
    a = f.astype(float).copy(); h = 1
    while h < len(a):
        for i in range(0, len(a), 2*h):
            x = a[i:i+h].copy(); y = a[i+h:i+2*h].copy(); a[i:i+h] = x+y; a[i+h:i+2*h] = x-y
        h *= 2
    return a/len(a)
plt.rcParams.update({"font.size": 10})
def load(*fs):
    d = {}
    for f in fs: d.update({int(k): v for k, v in json.load(open(f)).items()})
    n = np.array(sorted(d)); return n, np.array([d[k] for k in n])
def up_of(n): return np.array([bin(k).count("1") for k in range(2**n)])

# ---------- Fig A: the chain is ferromagnetic but lopsided ----------
n = 8; phi = chain(n).astype(float); H = np.log(phi); J = walsh(H); size = up_of(n)
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
for k in range(1, n+1):
    v = -J[size == k]; v = v[v > 1e-12]
    ax[0].scatter(np.full(v.size, k) + np.random.default_rng(k).uniform(-0.25, 0.25, v.size), v, s=8,
                  color="tab:red" if k % 2 else "tab:blue", alpha=0.6)
ax[0].set_yscale("log"); ax[0].set_xlabel("how many spins a link ties together (order)")
ax[0].set_ylabel("link strength  −J  (log scale)")
ax[0].set_title("A1. Every link pulls spins into line (all −J ≥ 0, n = 8)")
ax[0].text(0.98, 0.95, "blue: even-order links (flip-fair)\nred: odd-order links\n= a field built into the chain",
           transform=ax[0].transAxes, ha="right", va="top", fontsize=9)
for nn, c in ((8, "tab:green"), (12, "tab:purple")):
    p = chain(nn).astype(float); u = up_of(nn)
    ax[1].plot(range(nn+1), [np.log(p[u == k]).mean() for k in range(nn+1)], "o-", color=c, ms=4, label="n = %d" % nn)
ax[1].set_xlabel("number of up spins"); ax[1].set_ylabel("average energy  ln φ")
ax[1].set_title("A2. Flipping every spin changes the energy")
ax[1].annotate("all down: φ = 1", (0, 0), (1.5, 0.4), arrowprops=dict(arrowstyle="->"))
ax[1].annotate("all up: φ = n + 1", (12, np.log(13)), (8.6, 1.2), arrowprops=dict(arrowstyle="->"))
ax[1].legend(loc="upper left")
ax[1].text(0.98, 0.03, "a flip-fair chain would be\nmirror-symmetric here", transform=ax[1].transAxes, ha="right", fontsize=8.5)
fig.tight_layout(); fig.savefig("figures/figA_lopsided_ferromagnet.png", dpi=150); plt.close(fig)

# ---------- Fig B: field zeros, original vs the two fixes ----------
mp.mp.dps = 40
def roots_orig(n, b):
    p = chain(n).astype(float); u = up_of(n); c = np.bincount(u, weights=p**(-b), minlength=n+1); return np.roots(c[::-1])
def roots_sym(n, b):
    p = chain(n).astype(float); u = up_of(n); He = 0.5*(np.log(p)+np.log(p[::-1]))
    c = np.bincount(u, weights=np.exp(-b*He), minlength=n+1); return np.roots(c[::-1])
def roots_dbl(n, b):
    p = chain(n).astype(float); u = up_of(n); c = [mp.fsum([mp.mpf(x) for x in p[u == k]**(-b)]) for k in range(n+1)]
    full = [mp.mpf(0)]*(n+2)
    for k in range(n+1): full[k] += c[k]; full[n+1-k] += c[k]
    return np.array([complex(z) for z in mp.polyroots(full[::-1], maxsteps=500, extraprec=300)])
n = 12; t = np.linspace(0, 2*np.pi, 400)
panels = [("original chain, β = 1.5", roots_orig(n, 1.5), "tab:red"),
          ("flip-averaged chain, β = 1.0", roots_sym(n, 1.0), "tab:red"),
          ("flip-averaged chain, β = 1.5", roots_sym(n, 1.5), "tab:green"),
          ("doubled chain (compass spin), β = 1.5", roots_dbl(n, 1.5), "tab:green")]
fig, ax = plt.subplots(1, 4, figsize=(14, 4.4))
for a, (title, r, c) in zip(ax, panels):
    a.plot(np.cos(t), np.sin(t), color="0.6", lw=1); a.scatter(r.real, r.imag, s=18, color=c, zorder=3)
    dev = np.abs(np.abs(r)-1).max(); a.set_aspect("equal"); lim = max(1.6, 1.1*np.abs(r).max()); a.set_xlim(-lim, lim); a.set_ylim(-lim, lim)
    a.set_title(title, fontsize=9.5); a.set_xlabel("Re w   (w = magnetic field dial)")
    a.text(0.5, 0.03, "largest distance from circle: %.1e" % dev if dev > 1e-9 else "all on the circle", transform=a.transAxes, ha="center", fontsize=8.5)
ax[0].set_ylabel("Im w")
fig.suptitle("B. Field zeros of the n = 12 chain. Grey = Lee–Yang circle. Green panels pass, red fail.", y=1.0)
fig.tight_layout(); fig.savefig("figures/figB_field_zeros_fixes.png", dpi=150); plt.close(fig)

# ---------- Fig C: the temperature where the circle starts, against chain length ----------
fit = json.load(open("fitlim.json"))
ns, bs = load("threshold.json", "threshold_hi.json"); nd, bd = load("doubled_threshold.json")
dd = json.load(open("doubled_direct.json")) if __import__("os").path.exists("doubled_direct.json") else {}
fig, ax = plt.subplots(figsize=(8.4, 4.8))
ax.plot(ns, bs, "o", color="tab:orange", label="flip-averaged chain: circle exact above this β")
ax.plot(nd, bd, "s", color="tab:green", label="doubled chain: proved Lee–Yang above this β (Cohn test)")
if dd:
    k = sorted(int(x) for x in dd if dd[x]["last_fail"] is not None)
    ax.plot(k, [dd[str(x)]["last_fail"] for x in k], "x", color="k", label="doubled chain: last β where zeros left the circle")
xx = np.linspace(10, 60, 200)
for name, col, nn in (("symmetrized", "tab:orange", ns), ("doubled", "tab:green", nd)):
    f = fit[name]["power"]; ax.plot(xx, f["L"]-f["C"]*xx**(-f["a"]), "--", color=col, lw=1)
    ax.text(60.5, f["L"]-f["C"]*60**(-f["a"]), "fit → %.2f?" % f["L"], color=col, fontsize=8.5, va="center")
ax.axhline(2.0, color="tab:purple", lw=1); ax.text(2.5, 2.03, "β = 2: Knauf's phase transition (where Z stops converging)", color="tab:purple", fontsize=8.5)
ax.axhline(1.5, color="tab:red", lw=1); ax.text(27, 1.53, "β = 3/2: Knauf's Riemann line", color="tab:red", fontsize=8.5)
ax.set_xlabel("chain length n"); ax.set_ylabel("inverse temperature β"); ax.set_xlim(2, 70); ax.set_ylim(0, 2.3)
ax.set_title("C. The circle holds only on the cold side, and the edge climbs past the Riemann line")
ax.legend(loc="lower right", fontsize=8.5)
fig.tight_layout(); fig.savefig("figures/figC_threshold_vs_length.png", dpi=150); plt.close(fig)

# ---------- Fig D: does the fixed chain still turn into zeta? ----------
beta = 3.0; target = float(mp.zeta(beta-1)/mp.zeta(beta)); N = list(range(4, 23)); zo, zs = [], []
for n in N:
    p = chain(n).astype(float); zo.append((p**(-beta)).sum()); zs.append(((p*p[::-1])**(-beta/2)).sum())
fig, ax = plt.subplots(figsize=(8, 4.4))
ax.plot(N, zo, "o-", color="tab:blue", label="original chain  Z")
ax.plot(N, np.array(zo), "s", mfc="none", ms=10, color="tab:green", label="doubled chain  Z/2  (equals Z exactly)")
ax.plot(N, zs, "^-", color="tab:orange", label="flip-averaged chain")
ax.axhline(target, color="k", lw=1, ls=":"); ax.text(15, target-0.07, "ζ(β−1)/ζ(β) = %.4f" % target, fontsize=9)
ax.set_xlabel("chain length n"); ax.set_ylabel("partition function at β = 3")
ax.set_title("D. Flip-averaging drains the ζ limit away; the compass spin keeps it")
ax.legend(loc="center right")
fig.tight_layout(); fig.savefig("figures/figD_zeta_kept_or_lost.png", dpi=150); plt.close(fig)
print("figs done")
