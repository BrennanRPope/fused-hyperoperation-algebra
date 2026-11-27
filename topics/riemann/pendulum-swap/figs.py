"""Pendulum-swap figures. DRAFT."""
import numpy as np, json, sys, os, mpmath as mp, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, "../knauf-field"); from knauf import chain
from frames import walsh
plt.rcParams.update({"font.size": 10}); os.makedirs("figures", exist_ok=True)
bf = json.load(open("bestflip.json")); fr = json.load(open("frames.json"))

def figA():
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.8))
    a = ax[0]; n = 6
    for i in range(n):
        a.add_patch(plt.Circle((i, 1.0), 0.2, color="tab:blue")); a.text(i, 1.0, "s%d" % (i+1), ha="center", va="center", color="w", fontsize=8)
    for i in range(n-1):
        a.plot([i+0.2, i+0.8], [1.0, 1.0], color="k", lw=2.5)
        xb = i + 0.5; a.add_patch(plt.Circle((xb, 0.0), 0.2, color="tab:orange")); a.text(xb, 0.0, "b%d" % (i+1), ha="center", va="center", color="w", fontsize=8)
        a.annotate("", (xb, 0.24), (xb, 0.92), arrowprops=dict(arrowstyle="->", color="0.5"))
    a.add_patch(plt.Circle((n-1, 0.0), 0.2, color="tab:orange")); a.text(n-1, 0.0, "b6", ha="center", va="center", color="w", fontsize=8)
    a.text(-0.3, 1.38, "spin frame: energy sits in LINKS (black bars) between spins", fontsize=8.5)
    a.text(-0.3, -0.5, "dual frame: bond spin b_i = s_i·s_(i+1); each link is now a FIELD on one b", fontsize=8.5)
    a.text(2.5, -1.15, "pendulum swap on spins (Kramers–Wannier):  link (hot) ↔ field (cold),  e^(−β*) = tanh(β/2)",
           ha="center", fontsize=8.5, bbox=dict(fc="w", ec="0.7"))
    a.set_aspect("equal"); a.set_xlim(-0.5, n); a.set_ylim(-1.5, 1.7); a.axis("off"); a.set_title("A1. The swap rewrites links as fields")
    b = ax[1]; n = 12; H = np.log(chain(n).astype(float)); He = 0.5*(H+H[::-1]); J = walsh(He)
    sz = np.array([bin(k).count("1") for k in range(2**n)]); P = np.zeros((n, n))
    for k in np.nonzero(sz == 2)[0]:
        i, j = [q for q in range(n) if k >> q & 1]; P[i, j] = P[j, i] = -J[k]
    tot = (J[1:]**2).sum(); pairpow = (J[sz == 2]**2).sum()/tot
    im = b.imshow(P, cmap="Reds"); plt.colorbar(im, ax=b, label="pair link strength (ferromagnetic)")
    b.set_xlabel("spin j"); b.set_ylabel("spin i")
    b.set_title("A2. The chain's even half, n = 12: %.1f%% is pair links,\nstrongest between near neighbours" % (100*pairpow), fontsize=10)
    fig.tight_layout(); fig.savefig("figures/figA_swap_on_spins.png", dpi=150); plt.close(fig)

def figB():
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
    a = ax[0]; ns = [8, 10, 12, 14, 16]
    a.plot(ns, [bf[str(n)]["all"] for n in ns], "o-", color="tab:red", label="flip every spin (the old frame)")
    a.plot(ns, [bf[str(n)]["alternating"] for n in ns], "s-", color="tab:green", ms=9, label="flip every other spin = the dual frame")
    a.plot(ns, [bf[str(n)]["best"] for n in ns], "k+", ms=14, mew=2, label="best of all 2ⁿ−1 flip choices")
    a.axhline(1, color="0.6", ls=":"); a.text(8, 0.92, "100% = whole chain in the safe class (would carry the line)", fontsize=8.5, color="0.4")
    a.set_ylim(0, 1.05); a.set_xlabel("chain length n"); a.set_ylabel("share of energy in the odd (safe) class")
    a.legend(fontsize=8.5, loc="lower left"); a.set_title("B1. The dual frame puts 81% in the safe class, and it is the best")
    b = ax[1]; n = 12; H = np.log(chain(n).astype(float)); Jw = walsh(H); Jw[0] = 0; Pw = Jw**2; tt = Pw.sum()
    odd = 0.5*(tt - walsh(Pw)*len(Pw))/tt
    b.hist(odd[1:], bins=60, color="0.6"); alt = sum(1 << i for i in range(n-1, -1, -2))
    b.axvline(odd[alt], color="tab:green", lw=2); b.text(odd[alt]-0.01, b.get_ylim()[1]*0.8, "dual frame\n%.1f%%" % (100*odd[alt]), color="tab:green", ha="right", fontsize=9)
    b.axvline(odd[2**n-1], color="tab:red", lw=2); b.text(odd[2**n-1]+0.01, b.get_ylim()[1]*0.6, "flip all\n%.1f%%" % (100*odd[2**n-1]), color="tab:red", fontsize=9)
    b.set_xlabel("odd (safe) share"); b.set_ylabel("number of flip choices (n = 12)")
    b.set_title("B2. All 4095 flip choices; the last ~19% is frustrated\n(link triangles: no flip can make every link odd)", fontsize=10)
    fig.tight_layout(); fig.savefig("figures/figB_best_frame.png", dpi=150); plt.close(fig)

def figC():
    lc = json.load(open("linecount.json")); fig, ax = plt.subplots(figsize=(8.5, 4.4))
    for key, c, m in (("spin frame", "tab:blue", "o"), ("dual frame", "tab:orange", "s")):
        for n, ls in ((8, ":"), (10, "--"), (12, "-")):
            k = "%s n=%d" % (key, n)
            if k not in lc: continue
            off = np.zeros(91)
            for d, on, tot in lc[k]: off[d] = tot - on
            ax.plot(range(91), off + (0.04 if key == "dual frame" else 0), ls, color=c, marker=m, ms=3, markevery=5, label="%s, n = %d" % (key, n))
    ax.set_xlabel("rewind angle θ (degrees) — odd half only"); ax.set_ylabel("zeros off the line (Im ≤ 36)")
    ax.set_title("C. Rewinding the odd half: the spin frame never leaves the line;\nthe dual frame is clean at n = 8, but mirror pairs leave for a while (n = 10: 5°–23°, n = 12: 4°–27°)", fontsize=10)
    ax.legend(fontsize=8); fig.tight_layout(); fig.savefig("figures/figC_line_counts.png", dpi=150); plt.close(fig)

def g(x): return 1/(4*np.sinh(x/2)**2) - 1/x**2
def figD():
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.4)); x = np.linspace(0.25, 8, 400)
    a = ax[0]; a.plot(x, 1/x**2, color="tab:blue", label="1/x²  (classical oscillator)")
    a.plot(x, 1/(4*np.sinh(x/2)**2), color="tab:orange", label="1/(4 sinh²(x/2)) = Σ n e^(−nx)  (quantum)")
    a.fill_between(x, 1/(4*np.sinh(x/2)**2), 1/x**2, color="tab:purple", alpha=0.3, label="gap: total area = 1/2 = −ζ(0)")
    a.set_ylim(0, 4); a.set_xlabel("x = 1/temperature"); a.legend(fontsize=8); a.set_title("D1. The area finding (known: ζ(0) = −1/2)")
    b = ax[1]; X = np.linspace(0.05, 10, 400)
    b.plot(X, 1/(2*np.tanh(X/2)), color="tab:orange", label="quantum energy  ½coth(x/2)")
    b.plot(X, 1/X, color="tab:blue", label="classical energy  1/x")
    b.plot(X, 1/(2*np.tanh(X/2)) - 1/X, color="tab:purple", lw=2, label="difference → ½ as it cools")
    b.axhline(0.5, color="tab:purple", ls=":"); b.set_ylim(0, 3); b.set_xlabel("x = 1/temperature")
    b.legend(fontsize=8); b.set_title("D2. What the 1/2 is: the zero-point energy ½ħω\n(area = −[energy gap at cold end])", fontsize=10)
    c = ax[2]; xx = np.linspace(0.2, 10, 300)
    for K, col in ((1, "tab:red"), (3, "tab:green"), (3000, "k")):
        S = sum(2*(xx**2 - 4*np.pi**2*k**2)/(xx**2 + 4*np.pi**2*k**2)**2 for k in range(1, K+1))
        c.plot(xx, S, color=col, lw=1.2 if K < 50 else 3.5, alpha=0.9 if K < 50 else 0.3, label="first %d frequency pairs" % K)
    c.plot(xx, g(xx), "--", color="tab:purple", label="the gap 1/(4sinh²) − 1/x²")
    c.set_xlabel("x"); c.legend(fontsize=8)
    c.set_title("D3. The swap (Poisson): the gap is the sum over\nfrequencies 2πk, k ≠ 0; 1/x² is the k = 0 mode", fontsize=10)
    fig.tight_layout(); fig.savefig("figures/figD_area_half.png", dpi=150); plt.close(fig)

def figE():
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
    a = ax[0]; a.axvspan(0, 2, color="tab:purple", alpha=0.12)
    a.axvline(1.5, color="tab:green", lw=2.5); a.plot(1, 0, "o", color="tab:purple", ms=10); a.plot(2, 0, "x", color="k", ms=12, mew=3)
    a.text(0.15, 1.5, "s = 1: the area\n(value −1/2)", color="tab:purple", fontsize=9)
    a.text(1.55, 18, "Re s = 3/2:\nthe critical line\nof ζ(s−1)\n(Knauf's numerator)", color="tab:green", fontsize=9)
    a.text(2.05, -1.5, "pole of ζ(s−1)", fontsize=8.5); a.text(0.1, 27, "strip where the gap curve's\nMellin transform converges", color="tab:purple", fontsize=8.5)
    for t in (14.13, 21.02, 25.01): a.plot(1.5, t, "o", color="tab:green", ms=6)
    a.set_xlim(-0.3, 3); a.set_ylim(-3, 32); a.set_xlabel("Re s"); a.set_ylabel("Im s")
    a.set_title("E1. Same curve, two tilts: the area is a value\nat s = 1, the line is the spectrum at Re s = 3/2", fontsize=10)
    b = ax[1]; u = np.linspace(-6, 4, 600); K = np.exp(1.5*u)*g(np.exp(u))
    def Phi(u): return sum((2*np.pi**2*n**4*np.exp(4.5*u) - 3*np.pi*n**2*np.exp(2.5*u))*np.exp(-np.pi*n**2*np.exp(2*u)) for n in range(1, 8))
    up = np.linspace(0, 3, 300); uu = np.concatenate([-up[::-1], up[1:]]); Pp = Phi(up); Pv = np.concatenate([Pp[::-1], Pp[1:]])  # Phi is even (theta swap); series converges fast for u >= 0
    b.plot(u, K/abs(K).max(), color="tab:purple", lw=2, label="area kernel  e^(3u/2) g(e^u): one sign, LOPSIDED")
    b.plot(uu, Pv/Pv.max(), color="tab:green", lw=2, label="Riemann's Φ(u) (after the theta swap): EVEN")
    b.axhline(0, color="0.6", lw=0.8); b.axvline(0, color="0.6", lw=0.8, ls=":"); b.set_xlabel("u = log x"); b.legend(fontsize=8)
    b.set_title("E2. The pendulum swap is what makes the kernel even", fontsize=10)
    c = ax[2]; t = np.linspace(1, 40, 1500)
    FK = np.array([abs(complex(mp.gamma(mp.mpc(1.5, tt))*mp.zeta(mp.mpc(0.5, tt)))) for tt in t])
    xi = np.array([float(mp.re(mp.mpf(0.5)*(mp.mpc(0.5, tt))*(mp.mpc(-0.5, tt))*mp.pi**(-mp.mpc(0.5, tt)/2)*mp.gamma(mp.mpc(0.25, tt/2))*mp.zeta(mp.mpc(0.5, tt)))) for tt in t])
    c.semilogy(t, FK, color="tab:purple", label="|Fourier transform of the area kernel|")
    c.semilogy(t, np.abs(xi), color="tab:green", alpha=0.7, label="|ξ(½+it)| from the even kernel")
    for z in (14.1347, 21.0220, 25.0109, 30.4249, 32.9351, 37.5862): c.axvline(z, color="0.7", ls=":")
    c.set_xlabel("height t on the line"); c.legend(fontsize=8)
    c.set_title("E3. Both transforms vanish exactly at ζ's zeros\n(dotted); neither says why they are on the line", fontsize=10)
    fig.tight_layout(); fig.savefig("figures/figE_area_on_line.png", dpi=150); plt.close(fig)

if __name__ == "__main__":
    which = sys.argv[1:] or ["A", "B", "C", "D", "E"]
    for w in which: globals()["fig" + w](); print("fig", w, "done", flush=True)
