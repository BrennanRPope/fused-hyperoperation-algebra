"""Follow-up figures: long chains and the law of the swing. DRAFT."""
import numpy as np, json, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 10}); w = json.load(open("wide.json")); s = json.load(open("swing.json"))
Ns = range(16, 34)
def rms(name, N): return np.sqrt(np.mean([r[0]**2 for r in w["rows"]["%s N=%d" % (name, N)]]))
fig, ax = plt.subplots(figsize=(8.5, 4.6))
for name, c, lab in (("single", "tab:red", "one chain length"), ("k=8", "tab:green", "fixed window, 9 lengths"),
                     ("k=12", "tab:blue", "fixed window, 13 lengths"), ("k=N/2", "tab:purple", "growing window, half the chain length")):
    ax.plot(list(Ns), [rms(name, N) for N in Ns], "o-", color=c, ms=4, label=lab)
ax.set_yscale("log"); ax.set_xlabel("longest chain length N (2^N states; 33 means 8.6 billion)"); ax.set_ylabel("sideways wobble (RMS, first 5 zeros)")
ax.legend(fontsize=8.5); ax.set_title("C. Out to length 33: averaging takes the wobble from about 0.01 to about 0.001", fontsize=10)
fig.tight_layout(); fig.savefig("figures/figC_long_chains.png", dpi=150); plt.close(fig)
fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.6))
a = ax[0]; g = np.array(s["gam"]); t = np.array(s["turn"])
a.plot(g, -t, "o", color="tab:blue", ms=7, label="measured turn per length (lengths 16 to 30)")
gg = np.linspace(10, 45, 50); a.plot(gg, -(s["fit"][0]*gg + s["fit"][1]), "--", color="0.4", label="straight line: %.3f × height %+.2f" % (-s["fit"][0], -s["fit"][1]))
a.axhline(np.pi, color="tab:red", ls=":", lw=1); a.text(11, np.pi+0.05, "half a turn (swings faster than this look slower: aliasing)", fontsize=8, color="tab:red")
a.set_xlabel("height γ of ζ's zero"); a.set_ylabel("angle the chain's error turns per extra length (rad)"); a.set_ylim(0, 3.6); a.legend(fontsize=8.5, loc="lower right")
a.set_title("D1. The swing: the error turns by a steady angle each length,\nfaster for higher zeros", fontsize=10)
b = ax[1]; per = 2*np.pi/np.abs(t[:5]); best = s["best_window_per_zero"]
b.plot(per, best, "o", color="tab:green", ms=9)
for p, bb, gm in zip(per, best, g[:5]): b.annotate("γ = %.1f" % gm, (p, bb), xytext=(5, 4), textcoords="offset points", fontsize=8.5)
b.plot([2, 12], [2, 12], ":", color="0.5", label="window = one swing"); b.plot([2, 6], [4, 12], ":", color="0.75", label="window = two swings")
b.set_xlabel("swing period (lengths per full turn)"); b.set_ylabel("best averaging window (lengths)"); b.legend(fontsize=8.5)
b.set_xlim(2, 12); b.set_ylim(2, 12); b.set_title("D2. The best window roughly tracks the swing period\n(about one to two swings; windows tried: 1 to 13 lengths)", fontsize=10)
fig.tight_layout(); fig.savefig("figures/figD_swing_law.png", dpi=150); plt.close(fig); print("done")
