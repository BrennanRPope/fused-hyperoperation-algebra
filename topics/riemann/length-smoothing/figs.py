"""Length-smoothing figures. DRAFT."""
import numpy as np, json, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 10}); d = json.load(open("windows.json")); R = d["rows"]; Ns = range(16, 25)
ref = json.load(open("../knauf-field/refined.json"))
def rms(kind, k, N): return np.sqrt(np.mean([r[0]**2 for r in R["%s k=%d N=%d" % (kind, k, N)]]))
fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.8))
a = ax[0]
for key, c, lab in (("tri k=0", "tab:red", "one chain length (the raw wobble)"), ("tri k=8", "tab:green", "averaged over 9 lengths (triangle weights)")):
    pts = np.array([r for N in Ns for r in R["%s N=%d" % (key, N)]])
    a.plot(pts[:, 0], pts[:, 1], "o", color=c, ms=5, alpha=0.75, label=lab)
a.axvline(0, color="k", lw=0.8); a.axhline(0, color="k", lw=0.8)
a.set_xlabel("Re β − 3/2  (sideways miss from the line)"); a.set_ylabel("Im β − γ  (miss in height)")
a.set_title("A1. Where the chain's zeros land, relative to ζ's\n(first 5 zeros, chain lengths 16 to 24)", fontsize=10); a.legend(fontsize=8.5)
a.set_xlim(-0.06, 0.06); a.set_ylim(-0.06, 0.06); a.set_aspect("equal")
b = ax[1]
for k, c in ((0, "tab:red"), (4, "tab:orange"), (8, "tab:green"), (12, "tab:blue")):
    b.plot(list(Ns), [rms("tri", k, N) for N in Ns], "o-", color=c, label=("1 length (no averaging)" if k == 0 else "window of %d lengths" % (k+1)))
b.set_xlabel("longest chain length N"); b.set_ylabel("sideways wobble (root-mean-square, 5 zeros)"); b.set_yscale("log")
b.legend(fontsize=8.5); b.set_title("A2. Averaging over lengths cuts the wobble about 4×;\nwider windows keep improving as N grows", fontsize=10)
fig.tight_layout(); fig.savefig("figures/figA_wobble_shrinks.png", dpi=150); plt.close(fig)

fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.6))
a = ax[0]; ks = list(range(0, 13, 2))
for kind, c in (("flat", "tab:purple"), ("tri", "tab:green")):
    a.plot([k+1 for k in ks], [np.mean([rms(kind, k, N) for N in Ns]) for k in ks], "o-", color=c, label="%s weights" % ("equal" if kind == "flat" else "triangle"))
a.set_xlabel("window size (number of chain lengths averaged)"); a.set_ylabel("mean wobble, N = 16..24")
a.text(1.5, 0.017, "too narrow: the length-to-length\nswing is not cancelled", fontsize=8.5); a.text(8, 0.017, "too wide: short chains\ndrag the zeros off", fontsize=8.5)
a.legend(fontsize=8.5); a.set_title("B1. The best window is about 7 to 9 lengths (for N ≤ 24)", fontsize=10)
b = ax[1]; rows = ref["rows"]; ns = sorted(int(n) for n in rows)
for j, c in enumerate(("tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple")):
    b.plot(ns, [rows[str(n)][j][0]-1.5 for n in ns], "o-", color=c, ms=4, label="zero near height %.2f" % ref["gam"][j])
b.axhline(0, color="k", lw=0.8); b.set_xlabel("chain length n"); b.set_ylabel("Re β − 3/2 of one chain's zero")
b.legend(fontsize=7.5, ncol=2); b.set_title("B2. Why averaging helps: each zero swings to both sides\nof the line from one length to the next", fontsize=10)
fig.tight_layout(); fig.savefig("figures/figB_window_and_swing.png", dpi=150); plt.close(fig)
print("done")
