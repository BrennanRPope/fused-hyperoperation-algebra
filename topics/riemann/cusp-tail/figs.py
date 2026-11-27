"""Figures for cusp-tail. DRAFT."""
import numpy as np, json, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 10})
d = json.load(open("demod.json")); c = json.load(open("combo.json")); w = json.load(open("../length-smoothing/wide.json")); nl = json.load(open("null.json"))
ns = np.array(d["ns"]); gam = d["gam"]
# A: the swing is the cusp
fig, ax = plt.subplots(1, 3, figsize=(15, 4.8))
j = 3; g = gam[j]; E = np.array([complex(*e) for e in d["E"][j]]); sel = ns >= 14; D = E*ns**(1.5+1j*g-1)
for a, Z, ttl in ((ax[0], E, "A1. Chain error at ζ's 4th zero (γ = %.1f), lengths 14 to 33:\nit spins round 0" % g),
                  (ax[1], D, "A2. Same points times n^(ρ−1), the cusp's own spin:\nthey stop spinning and sit in one place")):
    z = Z[sel]; sc = a.scatter(z.real, z.imag, c=ns[sel], cmap="viridis", s=36, zorder=3); a.plot(z.real, z.imag, "-", color="0.8", lw=0.8, zorder=2)
    a.axhline(0, color="0.6", lw=0.6); a.axvline(0, color="0.6", lw=0.6); a.set_aspect("equal", "datalim"); a.set_title(ttl, fontsize=10)
    a.set_xlabel("real part"); a.set_ylabel("imaginary part"); plt.colorbar(sc, ax=a, fraction=0.046, label="chain length n")
b = ax[2]; x = np.arange(len(gam)); b.bar(x-0.2, d["raw"], 0.4, color="tab:red", label="raw error")
b.bar(x+0.2, d["demod"], 0.4, color="tab:blue", label="after removing the cusp's spin")
b.set_xticks(x); b.set_xticklabels(["%.1f" % v for v in gam]); b.set_xlabel("height γ of ζ's zero"); b.set_ylabel("how still the error's direction is (1 = fixed)")
b.set_ylim(0, 1.3); b.legend(fontsize=8.5, loc="upper center", ncol=2); b.set_title("A3. Lengths 18 to 33: the cusp's spin is\nmost of the swing for the first four zeros", fontsize=10)
fig.tight_layout(); fig.savefig("figures/figA_swing_is_cusp.png", dpi=150); plt.close(fig)
# B: removing the cusp tail
fig, ax = plt.subplots(figsize=(8.5, 4.6)); Ns = range(20, 34)
rms = lambda rows: np.sqrt(np.mean([r[0]**2 for r in rows]))
ax.plot(list(Ns), [rms(w["rows"]["single N=%d" % N]) for N in Ns], "o-", color="tab:red", ms=4, label="one chain length")
ax.plot(list(Ns), [rms(w["rows"]["k=N/2 N=%d" % N]) for N in Ns], "o-", color="tab:purple", ms=4, label="averaged over half the chain (length-smoothing)")
ax.plot(list(Ns), [rms(c["rows"]["cusp half N=%d" % N][:5]) for N in Ns], "o-", color="tab:blue", ms=4, label="cusp tail fitted out over half the chain")
ax.set_yscale("log"); ax.set_xlabel("longest chain length N"); ax.set_ylabel("sideways miss (RMS, first 5 zeros)"); ax.legend(fontsize=8.5)
ax.set_title("B. With the cusp tail removed, the miss drops to about 0.0005 by length 29 and stays there,\n2 to 3 times below the averaged chain on lengths 29 to 32 (at 33 the two meet)", fontsize=10)
fig.tight_layout(); fig.savefig("figures/figB_cusp_removed.png", dpi=150); plt.close(fig)
# C: null test
fig, ax = plt.subplots(figsize=(8.5, 4.2)); names = list(nl.keys()); vals = [np.mean(nl[k]) for k in names]
cols = ["tab:blue" if k.startswith("right") else "0.6" for k in names]
ax.barh(range(len(names)), vals, color=cols); ax.set_yticks(range(len(names))); ax.set_yticklabels(names, fontsize=8.5); ax.invert_yaxis()
ax.set_xscale("log"); ax.set_xlabel("mean sideways miss, lengths 27 to 33 (8 zeros)")
ax.set_title("C. The fit only works with the cusp's shape; wrong shapes do 1.4 to 20 times worse\n(a phase off by n^(0.5i) is too small to tell apart at these lengths)", fontsize=10)
fig.tight_layout(); fig.savefig("figures/figC_null_test.png", dpi=150); plt.close(fig); print("done")
