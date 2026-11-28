"""Figures for k-dial-cusp. DRAFT."""
import numpy as np, json, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from rewind import PH
plt.rcParams.update({"font.size": 10}); r = json.load(open("rewind.json")); cmp = json.load(open("../cusp-size/compare.json"))
fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
a = ax[0]; ks = np.array(r["ks"])
for key, col in (("5.0000", "tab:green"), ("14.1347", "tab:blue"), ("25.0109", "tab:purple"), ("40.0000", "tab:red")):
    a.loglog(ks, [v[1] for v in r["rows"][key]], "o-", color=col, ms=4, label="height %.1f%s" % (float(key), " (a ζ zero)" if key in ("14.1347", "25.0109") else ""))
a.loglog(ks, 3*ks, ":", color="0.4", label="slope 1: miss shrinks in step with k")
a.axhline(np.median(cmp["cusp2"]), color="tab:orange", ls="--", lw=1.2, label="cusp fit from chains to length 33 (median)")
a.invert_xaxis(); a.set_xlabel("k (dialled down toward 0 →)"); a.set_ylabel("distance from ζ(β−1)/ζ(β)"); a.legend(fontsize=8)
a.set_title("A1. Dialling k back to 0 converges, but only in step with k:\nslower than fitting out the cusp", fontsize=10)
b = ax[1]; d = np.load("/tmp/knauf_levels/levels_33.npz"); q, c = d["q"], d["c"]; M = 400000
cn = np.zeros(M); sel = q < M; cn[q[sel]] = c[sel]; share = cn[1:]/PH[1:M]; m = np.arange(1, M)
edges = np.unique(np.logspace(0, np.log10(M-1), 60).astype(int)); mid = []; avg = []
for lo, hi in zip(edges[:-1], edges[1:]): mid.append(np.sqrt(lo*hi)); avg.append(share[lo-1:hi-1].mean())
b.semilogx(mid, avg, "o-", color="tab:red", ms=3, lw=1, label="chain of length 33: share of each denominator kept (binned)")
for k, col in ((0.001, "tab:blue"), (0.0001, "tab:cyan")):
    L = np.log1p(k); x = L*m; b.semilogx(m, (x/np.expm1(x))**1.5, color=col, lw=2, label="k dial at k = %g: smooth weight (x/(eˣ−1))^β, x = m·log(1+k)" % k)
b.set_xlabel("denominator m"); b.set_ylabel("weight on φ(m)·m^(−β)"); b.set_ylim(0, 1.08); b.legend(fontsize=8, loc="lower left")
b.set_title("A2. Why: the chain keeps fractions by depth, a weight that thins out slowly and raggedly\n(the cusp); the k dial keeps them by size, a smooth Bose-shaped weight", fontsize=10)
fig.tight_layout(); fig.savefig("figures/figA_k_dial_rewind.png", dpi=150); plt.close(fig); print("done")
