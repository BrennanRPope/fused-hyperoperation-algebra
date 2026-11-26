"""Pictures for the flux-tube dial DRAFT. Run: python3 -I figures.py -> figures/ (uses data/ from make_data.py)"""
import json, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 11, "axes.titlesize": 12, "figure.dpi": 130})


def save(fig, name):
    fig.tight_layout(); fig.savefig(f"figures/{name}.png", bbox_inches="tight"); plt.close(fig)


# g1: flux maps
fig, axs = plt.subplots(1, 4, figsize=(15, 4.2))
N, box = 160, 20.0
x = (np.arange(N) - N / 2) * box / N
sel = (np.abs(x) < 5)
for ax, sv in zip(axs, ("2.0", "1.5", "1.25", "1.1")):
    Dm = np.load(f"data/Dm_s{sv}.npy")
    im = ax.imshow((Dm[np.ix_(sel, sel)] / Dm[np.ix_(sel, sel)].max()).T ** 0.5, origin="lower", cmap="magma",
                   extent=[x[sel][0], x[sel][-1], x[sel][0], x[sel][-1]])
    ax.plot([-2, 2], [0, 0], "c+", ms=12, mew=2)
    lab = {"2.0": "s = 2: Coulomb, flux spreads", "1.5": "s = 1.5", "1.25": "s = 1.25", "1.1": "s = 1.1: nearly a tube"}[sv]
    ax.set_title(lab); ax.set_xticks([]); ax.set_yticks([])
fig.suptitle("g1. The flux exponent s: same two charges (+ at left, − at right, 4 apart); brightness = flux density. "
             "As s falls toward 1 the flux pulls into a tube along the line between them")
save(fig, "g1_flux_maps")

# g2: effective dimension vs s and confinement threshold
s = np.linspace(0.9, 2.2, 300)
fig, ax = plt.subplots(figsize=(8.5, 5))
for D, col in ((3, "#d62728"), (2, "#1f77b4")):
    ax.plot(s, 1 + (D - 1) * (s - 1), color=col, lw=2.5, label=f"D = {D} space dimensions: D_eff = 1 + {D-1}(s − 1)")
    ax.axvline(D / (D - 1), color=col, ls=":", lw=1.5)
ax.axhspan(0, 2, color="#2ca02c", alpha=0.08)
ax.text(1.02, 0.5, "D_eff < 2: potential grows with distance (confines)", color="#2ca02c", fontsize=10)
ax.text(1.52, 3.3, "3D threshold s = 3/2\n(Pagels–Tomboulis δ = 1/4)", color="#d62728", fontsize=9)
ax.plot([1], [1], "ko", ms=9); ax.text(1.06, 0.72, "s = 1: straight tube,\nV = flux × distance", fontsize=10)
ax.plot([2], [3], "o", color="#d62728", ms=9); ax.text(1.75, 3.05, "s = 2: Maxwell", fontsize=10)
ax.set_xlabel("flux exponent s  (field energy density ∝ |D|^s)"); ax.set_ylabel("effective dimension of the flux")
ax.set_ylim(0, 4.2); ax.legend(loc="upper left", fontsize=9)
ax.set_title("g2. The flux exponent works as a dimension dial for the flux alone")
save(fig, "g2_effective_dimension")

# g3: bundle costs and junctions
fig, axs = plt.subplots(1, 2, figsize=(12, 4.6))
P = np.linspace(0, 3, 200)
for sv, col in ((0.6, "#2ca02c"), (1.0, "k"), (1.5, "#d62728")):
    axs[0].plot(P, P ** sv, color=col, lw=2.5,
                label=f"s = {sv}: " + ("merging saves (bunch)" if sv < 1 else ("neutral" if sv == 1 else "splitting saves (spread)")))
axs[0].set_xlabel("flux in one bundle Φ"); axs[0].set_ylabel("cost per length Φ^s"); axs[0].legend(fontsize=9)
axs[0].set_title("A: bundle cost; bending down means bundles want to merge")
ax = axs[1]
for i, sv in enumerate((1.0, 0.75, 0.5)):
    x0 = i * 3.2
    half = math.acos(2 ** (sv - 1))
    jx = x0 + 1.4
    for sgn in (1, -1):
        ax.plot([jx - 1.2 * math.cos(half), jx], [sgn * (0.02 + 1.2 * math.sin(half)), 0], color="#2ca02c", lw=2)
    ax.plot([jx, jx + 1.3], [0, 0], color="#2ca02c", lw=4)
    ax.text(jx - 0.3, -1.5, f"s = {sv}\nopening {2*math.degrees(half):.0f}°", ha="center", fontsize=10)
ax.set_xlim(-0.5, 9.5); ax.set_ylim(-2, 1.6); ax.set_aspect("equal"); ax.axis("off")
ax.set_title("B: two unit bundles join into one of 2\n(Gilbert junction: cos(half angle) = 2^(s−1))")
fig.suptitle("g3. s below 1 makes flux lines merge into one thick bundle: the switch 'keep the flux together' as a dial")
save(fig, "g3_bundles")

# g4: vortex tensions per unit
V = json.load(open("data/vortex_E.json"))
n = np.arange(1, 9)
fig, axs = plt.subplots(1, 2, figsize=(12, 4.6))
cm = plt.get_cmap("coolwarm")
lams = sorted(V, key=float)
for i, lv in enumerate(lams):
    Es = np.array(V[lv]); col = cm(i / (len(lams) - 1))
    kap = math.sqrt(float(lv) / 2)
    axs[0].plot(n, Es / n, "o-", color=col, lw=2, label=f"λ = {lv} (κ = {kap:.2f})")
axs[0].set_xlabel("flux units n in one tube"); axs[0].set_ylabel("tension per unit E(n)/(n π)")
axs[0].set_title("A: abelian Higgs tubes: flat at the border, falling (bunch) in type I, rising in type II")
axs[0].legend(fontsize=8, ncol=2, loc="upper left"); axs[0].set_ylim(0.35, 3.3)
F = json.load(open("results.json"))["vortex_fits"]
lv = np.array([float(k) for k in lams])
axs[1].semilogy(lv, [max(F[k]["err_lp"], 1e-6) * 100 for k in lams], "o-", color="#2ca02c", lw=2.5,
                label="ℓp sum (tension ∝ n^s), s from n = 1, 2")
axs[1].semilogy(lv, [max(F[k]["err_ok"], 1e-6) * 100 for k in lams], "s-", color="#9467bd", lw=2.5,
                label="our law ∘k, k from n = 1, 2")
axs[1].set_xscale("log"); axs[1].axvline(1, color="k", ls=":")
axs[1].set_xlabel("λ (border at 1)"); axs[1].set_ylabel("worst error, n = 1 to 8 (%)"); axs[1].legend(fontsize=9)
axs[1].set_title("B: which composition law fits real tubes")
axs[1].text(1.08, 1.3e-4, "both exact at the border\n(drawn at the floor)", fontsize=8)
fig.suptitle("g4. Real flux tubes: the border λ = 1 is plain addition; near it tensions combine like an ℓp sum, not like ∘k")
save(fig, "g4_vortex_tensions")

# g5: s reading vs kappa
fig, ax = plt.subplots(figsize=(8.5, 4.8))
kap = np.sqrt(lv / 2)
ax.plot(kap, [F[k]["s"] for k in lams], "o-", color="k", lw=2.5)
ax.axvline(1 / math.sqrt(2), color="k", ls=":"); ax.axhline(1, color="k", ls=":")
ax.axvspan(1.49, 1.70, color="#ff7f0e", alpha=0.25)
ax.text(1.5, 0.86, "lattice estimates for\nthe QCD vacuum (SU(2);\nconventions differ;\nnot settled)", fontsize=9)
ax.text(0.2, 1.02, "type I: s < 1, tubes merge", fontsize=10, color="#1f77b4")
ax.text(0.8, 0.82, "border κ = 1/√2: s = 1 exactly", fontsize=10)
ax.text(1.0, 1.17, "type II: s > 1, tubes split", fontsize=10, color="#d62728")
ax.set_xlabel("Ginzburg–Landau κ"); ax.set_ylabel("flux exponent read from tubes  s = log₂(E₂/E₁)")
ax.set_title("g5. The flux exponent as a reading of the medium: one number for 'bunch or spread'")
save(fig, "g5_s_reading")

# g6: summary ladder of s
fig, ax = plt.subplots(figsize=(11, 3.6))
ax.set_xlim(0.3, 2.3); ax.set_ylim(0, 1); ax.set_yticks([])
bands = [(0.3, 1.0, "#2ca02c", "s < 1: bundles merge;\njunctions (Y shapes);\nbranched transport"),
         (1.0, 1.5, "#bcbd22", "1 < s < 3/2: confines,\nbut the flux bundle\nwidens with distance"),
         (1.5, 2.3, "#1f77b4", "s > 3/2: no confinement;\ns = 2 is Maxwell")]
for x0, x1, c, t in bands:
    ax.axvspan(x0, x1, color=c, alpha=0.25); ax.text((x0 + x1) / 2, 0.55, t, ha="center", fontsize=10)
ax.axvline(1, color="k", lw=3); ax.text(1.0, 0.1, "s = 1: one straight tube,\nV = flux × distance", ha="center", fontsize=10,
                                         bbox=dict(fc="white", ec="k"))
ax.set_xlabel("flux exponent s (3D)")
ax.set_title("g6. The dial at a glance. A tube of fixed width also needs a length (the screening mass of the first page)")
save(fig, "g6_ladder")
print("ok")
