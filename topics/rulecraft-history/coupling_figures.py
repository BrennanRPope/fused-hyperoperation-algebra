"""Pictures for the coupling section of rulecraft-history (DRAFT). Every shape is computed from the formulas in
coupling_checks.py / coupling_polish.py; nothing fitted. Run: python3 coupling_figures.py -> figures/fig9..fig11."""
import itertools, math, random, numpy as np, matplotlib, networkx as nx
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
GREEN, ORANGE, BLUE, GREY, RED, PURPLE = "#2a9d55", "#e07b1f", "#2f6fb5", "#9a9a9a", "#c8382f", "#7a4fb5"
def save(fig, name):
    fig.savefig(f"figures/{name}.png", dpi=150, bbox_inches="tight"); plt.close(fig)
def chrom_at(G, q):
    G = nx.Graph(G)
    if G.number_of_edges() == 0: return q**G.number_of_nodes()
    u, v = next(iter(G.edges())); G1 = G.copy(); G1.remove_edge(u, v)
    G2 = nx.Graph(nx.contracted_nodes(G1, u, v, self_loops=False))
    return chrom_at(G1, q) - chrom_at(G2, q)
H = lambda G: abs(chrom_at(G, -1))

# fig9: how many different histories a rule package has
fig = plt.figure(figsize=(13, 6)); gs = fig.add_gridspec(2, 4, width_ratios=[1, 1, 1, 1.9], hspace=0.45)
gal = [("nothing clashes", nx.empty_graph(5)), ("a chain of clashes", nx.path_graph(5)),
       ("two separate groups", nx.disjoint_union(nx.path_graph(2), nx.complete_graph(3))),
       ("a cycle", nx.cycle_graph(5)), ("one rule clashes with all", nx.star_graph(4)),
       ("everything clashes", nx.complete_graph(5))]
for k, (ttl, G) in enumerate(gal):
    a = fig.add_subplot(gs[k // 3, k % 3]); pos = nx.circular_layout(G)
    nx.draw_networkx(G, pos, ax=a, with_labels=False, node_size=140, node_color=BLUE, edge_color=ORANGE, width=2)
    hh = H(G); a.set_title(f"{ttl}\n{hh} histor{'y' if hh == 1 else 'ies'} (of 120 orders)", fontsize=9); a.axis("off"); a.set_aspect("equal")
a = fig.add_subplot(gs[:, 3]); rng = random.Random(3); ns = range(2, 9)
for m, c in [(2, ORANGE), (4, GREEN), (8, BLUE)]:
    means = []
    for n in ns:
        vals = []
        for _ in range(60):
            G = nx.Graph(); G.add_nodes_from(range(n)); R = []
            for i in range(n):
                w = rng.randrange(m); R.append((w, {j for j in range(m) if j != w and rng.random() < 0.25}))
            for i, j in itertools.combinations(range(n), 2):
                (wi, ri), (wj, rj) = R[i], R[j]
                if wi == wj or wi in rj or wj in ri: G.add_edge(i, j)
            vals.append(H(G))
        means.append(np.mean(vals))
    a.semilogy(list(ns), means, "o-", color=c, lw=2, label=f"{m} sliders")
a.semilogy(list(ns), [math.factorial(n) for n in ns], "--", color=GREY, lw=1.5, label="every order differs (n!)")
a.set(xlabel="rules in the package", ylabel="different histories (average)",
      title="the more sliders the rules spread over,\nthe fewer histories really differ")
a.legend(fontsize=8)
fig.suptitle("Two rules clash if they write the same slider toward different targets, or one reads what the other writes.\n"
             "The number of genuinely different histories is the number of ways to point every clash one way with no loop: |χ(−1)|",
             fontsize=10, y=1.02)
save(fig, "fig9-histories-from-clashes")

# fig10: two rules that read each other, taking turns: five zones on the loop product p = b c
def zones(L, U):
    a1, a2 = math.sqrt(1 - L), math.sqrt(1 - U)
    return [-(2 - L)*(2 - U)/(L*U), -(a1 + a2)**2/(L*U), -(a1 - a2)**2/(L*U), 1.0]
ZN = ["growing swings", "overshoot each turn", "spiral in", "smooth settle", "runaway"]
ZC = [RED, PURPLE, BLUE, GREEN, ORANGE]
fig, ax = plt.subplots(1, 3, figsize=(14, 4.4), gridspec_kw={"width_ratios": [1.1, 1.3, 1.1]})
L, U = 0.6, 0.3; cuts = zones(L, U); lo, hi = -40, 2
edges = [lo] + cuts + [hi]
for z in range(5):
    ax[0].axvspan(edges[z], edges[z + 1], color=ZC[z], alpha=0.35)
ax[0].set_xscale("symlog", linthresh=1); ax[0].set(xlim=(lo, hi), ylim=(0, 1), yticks=[], xlabel="loop product b·c  (strengths 0.6 and 0.3)",
          title="(a) one number, b·c, sorts every\ntwo-way coupling into five behaviours")
for z, xx in enumerate([-25, -12.6, -3, 0.35, 1.5]):
    ax[0].text(xx, 0.55, ZN[z], ha="center", va="center", fontsize=8, rotation=90, clip_on=True)
picks = [(-14, 0), (-12.6, 1), (-2, 2), (0.4, 3), (1.15, 4)]
for p, z in picks:
    xs, ys = [0.9], [0.1]; x, y = 0.9, 0.1; a_, s_ = 0.2, 0.4; b_ = 1.0; c_ = p
    for t in range(12):
        x = (1 - L)*x + L*(a_ + b_*y); y = (1 - U)*y + U*(s_ + c_*x); xs.append(x)
    ax[1].plot(range(13), xs, "o-", color=ZC[z], ms=3, lw=1.6, label=f"b·c = {p}: {ZN[z]}")
    ax[0].plot([p], [0.15], "v", color=ZC[z], ms=8)
ax[1].set(ylim=(-3, 3), xlabel="turn", ylabel="slider X", title="(b) slider X turn by turn, one example per zone")
ax[1].legend(fontsize=7.5, loc="lower left")
S = np.linspace(0.02, 0.98, 300); P = np.linspace(-30, 2, 400); Z = np.zeros((len(S), len(P)))
for i, st in enumerate(S):
    cs = zones(st, st); Z[i] = [sum(p > q for q in cs) for p in P]
ax[2].imshow(Z, origin="lower", aspect="auto", extent=(P[0], P[-1], S[0], S[-1]), cmap=ListedColormap(ZC), alpha=0.55, vmin=-0.5, vmax=4.5)
ax[2].axvline(0, color="k", lw=0.8, ls=":")
ax[2].set(xlabel="loop product b·c", ylabel="strength of both rules",
          title="(c) equal strengths: spiral starts right at b·c = 0,\nso every balancing loop spirals")
ax[2].annotate("overshoot each turn\n(thin band)", xy=(-4.9, 0.6), xytext=(-21, 0.5), fontsize=7.5, ha="center",
               arrowprops=dict(arrowstyle="->", color="k", lw=0.8))
for z, (xx, yy) in [(0, (-25, 0.85)), (2, (-8, 0.25)), (3, (0.45, 0.6)), (4, (1.5, 0.5))]:
    ax[2].text(xx, yy, ZN[z], fontsize=7.5, ha="center", rotation=90 if z in (3, 4) else 0, clip_on=True)
save(fig, "fig10-five-zones")

# fig11: enact f, enact g, repeal f, repeal g: what is left over
def resid(L, U, p):
    F = np.array([[1 - L, L], [0, 1]]); G = np.array([[1, 0], [U*p, 1 - U]])
    return F @ G @ np.linalg.inv(F) @ np.linalg.inv(G)
fig, ax = plt.subplots(1, 2, figsize=(12, 4.4))
S = np.linspace(0.02, 0.98, 300); P = np.linspace(-2, 2, 400); Z = np.zeros((len(S), len(P)))
for i, st in enumerate(S):
    t = 2 + P*(P - 1)*st**4/(1 - st)**2
    Z[i] = np.where(t > 2, 0, np.where(t < -2, 2, 1))
ax[0].imshow(Z, origin="lower", aspect="auto", extent=(P[0], P[-1], S[0], S[-1]),
             cmap=ListedColormap([ORANGE, GREEN, PURPLE]), alpha=0.5, vmin=-0.5, vmax=2.5)
for xx, yy, tx in [(-1.2, 0.4, "squeeze:\nleftovers grow\nevery cycle"), (0.5, 0.35, "rotation:\nleftovers\nwheel round"),
                   (0.5, 0.93, "flip"), (1.55, 0.4, "squeeze")]:
    ax[0].text(xx, yy, tx, ha="center", fontsize=8.5)
ax[0].set(xlabel="loop product b·c", ylabel="strength of both rules",
          title="(a) pass two coupled rules and repeal them again:\nthe world does not come back, and the leftover has a type")
for p, c, lab in [(-0.5, ORANGE, "balancing loop b·c = −0.5: drifts off"), (0.5, GREEN, "reinforcing loop b·c = 0.5: wheels round")]:
    M = resid(0.5, 0.5, p); v = np.array([0.1, 0.05]); pts = [v]
    for _ in range(40): v = M @ v; pts.append(v)
    pts = np.array(pts); ax[1].plot(pts[:, 0], pts[:, 1], "o-", ms=2.5, lw=1, color=c, label=lab)
ax[1].set(xlim=(-1.5, 1.5), ylim=(-1.5, 1.5), xlabel="leftover in X", ylabel="leftover in Y",
          title="(b) the leftover after each pass-and-repeal cycle\n(strengths 0.5 and 0.5, 40 cycles)")
ax[1].legend(fontsize=8); ax[1].set_aspect("equal")
save(fig, "fig11-pass-and-repeal")
print("done")
