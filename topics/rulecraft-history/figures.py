"""Pictures for rulecraft-history (DRAFT). Every curve is computed from the formulas in checks.py; nothing fitted.
Run: python3 figures.py -> figures/fig1..fig8 (PNG)."""
import itertools, math, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
GREEN, ORANGE, BLUE, GREY, RED, PURPLE = "#2a9d55", "#e07b1f", "#2f6fb5", "#9a9a9a", "#c8382f", "#7a4fb5"
def save(fig, name):
    fig.savefig(f"figures/{name}.png", dpi=150, bbox_inches="tight"); plt.close(fig)
def run(order, L, T, x0):
    y = x0
    for i in order: y = (1 - L[i])*y + L[i]*T[i]
    return y

# fig1: a rule is a strength and a target
fig, ax = plt.subplots(1, 3, figsize=(13, 4))
xs = np.linspace(0, 1, 50)
for lam, t, c in [(0.3, 0.8, GREEN), (0.6, 0.2, ORANGE), (0.9, 0.5, BLUE)]:
    ax[0].plot(xs, (1-lam)*xs + lam*t, color=c, lw=2, label=f"strength {lam}, target {t}")
    ax[0].plot([t], [t], "o", color=c)
ax[0].plot(xs, xs, color=GREY, ls="--", lw=1, label="do nothing")
ax[0].set(xlabel="slider before", ylabel="slider after", title="(a) every rule on a 0-to-1 slider:\npull a share 'strength' of the way to a target")
ax[0].legend(fontsize=8, loc="upper left")
for x0 in np.linspace(0, 1, 11):
    ax[1].annotate("", xy=(0.7*x0 + 0.3*0.8, 1), xytext=(x0, 0), arrowprops=dict(arrowstyle="->", color=GREEN, lw=1))
ax[1].axvline(0.8, color=GREEN, ls=":"); ax[1].text(0.8, 1.05, "target 0.8", ha="center", color=GREEN)
ax[1].set(xlim=(-0.02, 1.02), ylim=(0, 1.12), yticks=[0, 1], yticklabels=["before", "after"], xlabel="slider",
          title="(b) strength 0.3: everything moves 30% of\nthe way to 0.8; the slider shrinks to 70%")
l1 = np.linspace(0, 1, 101)
for l2, c in [(0.2, GREEN), (0.5, ORANGE), (0.8, BLUE)]:
    ax[2].plot(l1, l1 + l2 - l1*l2, color=c, lw=2, label=f"second strength {l2}")
ax[2].set(xlabel="first strength", ylabel="combined strength", title="(c) two rules with the SAME target fuse:\nstrength a + b − ab, our law at k = −1")
ax[2].legend(fontsize=8)
save(fig, "fig1-strength-and-target")

# fig2: order changes the split, not the total
L = [0.30, 0.25, 0.40, 0.20, 0.35]; T = [0.9, 0.1, 0.6, 0.3, 0.75]; names = ["A", "B", "C", "D", "E"]
cols = [GREEN, ORANGE, BLUE, PURPLE, RED]
fig, ax = plt.subplots(1, 2, figsize=(12, 4.3), gridspec_kw={"width_ratios": [1.3, 1]})
orders = [(0, 1, 2, 3, 4), (4, 3, 2, 1, 0), tuple(sorted(range(5), key=lambda i: T[i])), tuple(sorted(range(5), key=lambda i: -T[i]))]
labels = ["A B C D E", "E D C B A", "by target, low→high\n(highest outcome)", "by target, high→low\n(lowest outcome)"]
for k, o in enumerate(orders):
    bottom = 0
    origin = np.prod([1 - l for l in L])
    ax[0].bar(k, origin, color="#ddd", edgecolor="k"); bottom = origin
    for pos, i in enumerate(o):
        w = L[i]*np.prod([1 - L[j] for j in o[pos+1:]])
        ax[0].bar(k, w, bottom=bottom, color=cols[i], edgecolor="k"); 
        if w > 0.05: ax[0].text(k, bottom + w/2, names[i], ha="center", va="center", color="w", weight="bold")
        bottom += w
    ax[0].text(k, 1.02, f"outcome {run(o, L, T, 0.5):.3f}", ha="center", fontsize=8)
ax[0].set_xticks(range(4)); ax[0].set_xticklabels(labels, fontsize=8)
ax[0].set(ylim=(0, 1.1), ylabel="share of the final value", title="(a) five rules, four orders: the grey share (the start, 0.5) and the total\npull are the same every time; order only changes who gets the pull")
outs = sorted(run(o, L, T, 0.5) for o in itertools.permutations(range(5)))
ax[1].plot(range(120), outs, "o", ms=3, color=BLUE)
ax[1].set(xlabel="the 120 orders, sorted by outcome", ylabel="final slider value", title="(b) 120 orders, 120 outcomes; spread %.3f\n(targets span 0.8)" % (outs[-1] - outs[0]))
save(fig, "fig2-split-not-total")

# fig3: the weak order of S4, coloured by outcome
L4 = [0.5, 0.3, 0.7, 0.4]; T4 = [0.2, 0.9, 0.5, 0.7]
asc = tuple(sorted(range(4), key=lambda i: T4[i]))
def inv(o):  # inversions relative to ascending-target order
    r = [asc.index(i) for i in o]; return sum(1 for a in range(4) for b in range(a+1, 4) if r[a] > r[b])
perms = list(itertools.permutations(range(4)))
byrank = {}
for o in perms: byrank.setdefault(inv(o), []).append(o)
pos = {}
for r, ps in byrank.items():
    ps.sort(key=lambda o: run(o, L4, T4, 0.5), reverse=True)
    for j, o in enumerate(ps): pos[o] = ((j - (len(ps)-1)/2)*1.25, -r)
vals = {o: run(o, L4, T4, 0.5) for o in perms}
fig, ax = plt.subplots(figsize=(10, 7))
for o in perms:
    for k in range(3):
        o2 = o[:k] + (o[k+1], o[k]) + o[k+2:]
        if inv(o2) == inv(o) + 1:
            ax.plot(*zip(pos[o], pos[o2]), color="#bbb", lw=1, zorder=1)
vmin, vmax = min(vals.values()), max(vals.values())
sc = ax.scatter([pos[o][0] for o in perms], [pos[o][1] for o in perms], c=[vals[o] for o in perms], cmap="viridis", s=520, zorder=2, edgecolor="k")
for o in perms:
    ax.text(*pos[o], "".join("ABCD"[i] for i in o), ha="center", va="center", fontsize=7, color="w" if vals[o] < (vmin+vmax)/2 else "k")
plt.colorbar(sc, ax=ax, label="final slider value")
ax.set_yticks([-r for r in range(7)]); ax.set_yticklabels([f"{r} swaps\nfrom sorted" for r in range(7)])
ax.set_xticks([])
ax.set_title("Four rules (targets A 0.2, B 0.9, C 0.5, D 0.7; strengths 0.5, 0.3, 0.7, 0.4).\nEvery line is one swap of neighbours; every swap that moves a higher target earlier lowers the outcome.\nTop: sorted by target (highest). Bottom: reversed (lowest). Strengths never change the ends.", fontsize=9)
save(fig, "fig3-orders-are-a-ladder")

# fig4: gradualism and the pendulum
fig, ax = plt.subplots(1, 2, figsize=(12, 4.2))
L2 = [0.6, 0.5]; T2 = [0.1, 0.9]
def inter(order, nn):
    y = 0.5
    for _ in range(nn):
        for i in order:
            li = 1 - (1 - L2[i])**(1/nn); y = (1 - li)*y + li*T2[i]
    return y
ns = np.array([1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 64])
g = [abs(inter((0, 1), n) - inter((1, 0), n)) for n in ns]
ax[0].loglog(ns, g, "o", color=BLUE, label="gap between the two orders")
ax[0].loglog(ns, g[0]/ns, "--", color=GREY, label="1/n")
ax[0].set(xlabel="phase-in steps n (each rule split into n interleaved pieces)", ylabel="order gap",
          title="(a) phase two rules in together: the order gap\nshrinks like 1/n (Lie–Trotter); the limit is adopting both at once")
ax[0].legend(fontsize=8)
y = 0.5; traj = []
for term in range(16):
    i = term % 2; y = (1 - L2[i])*y + L2[i]*T2[i]; traj.append(y)
ax[1].plot(range(1, 17), traj, "o-", color=PURPLE)
l1, l2 = L2; sw = l1*l2*(T2[1]-T2[0])/(l1 + l2 - l1*l2)
ax[1].set(xlabel="term (governments alternate)", ylabel="slider", title="(b) two governments alternate (strength 0.6 to 0.1, 0.5 to 0.9):\nthe slider settles into a swing of size %.3f = order gap / total pull" % abs(sw))
ax[1].axhline(T2[0], color=ORANGE, ls=":"); ax[1].axhline(T2[1], color=GREEN, ls=":")
save(fig, "fig4-phase-in-and-pendulum")

# fig5: coin-flip elections
fig, ax = plt.subplots(1, 4, figsize=(15, 3.5), sharey=False)
rng = np.random.default_rng(1)
for a, lam in zip(ax, [0.25, 1 - (math.sqrt(5)-1)/2, 0.5, 0.7]):
    N = 400000; x = np.full(N, 0.5)
    for _ in range(60):
        eps = rng.integers(0, 2, N); x = (1 - lam)*x + lam*eps
    a.hist(x, bins=400, color=BLUE, density=True)
    ttl = {0.25: "strength 0.25: smooth", 0.5: "strength 0.5: exactly flat", 0.7: "strength 0.7: a Cantor set;\nno position between 0.3 and 0.7"}.get(round(lam, 2), "strength 1 − 1/φ ≈ 0.382:\nspiky (Erdős's singular case)")
    a.set(title=ttl, xlabel="slider", yticks=[])
fig.suptitle("Each term a coin flip picks which side wins; the winner pulls the slider toward 0 or 1 with the same strength. Long-run spread of positions:", fontsize=10, y=1.13)
save(fig, "fig5-forbidden-middle")

# fig6: the memory exponent (decision gravity)
fig, ax = plt.subplots(1, 2, figsize=(12, 4.2))
rng = np.random.default_rng(4); N = 3000
targets = rng.random(N)
for alpha, c in [(0.0, GREY), (0.4, ORANGE), (0.75, GREEN), (1.0, BLUE), (1.6, RED)]:
    y = 0.05; tr = []
    for n in range(N):
        lam = 0.3 if alpha == 0 else min(1.0, 1/(n + 1)**alpha) if alpha <= 1 else 0.6/(n + 1)**alpha
        y = (1 - lam)*y + lam*targets[n]; tr.append(y)
    ax[0].plot(tr, color=c, lw=1, label={0: "α = 0: fixed strength (fading memory)", 0.4: "α = 0.4: never settles",
                                          0.75: "α = 0.75: settles, history decides", 1.0: "α = 1: plain average (Pólya, Arthur)", 1.6: "α = 1.6: freezes near the founders' 0.05"}[alpha])
ax[0].set(xscale="log", xlabel="rules adopted so far (log)", ylabel="slider", title="(a) strength of the n-th rule ∝ 1/n^α: the same random\ntargets, five memory exponents")
ax[0].legend(fontsize=7, loc="upper right")
al = np.linspace(0, 2.5, 200)
ax[1].axvspan(0, 0.5, color="#fde2c8"); ax[1].axvspan(0.5, 1, color="#d5ecd9"); ax[1].axvspan(1, 2.5, color="#f6d0cd")
for xx, txt in [(0.25, "keeps moving\n(sum of strength²\ndiverges)"), (0.75, "settles; start\nforgotten;\nhistory decides"), (1.75, "start never\nforgotten\n(sum of strengths\nconverges)")]:
    ax[1].text(xx, 0.5, txt, ha="center", va="center", fontsize=9)
ax[1].set(xlim=(0, 2.5), ylim=(0, 1), yticks=[], xlabel="memory exponent α", title="(b) three regimes, split at α = 1/2 and α = 1\n(Robbins–Monro conditions)")
save(fig, "fig6-memory-exponent")

# fig7: repeal window and the curved space of rules
fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
lam = np.linspace(0, 1, 200); t = 0.3
ax[0].fill_between(lam, lam*t, 1 - lam*(1 - t), color="#d5ecd9", label="slider values where 'undo' still works")
ax[0].plot(lam, lam*t, color=GREEN); ax[0].plot(lam, 1 - lam*(1-t), color=GREEN)
ax[0].plot([1], [t], "o", color=RED); ax[0].text(0.97, t + 0.05, "full strength:\nno way back", ha="right", color=RED, fontsize=8)
ax[0].set(xlabel="strength of the rule being repealed (target 0.3)", ylabel="slider value at repeal",
          title="(a) repeal by running the rule backward works only\ninside a window of width 1 − strength")
ax[0].legend(fontsize=8, loc="lower left")
# upper half-plane: reforms x -> A x + C as (C, A)
tri = plt.Polygon([[0, 0], [1, 0], [0, 1]], closed=True, color="#e3ecf8")
ax[1].add_patch(tri)
ax[1].text(0.03, 0.12, "rules that keep\na 0-to-1 slider\ninside 0 to 1", fontsize=8, color=BLUE)
for (l, tt, c, nm) in [(0.3, 0.8, GREEN, "f"), (0.6, 0.2, ORANGE, "g")]:
    ax[1].plot([l*tt], [1 - l], "o", color=c); ax[1].text(l*tt + 0.02, 1 - l + 0.02, nm, color=c, weight="bold")
ax[1].plot([0], [1], "ks"); ax[1].text(0.02, 1.02, "do nothing", fontsize=8)
for c0, r in [(0.5, 0.5), (0.2, 0.9), (0.9, 0.6)]:
    th = np.linspace(0, math.pi, 200); ax[1].plot(c0 + r*np.cos(th), r*np.sin(th), color=GREY, lw=0.8)
ax[1].axhline(0, color=RED, lw=2); ax[1].text(1.3, 0.03, "strength 1 ('set it to t'): the edge at infinity", color=RED, fontsize=8, ha="right")
ax[1].set(xlim=(-0.6, 1.4), ylim=(-0.05, 1.3), xlabel="shift C", ylabel="stretch A = 1 − strength",
          title="(b) rules x → Ax + C as points; with the distance that\ncomposing rules keeps, this is the hyperbolic plane (curvature −1);\ngrey arcs are its straight lines")
save(fig, "fig7-repeal-window-and-curved-space")

# fig8: anchors, counted by Delta = k^2 - 4d
fig, ax = plt.subplots(1, 3, figsize=(13, 4.3))
def circ(a):
    th = np.linspace(0, 2*math.pi, 300); a.plot(np.cos(th), np.sin(th), color="#ccc"); a.set_aspect("equal"); a.axis("off")
def ang(x):  # projective line to circle via x = tan(theta/2)
    return 2*np.arctan(x)
for a, (k, d, ttl) in zip(ax, [(1.0, -2.0, "Δ > 0: two anchors\n(odds-type rules; never reach either end)"),
                               (-1.0, 0.0, "Δ > 0 with d = 0: one anchor is ∞\n(nudges; our ∘₋₁ pulls to 1)"),
                               (0.0, 1.0, "Δ < 0: no anchor\n(rotation; the slider is a circle)")]):
    circ(a)
    xs = np.tan(np.linspace(-math.pi/2 + 0.05, math.pi/2 - 0.05, 22))
    for x0 in xs:
        x1 = (x0 + 0.25 + k*x0*0.25)/(1 - d*x0*0.25)
        a0, a1 = ang(x0), ang(x1)
        da = (a1 - a0 + math.pi) % (2*math.pi) - math.pi
        if abs(da) < 1e-3: continue
        thr = np.linspace(a0, a0 + da, 20)
        a.plot(0.9*np.cos(thr), 0.9*np.sin(thr), color=BLUE, lw=1)
        a.annotate("", xy=(0.9*np.cos(thr[-1]), 0.9*np.sin(thr[-1])), xytext=(0.9*np.cos(thr[-3]), 0.9*np.sin(thr[-3])), arrowprops=dict(arrowstyle="->", color=BLUE))
    disc = k*k - 4*d
    if disc >= 0:
        rts = [(-k + s*math.sqrt(disc))/(2*d) for s in (1, -1)] if d != 0 else [-1/k, math.inf]
        for r in rts:
            th = math.pi if r == math.inf else ang(r)
            a.plot(np.cos(th), np.sin(th), "o", color=RED, ms=9)
            a.text(1.15*np.cos(th), 1.15*np.sin(th), "∞" if r == math.inf else f"{r:g}", ha="center", va="center", color=RED)
    a.set_title(f"k = {k:g}, d = {d:g}: " + ttl, fontsize=9, pad=22)
fig.suptitle("The rule x → x ⊕ c of each (k, d) law, drawn on the projective line (a circle). Red: the anchors it never moves, the roots of 1 + kx + dx².", fontsize=10, y=1.02)
save(fig, "fig8-anchors")
print("done")
