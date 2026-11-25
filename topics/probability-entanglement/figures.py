"""Pictures for the DRAFT write-up. Run after checks.py (reads results.json)."""
import json, math
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle, FancyArrowPatch, Arc
from lib import *
R = json.load(open("results.json")); rng = np.random.default_rng(7)
BLUE, ORANGE, GREEN, RED, PURPLE, GRAY = "#2a6fdb", "#e8833a", "#2e9e5b", "#d64545", "#8a5cc2", "#888888"
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
def save(fig, name): fig.tight_layout(); fig.savefig(f"figures/{name}.png", dpi=140); plt.close(fig)

# ---------------- fig1: probability from entanglement
fig, ax = plt.subplots(1, 3, figsize=(14, 4.3))
p = np.array([.4, .25, .2, .1, .05])
ax[0].bar(range(5), p, color=BLUE); ax[0].set_xticks(range(5)); ax[0].set_xticklabels([f"|{i}⟩|{i}⟩" for i in range(5)])
ax[0].set_ylabel("weight"); ax[0].set_title("(a) Purification: any probability list p\nis the 'shadow' of one pure entangled whole")
ax[0].text(2.2, .33, "whole = Σ √pᵢ |i⟩_system |i⟩_outside\nsystem alone sees p\nShannon entropy = entanglement entropy", fontsize=9)
a = ax[1]; a.axis("off"); a.set_xlim(0, 10); a.set_ylim(0, 6)
a.set_title("(b) Zurek's envariance: a swap on the system\nis undone by a swap on the environment")
for y, lab in [(4.2, "|↑⟩_S|e₁⟩_E + |↓⟩_S|e₂⟩_E"), (2.6, "|↓⟩_S|e₁⟩_E + |↑⟩_S|e₂⟩_E"), (1.0, "|↓⟩_S|e₂⟩_E + |↑⟩_S|e₁⟩_E  = start")]:
    a.text(5, y, lab, ha="center", fontsize=10, bbox=dict(fc="#eef3fb", ec=BLUE))
a.annotate("", (5, 3.0), (5, 3.9), arrowprops=dict(arrowstyle="->", color=ORANGE, lw=2)); a.text(5.3, 3.35, "swap on S", color=ORANGE)
a.annotate("", (5, 1.4), (5, 2.3), arrowprops=dict(arrowstyle="->", color=GREEN, lw=2)); a.text(5.3, 1.75, "swap on E (far away)", color=GREEN)
a.text(5, 0.2, "nothing local changed, so ↑ and ↓ must be equally likely", ha="center", fontsize=9, style="italic")
a = ax[2]; Ma, Mb, aa, bb = 8, 5, 3, 2
for i in range(Ma):
    for j in range(Mb):
        hit = i < aa or j < bb
        a.add_patch(Rectangle((i, j), .95, .95, color=(ORANGE if i < aa and j < bb else BLUE) if hit else "#e5e5e5"))
a.set_xlim(-.2, Ma + .2); a.set_ylim(-.2, Mb + .2); a.set_aspect("equal"); a.set_xticks([]); a.set_yticks([])
a.set_xlabel(f"system 1: {Ma} equal branches, event A in {aa}"); a.set_ylabel(f"system 2: {Mb} equal branches, B in {bb}")
a.set_title(f"(c) Union rule by counting branches\n'A or B' = {aa*Mb+Ma*bb-aa*bb}/{Ma*Mb} = 3/8 ∘₋₁ 2/5 (orange counted once)")
save(fig, "fig1_probability_from_entanglement")

# ---------------- fig2: union rule = no correlation; typicality
fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
pts = np.array([(I, b, e, k) for I, b, e, k in R["B3"]], dtype=float)
cols = {0: BLUE, 1: GREEN, 2: PURPLE}; labs = {0: "random pure (entangled)", 1: "classical mixtures (no entanglement)", 2: "noisy entangled"}
for k in (1, 0, 2):
    m = pts[:, 3] == k; ax[0].scatter(pts[m, 0], pts[m, 1], s=5, alpha=.4, color=cols[k], label=labs[k])
I = np.linspace(0, pts[:, 0].max(), 200); ax[0].plot(I, np.sqrt(2 * I), color=RED, lw=2, label="ceiling √(2·I)")
ax[0].set_xlabel("mutual information I(A:B) between the two parts (nats)"); ax[0].set_ylabel("largest miss of the union rule\n|P(a or b) − p ∘₋₁ q| (12 random tests)")
ax[0].set_title("(a) The union rule holds exactly when the parts share nothing;\nany correlation, classical or quantum, shows up as a miss"); ax[0].legend(fontsize=8, markerscale=3)
C = np.array(R["C"])
ax[1].loglog(C[:, 0], C[:, 1], "o-", color=BLUE, label="mutual information of two small parts")
ax[1].loglog(C[:, 0], C[:, 2], "s-", color=ORANGE, label="largest union-rule miss")
ax[1].loglog(C[:, 0], C[:, 3], "^-", color=GREEN, label="distance of the pair from 'fair coins'")
ax[1].set_xlabel("size of the rest of the world (dimension of environment)"); ax[1].legend(fontsize=8)
ax[1].set_title("(b) Typicality: two small parts of a big entangled whole\nlook random AND independent, so they obey the union rule")
save(fig, "fig2_union_rule_means_no_correlation")

# ---------------- fig3: time from entanglement
fig = plt.figure(figsize=(14, 4.6))
a = fig.add_subplot(1, 3, 1, projection="3d")
t = np.linspace(0, 3, 300); ph = 2 * math.pi * t
a.plot(t, np.cos(ph), np.sin(ph), color=BLUE, lw=2)
a.scatter(np.arange(0, 3.01, 1 / 8), np.cos(2 * math.pi * np.arange(0, 3.01, 1 / 8)), np.sin(2 * math.pi * np.arange(0, 3.01, 1 / 8)), color=ORANGE, s=12)
a.set_xlabel("clock reading"); a.set_ylabel("Re"); a.set_zlabel("Im")
a.set_title("(a) Page–Wootters: the whole is frozen, but the\nsystem's phase, read against the clock, is a helix")
D = np.array(R["D"])
ax2 = fig.add_subplot(1, 3, 2); ax2.plot(np.degrees(D[:, 0]), D[:, 1], "o-", color=BLUE)
ax2.set_xlabel("qubit's tilt between two energy levels (degrees)"); ax2.set_ylabel("clock–system entanglement (bits)")
ax2.set_title("(b) No change, no time: entanglement with\nthe clock is zero exactly when nothing evolves")
ax2.annotate("stationary\n(energy level)", (0, 0), (20, .35), arrowprops=dict(arrowstyle="->")); ax2.annotate("stationary", (180, 0), (130, .35), arrowprops=dict(arrowstyle="->"))
a3 = fig.add_subplot(1, 3, 3); n = 7
for k in range(n):
    x, y = math.cos(2 * math.pi * k / n), math.sin(2 * math.pi * k / n)
    a3.add_patch(Circle((x - 1.2, y), .12, color=BLUE)); a3.text(x * 1.3 - 1.2, y * 1.3, str(k), ha="center", va="center", fontsize=8)
    a3.plot([1.2, 1.8], [k * .35 - 1.05] * 2, color=ORANGE, lw=2); a3.text(1.9, k * .35 - 1.05, f"E = 2π·{k}/{n}", va="center", fontsize=8)
a3.annotate("", (-1.2 + .3, .95), (-1.2 + .95, .3), arrowprops=dict(arrowstyle="->", color=GRAY, connectionstyle="arc3,rad=.3"))
a3.text(-1.2, -1.5, "deterministic clock\nk → k+1 mod 7", ha="center"); a3.text(1.6, -1.6, "same thing as a quantum\nsystem with even energy steps", ha="center")
a3.set_xlim(-2.6, 3); a3.set_ylim(-1.9, 1.6); a3.axis("off"); a3.set_title("(c) 't Hooft's cogwheel: a 'clock mod n'\nIS a quantum system")
save(fig, "fig3_time_from_entanglement")

# ---------------- fig4: small loop read coarsely
fig, ax = plt.subplots(1, 4, figsize=(17, 4.3))
E1 = np.array(R["E1"]); ax[0].semilogy(E1[:, 0], np.maximum(E1[:, 1], 1e-4), "o-", color=BLUE)
ax[0].set_xlabel("loop speed (turns per unit of reading-time spread ×2π)"); ax[0].set_ylabel("|P(phase in 1/3 sector) − 1/3|")
ax[0].set_title("(a) A fast loop read at fuzzy times\ngives fair odds (method of arbitrary functions)")
T = np.linspace(0, 60, 6000); gold = (1 + 5 ** .5) / 2
for r, c, lab in [(gold, BLUE, "speeds 1 : golden ratio"), (2, RED, "speeds 1 : 2 (locked)")]:
    ax[1].plot(np.mod(T, 1), np.mod(r * T, 1), ",", color=c, alpha=.6)
    ax[1].plot([], [], color=c, label=lab)
ax[1].add_patch(Rectangle((0, 0), .3, 1, fc=ORANGE, alpha=.15)); ax[1].add_patch(Rectangle((0, 0), 1, .4, fc=GREEN, alpha=.15))
ax[1].set_xlabel("loop 1 phase (A = first 30%)"); ax[1].set_ylabel("loop 2 phase (B = first 40%)"); ax[1].legend(fontsize=8, loc="upper right")
ax[1].set_title(f"(b) Two loops: unlocked fill the torus → union rule\n(miss {R['E3']['golden']:+.4f}); locked don't (miss {R['E3']['two']:+.3f})")
orb = R["E6"]["orbit"]; ax[2].plot(range(1, len(orb) + 1), orb, ".", color=PURPLE, ms=4)
ax[2].set_xlabel("k (number of ∘₁ steps)"); ax[2].set_ylabel(f"1 ∘₁ 1 ∘₁ … ∘₁ 1 = 2ᵏ − 1  mod {R['E6']['p']}")
ins = ax[2].inset_axes([.62, .62, .36, .33]); ins.bar(range(10), R["E6"]["lastdigit"], color=PURPLE); ins.set_ylim(0, .15); ins.set_title("last digit", fontsize=8); ins.tick_params(labelsize=7)
ax[2].set_title("(c) Our finite ∘₁ clock: one loop through\nevery value but the wall; last digit is fair")
hd = np.array(R["E7"]["clock"]); hr = np.array(R["E7"]["rotation"]); eps = R["E7"]["eps"]
bins = np.linspace(0, 2000, 41)
ax[3].hist(hd, bins, density=True, alpha=.6, color=BLUE, label="chaotic ∘₁ clock")
ax[3].hist(hr, bins, density=True, alpha=.5, color=RED, label="plain rotation")
xx = np.linspace(0, 2000, 200); ax[3].plot(xx, eps * np.exp(-eps * xx), color="k", lw=2, label="exponential (Poisson clock)")
ax[3].set_xlabel("steps until the loop first lands in a 0.2% window"); ax[3].legend(fontsize=8)
ax[3].set_title("(d) A chaotic loop turns into a decay clock:\nP(ticked by t) = 1 − e^(−t/τ)")
save(fig, "fig4_small_loop_probability")

# ---------------- fig5: Born from a loop and the Bell wall
fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))
a = ax[0]; phs = np.linspace(-math.pi / 2, math.pi / 2, 200)
a.fill_between(np.degrees(phs), 0, np.cos(phs), color=BLUE, alpha=.3, label="loop dwell time ∝ cos(phase)")
B = np.array(R["F1"]); a2 = a.twinx()
a2.plot(np.degrees(B[:, 0]), B[:, 1], "o", color=ORANGE, label="loop model P(+)"); th = np.linspace(0, math.pi, 100)
a2.plot(np.degrees(th), np.cos(th / 2) ** 2, color="k", lw=1, label="Born cos²(θ/2)")
a.set_xlabel("loop phase  /  measurement angle θ (degrees)"); a.set_ylabel("dwell density"); a2.set_ylabel("P(+)")
a.legend(loc="lower left", fontsize=8); a2.legend(loc="upper right", fontsize=8); a.set_title("(a) One qubit: a lopsided loop gives the Born rule\n(Kochen–Specker 1967 model)")
D = np.linspace(-math.pi, math.pi, 200); Esh = [-(1 - 2 * abs(d) / math.pi) for d in D]
ax[1].plot(np.degrees(D), Esh, color=RED, lw=2, label="shared loop, settings independent")
ax[1].plot(np.degrees(D), -np.cos(D), color="k", lw=2, label="quantum (singlet)")
G2 = np.array(R["G2"]); ax[1].plot(np.degrees(G2[:, 0]), G2[:, 1], "o", color=GREEN, label="whole-path loop model")
ax[1].set_xlabel("angle between the two detectors (degrees)"); ax[1].set_ylabel("correlation E"); ax[1].legend(fontsize=8)
ax[1].set_title("(b) Straight lines vs the cosine:\nthe gap is what Bell's theorem measures")
names = ["local loop\n(circle)", "local helix\n(phase+winding)", "whole-path\nloop weight", "path sum\n|Σ phases|²", "whole-path,\nunconstrained"]
vals = [2, 2, 2 * math.sqrt(2), 2 * math.sqrt(2), 4]; cs = [RED, RED, GREEN, GREEN, GRAY]
ax[2].bar(names, vals, color=cs); ax[2].axhline(2, ls="--", color=RED); ax[2].axhline(2 * math.sqrt(2), ls="--", color=GREEN)
ax[2].text(4.45, 2.05, "Bell 2", color=RED, fontsize=8, ha="right"); ax[2].text(4.45, 2.88, "Tsirelson 2√2", color=GREEN, fontsize=8, ha="right")
ax[2].set_ylabel("CHSH score"); ax[2].tick_params(axis="x", labelsize=8)
ax[2].set_title("(c) Only whole-path loop weights reach quantum;\nonly |Σ phases|² stops at exactly 2√2")
save(fig, "fig5_born_and_bell")

# ---------------- fig6: AdS3 cross-ratio as a modular-time Poisson clock
H = R["H3"]; A_, B1, B2 = H["A"], H["B1"], H["B2"]
fig = plt.figure(figsize=(17, 4.6)); a = fig.add_subplot(1, 4, 1)
a.add_patch(Circle((0, 0), 1, fill=False, color="k"))
def geo(u, v, **kw):
    pu, pv = np.array([math.cos(u), math.sin(u)]), np.array([math.cos(v), math.sin(v)])
    s = (v - u) % (2 * math.pi)
    if abs(s - math.pi) < 1e-3: a.plot(*zip(pu, pv), **kw); return
    mid = (u + s / 2); d = 1 / math.cos(s / 2); c = d * np.array([math.cos(mid), math.sin(mid)]); r = abs(math.tan(s / 2))
    angs = np.linspace(0, 2 * math.pi, 400); P = c[:, None] + r * np.array([np.cos(angs), np.sin(angs)]); keep = (P ** 2).sum(0) <= 1.0001
    a.plot(P[0, keep], P[1, keep], **kw)
def inarc(x, lo, hi): return lo <= x <= hi
n = rng.poisson(40)
for _ in range(n):
    u = rng.uniform(0, 2 * math.pi); s = 2 * math.atan(1 / (1 / math.tan(.15) - rng.uniform() / math.tan(.15) * .999)); v = (u + s) % (2 * math.pi)
    link = (inarc(u, *A_) and inarc(v, B1[0], B2[1])) or (inarc(v, *A_) and inarc(u, B1[0], B2[1]))
    geo(u, v, color=ORANGE if link else "#cccccc", lw=2 if link else .7)
for (lo, hi), c, lab in [(A_, BLUE, "A"), (B1, GREEN, "B₁"), (B2, PURPLE, "B₂")]:
    th = np.linspace(lo, hi, 50); a.plot(1.03 * np.cos(th), 1.03 * np.sin(th), color=c, lw=5); m = (lo + hi) / 2; a.text(1.2 * math.cos(m), 1.2 * math.sin(m), lab, color=c, fontsize=12, ha="center")
a.set_aspect("equal"); a.axis("off"); a.set_title("(a) A random spray of geodesics\n(orange: ones linking A to B₁B₂)")
a = fig.add_subplot(1, 4, 2)
z1 = abs(zeta(B1[1], *A_) - zeta(B1[0], *A_)); z2 = abs(zeta(B2[1], *A_) - zeta(B2[0], *A_))
a.add_patch(Rectangle((0, 0), z1, 1, color=GREEN, alpha=.3)); a.add_patch(Rectangle((z1, 0), z2, 1, color=PURPLE, alpha=.3))
ticks = np.array([0.21, 0.62, 1.05])  # one sample universe: rate-1 ticks
for tk in ticks: a.axvline(tk, color=ORANGE, lw=2)
a.text(z1 / 2, 1.05, f"B₁: width t₁={z1:.2f}", ha="center", color=GREEN); a.text(z1 + z2 / 2, 1.05, f"B₂: t₂={z2:.2f}", ha="center", color=PURPLE)
a.set_xlim(-.1, z1 + z2 + .5); a.set_ylim(0, 1.25); a.set_yticks([]); a.set_xlabel("A's own (modular / Rindler) time ζ")
a.set_title("(b) Seen from A, each region is a time window;\nlinks are ticks of a rate-1 Poisson clock (orange)")
a = fig.add_subplot(1, 4, 3); ks = ["B1", "B2", "B12"]; lab = ["A–B₁", "A–B₂", "A–B₁B₂"]
a.bar(np.arange(3) - .18, [H["x"][k] for k in ks], .36, color=BLUE, label="cross-ratio x (formula)")
a.bar(np.arange(3) + .18, [H["sim"][k] for k in ks], .36, color=ORANGE, label="P(at least one link), 3000 sprays")
a.set_xticks(range(3)); a.set_xticklabels(lab); a.legend(fontsize=8)
a.set_title(f"(c) x = P(≥1 link) = 1 − e^(−t);\nx(A,B₁B₂) = {H['x']['B1']:.3f} ∘₋₁ {H['x']['B2']:.3f} = {H['x']['B12']:.3f}")
ff = np.array(R["H5"]); a = fig.add_subplot(1, 4, 4)
a.plot([0, .32], [0, .32], color="k", lw=1); a.plot(ff[:, 4], ff[:, 3], "o", color=RED, ms=7)
a.set_xlabel("t / 3  (clock-window width, c = 1)"); a.set_ylabel("conditional mutual information I(A:C|B) (nats)")
a.set_title("(d) A real quantum state (free fermions):\nI(A:C|B) = (c/3)·t, so x = 1 − e^(−3I/c)")
save(fig, "fig6_cross_ratio_is_a_clock")

# ---------------- fig7: helix
fig = plt.figure(figsize=(15, 4.6)); a = fig.add_subplot(1, 3, 1, projection="3d")
th = np.linspace(0, 6 * math.pi, 600); a.plot(np.cos(th), np.sin(th), th / (2 * math.pi), color=BLUE, lw=2, label="helix (θ kept on a line)")
a.plot(np.cos(th), np.sin(th), np.zeros_like(th), color=ORANGE, lw=3, label="circle (θ mod 2π)")
for z0 in [0, 1, 2, 3]: a.scatter([1], [0], [z0], color=RED, s=20)
a.set_zlabel("turns (winding)"); a.legend(fontsize=8, loc="upper left"); a.set_title("(a) The helix is the loop with its turns kept:\nevery red dot projects to the same clock reading")
a = fig.add_subplot(1, 3, 2); s1 = np.linspace(0, 3.0, 300); a2c = 1.9
tot = s1 + a2c; principal = np.angle(np.exp(1j * tot))
a.plot(s1, tot, color=BLUE, lw=2, label="helix: θ₁ + θ₂ (always additive)"); a.plot(s1, principal, color=ORANGE, lw=2, label="circle: log of E(iθ₁) ∘₁ E(iθ₂), principal")
a.axhline(math.pi, ls=":", color=GRAY); a.text(0.05, math.pi + .1, "π: the circle drops a full turn here", fontsize=8)
a.set_xlabel("θ₁ (with θ₂ = 1.9)"); a.set_ylabel("angle"); a.legend(fontsize=8)
a.set_title("(b) Our lift ℓ₁ = log(1+x) unrolls ∘₁ to +;\non the circle it loses count, the helix keeps it")
a = fig.add_subplot(1, 3, 3)
for n_, c in [(1, BLUE), (2, GREEN)]:
    s = np.linspace(0, 2 * math.pi * n_, 2000); r = 1 + .3 * n_ * 0.6 + .25 * np.sin(5 * s + n_)
    a.plot(r * np.cos(s), r * np.sin(s), color=c, lw=1.2, label=f"winds {n_}×: phase {n_}×Φ, any shape")
a.add_patch(Circle((0, 0), .25, color=RED)); a.text(0, 0, "Φ", color="w", ha="center", va="center", fontsize=12)
a.set_aspect("equal"); a.axis("off"); a.legend(fontsize=8, loc="lower center", bbox_to_anchor=(.5, -.12))
a.set_title("(c) Aharonov–Bohm: the hidden flux is felt only\nthrough the winding count, which noise can't change")
save(fig, "fig7_helix")

# ---------------- fig8: the merged chain
fig, a = plt.subplots(figsize=(15, 5.2)); a.axis("off"); a.set_xlim(0, 15); a.set_ylim(0, 5.2)
boxes = [(0.2, 3.2, "fast loop / helix\n(clock mod n,\ncogwheel)", BLUE, "known"),
         (3.2, 3.2, "read at fuzzy times\n→ fair odds", BLUE, "known (Poincaré, Hopf)"),
         (6.2, 3.2, "chaotic loop + small\nwindow → decay clock\nx = 1 − e^(−t)", BLUE, "known (Poisson law)"),
         (9.2, 3.2, "adjacent windows glue\nby ∘₋₁ (union rule)\nodds X = E(t), ∘₁", ORANGE, "our reading"),
         (12.2, 3.2, "AdS₃ / CFT: t = width in\nA's modular time,\nI(A:C|B) = (c/3) t", ORANGE, "checked; reading ours"),
         (3.2, 0.5, "entanglement → each part's\nprobabilities (purification,\nZurek, typicality)", GREEN, "known"),
         (6.2, 0.5, "entanglement with a clock\n→ time (Page–Wootters);\nstate → modular time", GREEN, "known"),
         (9.2, 0.5, "local hidden loop for\nentangled pairs: CHSH ≤ 2\n(quantum reaches 2√2)", RED, "Bell: the wall"),
         (12.2, 0.5, "whole-path loop weight:\nreaches 2√2, also 4;\nonly |Σ phases|² stops at 2√2", PURPLE, "toy model, ours")]
for x, y, txt, c, tag in boxes:
    a.add_patch(Rectangle((x, y), 2.6, 1.55, fc=c, alpha=.15, ec=c, lw=2)); a.text(x + 1.3, y + .82, txt, ha="center", va="center", fontsize=9)
    a.text(x + 1.3, y + 1.68, tag, ha="center", fontsize=8, color=c, style="italic")
for x0, x1, y in [(2.8, 3.2, 4.0), (5.8, 6.2, 4.0), (8.8, 9.2, 4.0), (11.8, 12.2, 4.0), (5.8, 6.2, 1.3), (8.8, 9.2, 1.3), (11.8, 12.2, 1.3)]:
    a.annotate("", (x1, y), (x0, y), arrowprops=dict(arrowstyle="->", lw=2))
a.annotate("", (13.5, 3.2), (7.5, 2.05), arrowprops=dict(arrowstyle="->", lw=1.5, ls="--", color=GREEN, connectionstyle="arc3,rad=-.2"))
a.text(10.6, 2.45, "same object: a part's reduced state\ngives its probabilities AND its time", fontsize=8, color=GREEN, ha="center")
a.set_title("Where 'probability', 'time', 'the small loop' and our union rule genuinely meet (top row), and where a hidden loop hits Bell's wall (bottom right)")
save(fig, "fig8_merged_picture")
print("figures done")
