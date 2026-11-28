"""DRAFT physics bridges: pictures pb1-pb8, from results.json (run checks.py first)."""
import json, os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle, Rectangle
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "figures")
R = json.load(open(os.path.join(HERE, "results.json")))
GREEN, BLUE, GREY, RED, GOLD, PURPLE, ORANGE = "#2e9e4f", "#1f6fc5", "#8a8a8a", "#c8372d", "#d9a400", "#7b4fa8", "#e08a1e"
plt.rcParams.update({"font.size": 12})
def box(ax, x, y, text, fc, w=2.6, h=0.9, fs=11, tc="white"):
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle="round,pad=0.08", fc=fc, ec="none"))
    ax.text(x, y, text, ha="center", va="center", fontsize=fs, color=tc, wrap=True)
def save(fig, name): fig.savefig(os.path.join(OUT, name), dpi=130, bbox_inches="tight"); plt.close(fig)

# ---------------------------------------------------------------- pb1: the bridge map
fig, ax = plt.subplots(figsize=(18, 13)); ax.set_xlim(-10, 10); ax.set_ylim(-7.6, 7.6); ax.axis("off")
ax.set_title("pb1 (DRAFT). Bridges from Affine Composition Theory to physics theories", fontsize=17)
box(ax, 0, 0, "Affine Composition Theory\nlaw  a ∘ₖ b = a + b + k·ab", PURPLE, w=4.2, h=1.3, fs=14)
parts = {  # the law's pieces that do the bridging (inner ring)
    "S10": (-3.3, 2.6, "law on a 10-sphere\n→ symmetry SO(10)"),
    "f": (3.3, 2.6, "twist family f\n(f = 1 is the law)"),
    "S3": (3.6, -2.4, "law's 3-sphere\n= SU(2); swap rule q"),
    "time": (0, -3.4, "one boosting direction;\nfused shift a ↦ 1 + ka"),
    "mirror": (-3.6, -2.4, "law is its own mirror;\ntwisted shapes"),
}
for kx, (x, y, t) in parts.items():
    box(ax, x, y, t, "#444444", w=3.3, h=1.1, fs=11)
    ax.plot([0, x * 0.72], [0, y * 0.72], color="#444444", lw=2, zorder=0)
theories = [  # name, pos, inner part, kind, lead
    ("SO(10) grand\nunification", (-7.5, 5.6), "S10", "same", True),
    ("Superstrings\n(heterotic E8×E8)", (-2.6, 6.4), "S10", "same_end", False),
    ("M-theory\n(11D, Hořava–Witten)", (3.0, 6.4), "f", "family", True),
    ("Einstein–Cartan\ntorsion", (7.6, 5.0), "f", "same", False),
    ("Regge, loop quantum\ngravity, spin foams", (7.6, -0.2), "S3", "same3d", True),
    ("Kaluza–Klein", (7.2, -5.8), "S3", "mixed", False),
    ("Page–Wootters, thermal\ntime, neo-Lorentzian", (0, -6.5), "time", "resemb", False),
    ("Pin± and topology\n(Z16)", (-7.4, -5.0), "mirror", "arith", False),
    ("10D spinors\n(16 of SO(9,1))", (-6.9, 0.9), "S10", "arith", False),
]
col = dict(same=GREEN, same_end=GREEN, family=GREEN, same3d=GREEN, mixed=BLUE, arith=BLUE, resemb=GREY)
lab = dict(same="same thing", same_end="same group and 16;\ndifferent road", family="same family,\nother member",
           same3d="same group; exact in 3D,\nopen in 4D", mixed="same mechanism (known);\n14D, not 10 or 11",
           arith="same arithmetic", resemb="resemblance")
for name, (x, y), part, kind, lead in theories:
    px, py, _ = parts[part]
    ax.plot([px, x], [py, y], color=col[kind], lw=4 if kind in ("same", "same_end", "family", "same3d") else 2.5,
            ls="-" if kind != "resemb" else "--", zorder=0)
    box(ax, x, y, name, col[kind], w=3.4, h=1.15, fs=11.5)
    ax.text((px + x) / 2, (py + y) / 2 + 0.25, lab[kind], fontsize=9, color=col[kind], ha="center",
            bbox=dict(fc="white", ec="none", alpha=0.85))
    if lead: ax.plot(x + 1.75, y + 0.55, marker="*", ms=26, color=GOLD, mec="k", zorder=5)
# Pin to SO(10): the 16 ties
ax.annotate("", xy=(-9.4, 5.0), xytext=(-9.4, -4.4), arrowprops=dict(arrowstyle="<->", color=BLUE, lw=2, ls=":"))
ax.text(-9.3, -2.6, "16 per generation\n= Z16-safe\n(Spin(10) centre Z4)", color=BLUE, fontsize=10)
handles = [plt.Line2D([], [], color=GREEN, lw=4, label="same thing (or same family / same group)"),
           plt.Line2D([], [], color=BLUE, lw=2.5, label="same arithmetic or known mechanism, different object"),
           plt.Line2D([], [], color=GREY, lw=2.5, ls="--", label="resemblance only"),
           plt.Line2D([], [], color=GOLD, marker="*", ms=18, lw=0, mec="k", label="opens a lead on its own merits")]
ax.legend(handles=handles, loc="upper center", fontsize=12, frameon=True, ncol=2, bbox_to_anchor=(0.5, 0.0))
save(fig, "pb1_bridge_map.png")

# ---------------------------------------------------------------- pb2: two roads to SO(10)
fig, axs = plt.subplots(1, 3, figsize=(21, 7.6)); fig.subplots_adjust(top=0.82); fig.suptitle("pb2 (DRAFT). Two roads to the same SO(10) and the same 16 (checks B2, B3)", fontsize=15)
ax = axs[0]; ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.set_title("A. the roads")
box(ax, 2.5, 9, "law on a 10-sphere\n(round 36)", PURPLE, w=4, h=1.2)
box(ax, 7.5, 9, "heterotic string:\nE8 × E8 in 10D", ORANGE, w=4, h=1.2)
box(ax, 2.5, 6.2, "fix the point 1:\nSO(11) → SO(10)", "#444444", w=4, h=1.2)
box(ax, 7.5, 6.2, "SU(4) bundle on 6 hidden\ndims: E8 → SO(10) × SU(4)", "#444444", w=4, h=1.4, fs=10)
box(ax, 5, 3.4, "SO(10), matter in 16s\nsin²θ_W = 3/8 where it holds", GREEN, w=5.6, h=1.3)
box(ax, 2.5, 0.9, "middle stage: missing\n(round 36)", RED, w=4, h=1.2)
box(ax, 7.5, 0.9, "Wilson lines → SM × U(1)B−L\n(Braun–He–Ovrut–Pantev 2005)", BLUE, w=4.4, h=1.3, fs=10)
for (x0, y0, x1, y1) in [(2.5, 8.35, 2.5, 6.85), (7.5, 8.35, 7.5, 6.95), (2.5, 5.55, 4.2, 4.05), (7.5, 5.45, 5.8, 4.05), (4.2, 2.75, 2.5, 1.55), (5.8, 2.75, 7.5, 1.55)]:
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0), arrowprops=dict(arrowstyle="->", lw=2))
ax = axs[1]; d = R["B3"]
parts_ = [("(45,1) SO(10)", d["adj_SO10"], GREEN), ("(1,15) SU(4)", d["adj_SU4"], GREY), ("(10,6)", d["rep_10_6"], BLUE),
          ("(16,4)", d["rep_16_4"], PURPLE), ("(16bar,4bar)", d["rep_16b_4b"], ORANGE)]
left = 0
for n, v, c in parts_:
    ax.barh(0, v, left=left, color=c); ax.text(left + v / 2, 0, f"{n}: {v}", ha="center", va="center", color="white", fontsize=11, rotation=90); left += v
ax.set_xlim(0, 248); ax.set_yticks([]); ax.set_xlabel("dimension (of 248)")
ax.set_title("B. E8 split by its 240 roots: 248 = 45 + 15 + 60 + 64 + 64")
ax = axs[2]
st = R["B2"]["states"]; from fractions import Fraction as Fr
pts = {}
for s in st: pts.setdefault((float(Fr(s["Y"])), float(Fr(s["Q"]))), []).append(s["name"])
cmap = {"Q": BLUE, "u^c": ORANGE, "d^c": ORANGE, "L": GREEN, "e^c": GREEN, "nu^c": GREEN}
for (y, q), names in pts.items():
    ax.scatter(y, q, s=260 * len(names), color=cmap[names[0]], alpha=0.8, edgecolor="k")
    ax.text(y + 0.06, q + 0.08, f"{names[0]} ×{len(names)}", fontsize=10)
ax.axhline(0, color=GREY, lw=0.6); ax.axvline(0, color=GREY, lw=0.6)
ax.set_xlabel("hypercharge Y"); ax.set_ylabel("electric charge Q")
ax.set_title("C. one 16 = one generation + ν_R\n(blue triplets, orange antitriplets, green singlets)\nanomalies all 0; X = 5(B−L) − 4Y = 1 mod 4", fontsize=11); ax.set_xlim(-0.85, 1.25); ax.set_ylim(-1.2, 1.25)
save(fig, "pb2_two_roads_to_SO10.png")

# ---------------------------------------------------------------- pb3: the three 16s
fig, axs = plt.subplots(1, 3, figsize=(20, 6.5)); fig.suptitle("pb3 (DRAFT). Three 16s: same number, three different objects (checks B1, B2, B7)", fontsize=15); fig.subplots_adjust(top=0.8)
ax = axs[0]; ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(0, 10)
ax.set_title("A. spinors: 32 = 16 + 16bar")
for i, (lab_, sig, real, c) in enumerate([("SO(10): the law's hidden symmetry", "(10,0)", "complex: conjugate is the other hand", PURPLE),
                                           ("SO(9,1): superstring spacetime", "(9,1)", "can be real: Majorana–Weyl, 16 real numbers", ORANGE)]):
    y = 7.5 - 4 * i
    box(ax, 5, y + 1, lab_, c, w=8.5, h=1.0)
    for j in range(16): ax.add_patch(Rectangle((0.8 + j * 0.27, y - 0.5), 0.24, 0.6, color=c, alpha=0.9))
    for j in range(16): ax.add_patch(Rectangle((5.3 + j * 0.27, y - 0.5), 0.24, 0.6, color=c, alpha=0.35))
    ax.text(2.9, y - 1.0, "16 (left)", ha="center"); ax.text(7.4, y - 1.0, "16 (right)", ha="center")
    ax.text(5, y - 1.7, real, ha="center", fontsize=10.5, color=c)
ax = axs[1]
ax.bar(["with ν_R\n3 × 16", "without ν_R\n3 × 15"], [48, 45], color=[GREEN, RED])
for x, v, m in [(0, 48, 0), (1, 45, 13)]: ax.text(x, v + 1, f"{v} = {m} mod 16", ha="center", fontsize=12)
ax.set_ylim(0, 60); ax.set_ylabel("left-handed Weyl fermions with X = 1 mod 4")
ax.set_title("B. the Z4 = centre of Spin(10) anomaly is counted mod 16\n(García-Etxebarria–Montero 2019)")
ax = axs[2]; ax.axis("off"); ax.set_xlim(-1.6, 1.6); ax.set_ylim(-1.6, 1.6); ax.set_aspect("equal")
ax.set_title("C. Pin+ and Pin−: a mirror's spinor map squares to +1 or −1\nPin+ bordism in 4D is Z16 (Kirby–Taylor 1990)")
t = np.linspace(0, 2 * np.pi, 200); ax.plot(np.cos(t), np.sin(t), color=GREY)
for ang, lab_, c in [(0, "+1", GREEN), (np.pi, "−1", RED), (np.pi / 2, "+i", BLUE), (-np.pi / 2, "−i", BLUE)]:
    ax.plot(np.cos(ang), np.sin(ang), "o", ms=14, color=c); ax.text(1.25 * np.cos(ang) - 0.1, 1.25 * np.sin(ang) - 0.05, lab_, fontsize=14, color=c)
ax.annotate("", xy=(0, 1), xytext=(1, 0), arrowprops=dict(arrowstyle="->", color=BLUE, connectionstyle="arc3,rad=0.3", lw=2))
ax.annotate("", xy=(-1, 0), xytext=(0, 1), arrowprops=dict(arrowstyle="->", color=BLUE, connectionstyle="arc3,rad=0.3", lw=2))
ax.text(-1.5, -1.55, "Pin−: once round a Klein twist ±i, twice −1 (forced)\nPin+: once round ±1, twice +1", fontsize=10.5)
save(fig, "pb3_three_16s.png")

# ---------------------------------------------------------------- pb4: dimensions
fig, ax = plt.subplots(figsize=(14, 6.5)); ax.set_title("pb4 (DRAFT). Dimension bookkeeping (check B4): visible 4 + hidden", fontsize=15)
rows = [("our spacetime", 4, 0, GREY), ("law, 3-sphere hidden (weak SU(2))", 4, 3, PURPLE), ("superstring", 4, 6, ORANGE),
        ("M-theory / 11D supergravity", 4, 7, ORANGE), ("Witten's minimal KK for SU(3)×SU(2)×U(1)", 4, 7, BLUE),
        ("law, 10-sphere hidden (SO(10))", 4, 10, PURPLE), ("bosonic string", 4, 22, GREY)]
for i, (n, v, h, c) in enumerate(rows):
    ax.barh(i, v, color="#bbbbbb"); ax.barh(i, h, left=v, color=c)
    ax.text(v + h + 0.3, i, f"{v + h}", va="center", fontsize=12)
ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[0] for r in rows]); ax.invert_yaxis(); ax.set_xlabel("dimensions (grey = visible 4, colour = hidden)")
ax.axvline(11, color=ORANGE, ls=":"); ax.axvline(10, color=ORANGE, ls=":")
ax.text(15.2, 5.35, "the law's 10 is a sphere of hidden directions,\nnot spacetime: 4 + 10 = 14 matches no string theory", fontsize=10.5, color=PURPLE)
save(fig, "pb4_dimensions.png")

# ---------------------------------------------------------------- pb5: Cartan-Schouten family
fig, axs = plt.subplots(1, 2, figsize=(18, 6.5)); fig.suptitle("pb5 (DRAFT). The twist family is the Cartan–Schouten family (checks B5, B6)", fontsize=15); fig.subplots_adjust(top=0.8)
ax = axs[0]; fv = np.linspace(0, 2, 200)
ax.plot(fv, 1 - fv, color=GREY, ls="--", lw=2.5, label="torsion ∝ 1 − f (straight line)")
ax.plot(fv, (2 * fv - fv ** 2) / 4 * 4, color=BLUE, lw=2.5, label="curvature ∝ 2f − f² (arch)")
for x, y, t, c in [(1, 0, "law f = 1: torsion 0,\nLevi-Civita\n(Einstein–Cartan in vacuum)", GREEN), (0, 1, "f = 0: flat, torsion\n(teleparallel / Weitzenböck;\non S⁷: Englert's torsion)", ORANGE),
                   (2, -1, "f = 2: flat, torsion\n(the mirror image)", ORANGE)]:
    ax.plot(x, y, "o", ms=14, color=c); ax.text(x + (0.08 if x < 2 else -0.62), y + (-0.33 if x == 0 else (0.1 if y >= 0 else 0.12)), t, fontsize=10.5, color=c)
ax.set_xlabel("f  (member x^(1−f/2) · y · x^(f/2))"); ax.legend(loc="lower left"); ax.axhline(0, color="k", lw=0.5)
ax.set_title("A. on the 3-sphere = SU(2): exact identities (B5)")
ax = axs[1]; g = R["B6"]["order_gap_rel_err"]; g5 = R["B5"]["order_gap_rel_err"]
ax.semilogy([float(k) for k in g5], list(g5.values()), "o-", color=PURPLE, ms=10, label="3-sphere (quaternions)")
ax.semilogy([float(k) for k in g], list(g.values()), "s-", color=ORANGE, ms=10, label="7-sphere (octonions)")
ax.axhline(5e-3, color=RED, ls=":", label="pass line")
ax.set_xlabel("f"); ax.set_ylabel("relative error of  order gap = (1 − f)[a, b]")
ax.set_title("B. measured order gap vs (1 − f)[a, b] at step 10⁻³;\non S⁷, [e_i, e_j] = 2c_ijk e_k: the octonion constants (7 Fano lines)")
ax.legend(); save(fig, "pb5_cartan_schouten.png")

# ---------------------------------------------------------------- pb6: SU(2) to quantum gravity
fig, ax = plt.subplots(figsize=(17, 6.5)); ax.axis("off"); ax.set_xlim(0, 17); ax.set_ylim(0, 6.5)
ax.set_title("pb6 (DRAFT). The law's 3-sphere is SU(2), the group of loop quantum gravity and spin foams", fontsize=15)
steps = [(1.8, "law's 3-sphere\n= unit quaternions\n= SU(2)", PURPLE), (5.2, "spin labels j on edges\n(spin networks;\narea ∝ √(j(j+1)))", GREEN),
         (8.6, "3D: Ponzano–Regge sum\n→ Regge action\n(exact; repo action/)", GREEN), (12.0, "3D: q face weight\n→ Turaev–Viro, Λ > 0\n(repo action/, E3)", GREEN),
         (15.3, "4D: q-deformed EPRL\n(Han 2011; Fairbairn–\nMeusburger 2012): Λ > 0", GOLD)]
for x, t, c in steps: box(ax, x, 3.6, t, c, w=2.7, h=2.0, fs=10.5, tc="white" if c != GOLD else "black")
for i in range(4): ax.annotate("", xy=(steps[i + 1][0] - 1.45, 3.6), xytext=(steps[i][0] + 1.45, 3.6), arrowprops=dict(arrowstyle="->", lw=2.5))
ax.text(15.3, 1.7, "open: does the swap rule's\nface weight give this 4D\ndeformation? (the lead)", ha="center", fontsize=11, color="#8a6d00")
ax.text(8.6, 1.4, "what the law does not set: the Immirzi parameter γ;\nwhy 4D Einstein gravity is not renormalisable", ha="center", fontsize=10.5, color=RED)
save(fig, "pb6_su2_quantum_gravity.png")

# ---------------------------------------------------------------- pb7: clocks
fig, axs = plt.subplots(1, 2, figsize=(17, 7)); fig.suptitle("pb7 (DRAFT). Two clocks: Page–Wootters and the law's fused shift (check B8)", fontsize=15); fig.subplots_adjust(top=0.8)
ax = axs[0]; N = R["B8"]["N"]
t = np.arange(N); ph = (2 * np.pi * 3 * t / N) % (2 * np.pi)
ax.plot(t, ph, "o", color=BLUE, ms=6, label="relative phase of the qubit, given clock reading t")
ax.plot(t, 2 * np.pi * 3 * t / N - 2 * np.pi * np.floor(3 * t / N), color=GREY, lw=0.8)
ax.set_xlabel("clock reading t"); ax.set_ylabel("phase (radians, mod 2π)"); ax.legend(loc="upper left", fontsize=10)
ax.set_title("A. Page–Wootters: the joint state never changes (frozen),\nyet the qubit read at each tick t turns at a steady rate (Schrödinger)")
ax = axs[1]
k = np.exp(2j * np.pi / N); a = 0.3 + 0.1j; astar = 1 / (1 - k); tr = [a]
for _ in range(N): tr.append(1 + k * tr[-1])
tr = np.array(tr); ax.plot(tr.real, tr.imag, "o-", color=PURPLE, ms=5, lw=1)
ax.plot(astar.real, astar.imag, "*", ms=18, color=GOLD, mec="k"); ax.text(astar.real + 0.3, astar.imag, "a* = 1/(1 − k)", fontsize=11)
ax.set_aspect("equal"); ax.set_title(f"B. a ↦ 1 + k·a with k = e^(2πi/{N}): a clock hand about a*,\nback after {N} ticks (same cycle; resemblance, not derivation)")
save(fig, "pb7_clocks.png")

# ---------------------------------------------------------------- pb8: mirrors
fig, axs = plt.subplots(1, 2, figsize=(17, 6)); fig.suptitle("pb8 (DRAFT). Two ways a mirror shapes space: Hořava–Witten walls and the Klein glide", fontsize=15); fig.subplots_adjust(top=0.8)
ax = axs[0]; ax.axis("off"); ax.set_xlim(-1.5, 1.5); ax.set_ylim(-1.5, 1.5); ax.set_aspect("equal")
t = np.linspace(0, 2 * np.pi, 200); ax.plot(np.cos(t), np.sin(t), color=GREY)
ax.plot([0, 0], [-1.35, 1.35], color="k", ls=":"); ax.plot([0, 0], [1, -1], "o", ms=16, color=RED)
ax.text(0.12, 1.08, "E8 wall (fixed by the mirror)", color=RED, fontsize=12); ax.text(0.12, -1.18, "E8 wall (fixed by the mirror)", color=RED, fontsize=12)
ax.annotate("", xy=(-0.7, 0.5), xytext=(0.7, 0.5), arrowprops=dict(arrowstyle="<->", color=GREY)); ax.text(-0.35, 0.58, "mirror pairs", color=GREY)
ax.set_title("A. M-theory on a circle folded by a mirror (S¹/Z2):\nthe mirror's two fixed points are walls carrying E8 each\n(Hořava–Witten 1996)")
ax = axs[1]; ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(0, 6)
for i in range(5):
    y = 0.7 + i * 1.1; ax.text(2 + (i % 2) * 0.0, y, "R" if i % 2 == 0 else "Я", fontsize=26, color=ORANGE if i % 2 else BLUE, ha="center")
    ax.annotate("", xy=(2, y + 0.95), xytext=(2, y + 0.45), arrowprops=dict(arrowstyle="->", color=GREY))
ax.text(3.2, 3.0, "Klein glide: each step mirrors, with no fixed point,\nso no walls; but it flips the hand of matter,\nwhich the weak force forbids (klein-universe U4).\n\nThe law itself reads the same in the mirror\n(round 32): the mirror is a symmetry of the law,\nnot of the matter.", fontsize=11.5)
ax.set_title("B. a free mirror (Klein space): no walls, hand flipped")
save(fig, "pb8_mirrors.png")
print("figures written")
