"""DRAFT figures for the Klein-universe round. Reuses checks.py's objects (re-running it is deterministic, about 2 s).
Colour key: orange = the Klein glide (flips the hand); blue = slides and the half-turn (keep the hand); green = allowed;
red = ruled out; grey = reference."""
import runpy, io, contextlib, json
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle, FancyArrowPatch

with contextlib.redirect_stdout(io.StringIO()):
    g = runpy.run_path("checks.py")
R = json.load(open("results.json"))
OR, BL, GR, RD, GY = "#e07b00", "#1f6fb4", "#2a9d4b", "#c62828", "#777777"
def draft(fig): fig.text(0.995, 0.005, "DRAFT · Nov 27 2026 · Affine Composition Theory", ha="right", va="bottom", fontsize=8, color=GY)

F = np.array([[0, 0], [0.18, 0], [0.18, 0.5], [0.45, 0.5], [0.45, 0.66], [0.18, 0.66], [0.18, 0.84], [0.55, 0.84], [0.55, 1.0], [0, 1.0]])
F = F - [0.27, 0.5]
def glyph(ax, cx, cy, s, mirror=False, rot=False, color="k"):
    P = F.copy()
    if mirror: P[:, 0] *= -1
    if rot: P *= -1
    ax.add_patch(Polygon(P * s + [cx, cy], closed=True, fc=color, ec="none", alpha=0.9))

# ------------------------------------------------------------------ ku1 the gluing
fig, axs = plt.subplots(1, 3, figsize=(15, 5.2), gridspec_kw=dict(width_ratios=[1, 1, 1.5]))
for ax, title, kw, col in ((axs[0], "A. Klein space (E7): the glide\nslide up by c and mirror left-right", dict(mirror=True), OR),
                           (axs[1], "B. Half-turn space (E2): the screw\nslide up by c and turn 180°", dict(rot=True), BL)):
    for j in range(3):
        y = j * 1.25
        ax.add_patch(Rectangle((0, y), 1.6, 1.0, fc="#f4f4f4", ec=GY))
        glyph(ax, 0.8, y + 0.5, 0.75, color=col if j % 2 else "k", **(kw if j % 2 else {}))
        ax.text(1.7, y + 0.5, ["layer 0", "layer c", "layer 2c"][j], va="center", fontsize=9, color=GY)
    for j in range(2):
        ax.add_patch(FancyArrowPatch((-0.15, j * 1.25 + 0.5), (-0.15, (j + 1) * 1.25 + 0.5), arrowstyle="-|>", mutation_scale=14, color=col,
                                     connectionstyle="arc3,rad=0.35"))
    ax.set_xlim(-0.6, 2.3); ax.set_ylim(-0.2, 3.7); ax.set_aspect("equal"); ax.axis("off"); ax.set_title(title, fontsize=11)
axs[0].text(0.8, -0.15, "two steps up: back to the start,\nbut one step up: the hand is flipped", ha="center", va="top", fontsize=9, color=OR)
axs[1].text(0.8, -0.15, "every step keeps the hand:\ntwisted, but orientable", ha="center", va="top", fontsize=9, color=BL)
ax = axs[2]; a, c = 1.0, 0.6
for i in range(-2, 3):
    for j in range(-2, 3):
        ax.add_patch(Rectangle((i * a, j * c), a, c, fc="#fff3e0" if j % 2 else "#f4f4f4", ec="#bbbbbb", lw=0.6))
        glyph(ax, i * a + a / 2, j * c + c / 2, 0.38, mirror=bool(j % 2), color=OR if j % 2 else "k")
for i in range(-4, 5):
    ax.axvline(i * a / 2, color=OR, lw=0.8, ls=":", alpha=0.7)
ax.add_patch(Rectangle((0, 0), a, c, fill=False, ec="k", lw=2))
x0, z0 = 0.14, 0.1
clones = [(x0, z0, "us")]
for (sx, sz) in ((-1, 1), (-1, -1)):
    for n in (0, 1):
        clones.append((-x0 + n * a, z0 + sz * c, "glide copy"))
for n in (-1, 1): clones.append((x0 + n * a, z0, "slide copy"))
clones.append((x0, z0 + 2 * c, "slide copy")); clones.append((x0, z0 - 2 * c, "slide copy"))
for (x, z, lab) in clones:
    ax.plot(x, z, "o", ms=8 if lab == "us" else 6, color="k" if lab == "us" else (OR if "glide" in lab else BL), zorder=5)
ax.annotate("us", (x0, z0), (x0 + 0.1, z0 + 0.12), fontsize=10, fontweight="bold")
ax.annotate("nearest glide copy:\ndistance √(4d² + c²)", (-x0, z0 + c), (-1.95, 1.05), fontsize=9, color=OR, bbox=dict(fc="white", ec="none", alpha=0.9), arrowprops=dict(arrowstyle="->", color=OR))
ax.annotate("mirror planes\n(dotted, every a/2)", (a / 2, -1.2), (0.9, -1.45), fontsize=9, color=OR, arrowprops=dict(arrowstyle="->", color=OR))
ax.annotate("", (0, z0 - 0.02), (x0, z0 - 0.02), arrowprops=dict(arrowstyle="<->", color="k"))
ax.text(x0 / 2, z0 - 0.1, "d", ha="center", fontsize=9)
ax.set_xlim(-2, 2); ax.set_ylim(-1.5, 1.5); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
ax.set_title("C. A slice through Klein space: every row up is mirrored\nblack box = one cell (a by c); dots = copies of us", fontsize=11)
fig.suptitle("ku1 · Two ways to twist flat space: the Klein glide flips the hand, the half-turn keeps it (U1)", fontsize=13)
draft(fig); fig.tight_layout(rect=(0, 0.02, 1, 0.95)); fig.savefig("figures/ku1_gluing.png", dpi=140); plt.close(fig)

# ------------------------------------------------------------------ ku2 the sky
O, CHI = g["O"], g["CHI"]; fK = g["fK"]; GT, TA = g["GT"], g["TA"]; ginv, app, basis = g["ginv"], g["app"], g["basis"]
nl, nb = 360, 180
lon = np.linspace(-np.pi, np.pi, nl); lat = np.linspace(-np.pi / 2, np.pi / 2, nb)
LO, LA = np.meshgrid(lon, lat)
dirs = np.stack([np.cos(LA) * np.cos(LO), np.cos(LA) * np.sin(LO), np.sin(LA)], -1).reshape(-1, 3)
T = fK(O + CHI * dirs).reshape(nb, nl)
fig = plt.figure(figsize=(15, 8.4))
ax = fig.add_subplot(2, 2, (1, 2), projection="mollweide")
ax.pcolormesh(LO, LA, T, cmap="RdBu_r", shading="auto", alpha=0.75)
def circ_pts(gg, N=400):
    gO = app(gg, O); s = np.linalg.norm(gO - O); al = np.arccos(s / (2 * CHI)); n1 = (gO - O) / s; e1, e2 = basis(n1)
    phi = np.linspace(0, 2 * np.pi, N, endpoint=False)
    P = O + CHI * (np.cos(al) * n1 + np.sin(al) * (np.outer(np.cos(phi), e1) + np.outer(np.sin(phi), e2)))
    gi = ginv(gg); Q = (gi[0] @ P.T).T + gi[1]
    return P, Q
def ll(X):
    d = (X - O) / CHI; return np.arctan2(d[:, 1], d[:, 0]), np.arcsin(np.clip(d[:, 2], -1, 1))
for gg, col, name in ((GT, OR, "glide pair"), (TA, BL, "slide pair")):
    P, Q = circ_pts(gg)
    for X, ls in ((P, "-"), (Q, "--")):
        lo, la = ll(X)
        brk = np.where(np.abs(np.diff(lo)) > np.pi)[0] + 1
        for seg in np.split(np.arange(len(lo)), brk):
            ax.plot(lo[seg], la[seg], ls, color=col, lw=2)
    for k, idx in enumerate((0, 100, 200, 300)):
        for X in (P, Q):
            lo, la = ll(X[idx:idx + 1])
            ax.text(lo[0], la[0], str(k + 1), color="white", fontsize=9, ha="center", va="center", fontweight="bold",
                    bbox=dict(boxstyle="circle,pad=0.15", fc=col, ec="none"))
ax.set_xticklabels([]); ax.set_yticklabels([])
ax.set_title("A. The whole sky of a small toy Klein universe, seen from us. Matched circles: orange = glide pair, blue = slide pair.\n"
             "Points with the same number show the same temperature. Glide pair: 1→2→3→4 run the same way round both circles. "
             "Slide pair: opposite ways.", fontsize=10)
pr = R["U2"]["pairs"]
for pos, gg, key, col, good in ((3, GT, "Klein glide", OR, +1), (4, TA, "Klein slide a", BL, -1)):
    ax = fig.add_subplot(2, 2, pos)
    P, Q = circ_pts(gg, 720)
    gi = ginv(gg); n2 = (app(gi, O) - O) / np.linalg.norm(app(gi, O) - O); f1, f2 = basis(n2)
    s = np.linalg.norm(app(gg, O) - O); al = np.arccos(s / (2 * CHI)); phi = np.linspace(0, 2 * np.pi, 720, endpoint=False)
    C2 = O + CHI * (np.cos(al) * n2 + np.sin(al) * (np.outer(np.cos(phi), f1) + np.outer(np.sin(phi), f2)))
    T1, T2 = fK(P), fK(C2)
    best = {}
    for sg in (+1, -1):
        sc = []
        for d in range(720):
            T2r = T2[(sg * np.arange(720) + d) % 720]; sc.append(2 * np.mean(T1 * T2r) / (np.mean(T1**2) + np.mean(T2r**2)))
        d = int(np.argmax(sc)); best[sg] = T2[(sg * np.arange(720) + d) % 720]
    deg = np.degrees(phi)
    ax.plot(deg, T1, color="k", lw=2, label="circle 1")
    ax.plot(deg, best[good], color=col, lw=1.6, ls="--", label=f"circle 2 read {'the same way' if good > 0 else 'the opposite way'} round: S = 1.00")
    ax.plot(deg, best[-good], color=GY, lw=1, alpha=0.8,
            label=f"circle 2 read {'the opposite way' if good > 0 else 'the same way'} round, best shift: S = {pr[key]['S_opposite' if good > 0 else 'S_same']:.2f}")
    ax.set_xlim(0, 360); ax.set_xlabel("angle round the circle (degrees)"); ax.set_ylabel("temperature (toy units)")
    ax.legend(fontsize=8, loc="upper right")
    ax.set_title(("B. Glide pair (flipped): matches only when both circles are read the same way round"
                  if good > 0 else "C. Slide pair (ordinary): matches only when read opposite ways round"), fontsize=10, color=col)
fig.suptitle("ku2 · What a Klein-shaped universe would show: matched circle pairs, and the glide's pairs run the same way round (U2)", fontsize=13)
draft(fig); fig.tight_layout(rect=(0, 0.02, 1, 0.95)); fig.savefig("figures/ku2_sky_circles.png", dpi=140); plt.close(fig)

# ------------------------------------------------------------------ ku3 the size window
D = R["U0"]["lss_diameter_Gpc"]; af = g["allowed_fraction"]
fig, axs = plt.subplots(1, 3, figsize=(19, 5.8), gridspec_kw=dict(width_ratios=[1.2, 1.3, 1]))
ax = axs[0]
A_ = np.linspace(0.3, 2.6, 300); C_ = np.linspace(0.2, 1.4, 300)
Z = np.array([[af(a, c) for a in A_] for c in C_])
m = ax.pcolormesh(A_, C_, np.ma.masked_where(Z <= 0, Z), cmap="Greens", vmin=0, vmax=1, shading="auto")
ax.pcolormesh(A_, C_, np.ma.masked_where(Z > 0, np.ones_like(Z)), cmap="Reds", vmin=0, vmax=2, shading="auto", alpha=0.5)
plt.colorbar(m, ax=ax, label="share of positions in the universe where we could sit")
cm, am = R["U3"]["shortest_loop_min_D"], R["U3"]["a_needed_at_shortest_D"]
ax.plot(am, cm, "k*", ms=14); ax.annotate(f"shortest allowed loop:\nc = {cm:.3f} D = {cm*D:.1f} Gpc\n(needs a ≥ {am:.2f} D)", (am, cm), (1.05, 0.24), fontsize=9, bbox=dict(fc="white", ec="none", alpha=0.85),
                                         arrowprops=dict(arrowstyle="->"))
ax.axvline(0.985, color=RD, lw=1); ax.axhline(0.985 / 2, color=RD, lw=1)
ax.text(0.33, 1.3, "red: ruled out for every\nposition (circles would\nalready have been seen)", fontsize=9, color=RD, va="top")
ax.set_xlabel("a, side of the cell across the mirrors (units of D)"); ax.set_ylabel("c, length of the glide (units of D)")
ax.set_title(f"A. Which Klein universes survive (b large)\nD = last-scattering diameter = {D:.2f} Gpc", fontsize=11)
ax = axs[1]
items = [("shortest Klein loop still allowed", cm * D, OR), ("radius of last-scattering sphere", D / 2, GY),
         ("radius of observable universe", R["U0"]["particle_horizon_Gpc"], GY), ("circle bound: every copy of us beyond", 0.985 * D, RD),
         ("diameter of last-scattering sphere", D, GY), ("diameter of observable universe", R["U0"]["observable_diameter_Gpc"], GY),
         ("future tests can still see it up to (COMPACT IIb)", 1.2 * D, GR)]
for i, (lab, v, col) in enumerate(items):
    ax.barh(i, v, color=col, alpha=0.8); ax.text(v + 0.4, i, f"{v:.1f} Gpc", va="center", fontsize=9)
ax.set_yticks(range(len(items))); ax.set_yticklabels([x[0] for x in items], fontsize=9); ax.invert_yaxis()
ax.axvspan(0.985 * D, 1.2 * D, color=GR, alpha=0.12); ax.text(30.3, 2.5, "allowed and testable", ha="center", va="center", rotation=90, fontsize=9, color=GR)
ax.set_xlim(0, 40); ax.set_xlabel("comoving distance (Gpc)")
ax.set_title("B. The numbers (Planck 2018 distances)", fontsize=11)
ax = axs[2]; a, c = 2.0, 0.6
dgrid = np.linspace(0, a, 801)
dm = np.abs(dgrid - np.round(dgrid / (a / 2)) * (a / 2))
near = np.minimum.reduce([np.full_like(dgrid, a), np.full_like(dgrid, 2 * c), np.sqrt(4 * dm**2 + c**2)])
ax.plot(dgrid, near, color=OR, lw=2, label="distance to nearest copy of us")
ax.axhline(0.985, color=RD, ls="--", label="circle bound 0.985 D")
ax.fill_between(dgrid, 0, 2.2, where=near > 0.985, color=GR, alpha=0.15, label="allowed places to sit")
for k in range(5): ax.axvline(k * a / 2, color=OR, ls=":", lw=1)
ax.set_xlim(0, a); ax.set_ylim(0, 2.2); ax.set_xlabel("where we sit across the cell (units of D); dotted = mirror planes")
ax.set_ylabel("units of D"); ax.legend(fontsize=8, loc="upper right")
ax.set_title(f"C. Example a = 2 D, c = 0.6 D: we must sit away\nfrom the mirrors ({100*af(a, c):.0f}% of places allowed)", fontsize=11)
fig.suptitle("ku3 · A Klein-shaped universe is not ruled out by the sky: its shortest loop can be about half the last-scattering diameter (U0, U3)", fontsize=13)
draft(fig); fig.tight_layout(rect=(0, 0.02, 1, 0.94)); fig.savefig("figures/ku3_size_window.png", dpi=140); plt.close(fig)

# ------------------------------------------------------------------ ku4 handedness, walls and torsion
fig, axs = plt.subplots(1, 3, figsize=(16, 5.2), gridspec_kw=dict(width_ratios=[1.3, 1, 1]))
ax = axs[0]
t = np.linspace(0, 4 * np.pi, 400)
for x0, hand, col, lab in ((0.0, +1, BL, "left-handed\n(as the weak force needs)"), (3.2, -1, OR, "comes back\nright-handed")):
    ax.plot(x0 + 0.45 * np.cos(t), t / (4 * np.pi) * 3 + 0.3 * hand * np.sin(t) * 0, color=col, lw=0)
    xs = x0 + 0.45 * np.cos(hand * t); ys = t / (4 * np.pi) * 3; zs = np.sin(hand * t)
    for i in range(len(t) - 1):
        ax.plot(xs[i:i + 2], ys[i:i + 2], color=col, lw=2.5 if zs[i] > 0 else 1, alpha=1 if zs[i] > 0 else 0.35)
    ax.annotate("", (x0, 3.35), (x0, 3.0), arrowprops=dict(arrowstyle="-|>", color=col, lw=2))
    ax.text(x0, -0.45, lab, ha="center", va="top", fontsize=10, color=col)
ax.add_patch(FancyArrowPatch((0.7, 1.5), (2.5, 1.5), arrowstyle="-|>", mutation_scale=18, color="k", connectionstyle="arc3,rad=-0.5"))
ax.text(1.6, 2.35, "once round the glide loop\nS γ5 S⁻¹ = −γ5 (exact)", ha="center", fontsize=10)
ax.text(1.6, 0.55, "a turn (half-turn space)\nkeeps the hand: [S, γ5] = 0", ha="center", fontsize=9, color=BL)
ax.set_xlim(-0.9, 4.1); ax.set_ylim(-1.3, 3.6); ax.axis("off")
ax.set_title(f"A. Why the weak force objects (U4)\nthe Standard Model has no mirror symmetry: J = {R['U4']['jarlskog']:.2e} ≠ 0", fontsize=11)
ax = axs[1]
vals = [("hand chosen at the weak scale\n(wall tension ≈ v³)", R["U5"]["sigma_weak_GeV3"], RD),
        ("most a horizon wall can carry\n(below critical density)", R["U5"]["sigma_max_energy_GeV3"], GY),
        ("most the smooth CMB allows\n(δT/T < 10⁻⁵)", R["U5"]["sigma_max_cmb_GeV3"], GR)]
for i, (lab, v, col) in enumerate(vals):
    ax.barh(i, np.log10(v) + 12, left=-12, color=col, alpha=0.8); ax.text(np.log10(v) + 0.3, i, f"{v:.0e} GeV³", va="center", fontsize=9)
ax.set_yticks(range(3)); ax.set_yticklabels([v[0] for v in vals], fontsize=9); ax.invert_yaxis()
ax.set_xlim(-12, 12); ax.set_xlabel("log₁₀ wall tension (GeV³)")
ax.set_title(f"B. The loophole fails (U5): a forced mirror wall\nwould be ~{R['U5']['ratio_weak_over_cmb']:.0e} times too heavy", fontsize=11)
ax = axs[2]
dens = [("water", 1e3), ("atomic nucleus", 2.3e17), ("neutron star core", 1e18), ("Cartan density (neutrons):\nwhere local torsion matters", R["U6"]["cartan_density_neutron_kg_m3"])]
for i, (lab, v) in enumerate(dens):
    ax.barh(i, np.log10(v), color=(GY if i < 3 else "#7b4fa0"), alpha=0.85); ax.text(np.log10(v) + 0.8, i, f"10^{np.log10(v):.0f}", va="center", fontsize=9)
ax.set_yticks(range(4)); ax.set_yticklabels([d[0] for d in dens], fontsize=9); ax.invert_yaxis(); ax.set_xlim(0, 70)
ax.set_xlabel("log₁₀ density (kg/m³)")
ax.set_title("C. Shape twist ≠ torsion (U6): the shape's twist has\nzero torsion everywhere; Einstein–Cartan torsion\nneeds densities far past any star", fontsize=11)
fig.suptitle("ku4 · The catch: a Klein-shaped space flips left and right, and the weak force can't survive that (U4, U5, U6)", fontsize=13)
draft(fig); fig.tight_layout(rect=(0, 0.02, 1, 0.93)); fig.savefig("figures/ku4_handedness.png", dpi=140); plt.close(fig)
print("wrote figures/ku1_gluing.png, ku2_sky_circles.png, ku3_size_window.png, ku4_handedness.png")

# ------------------------------------------------------------------ ku5 phases round a twist (for the clock question)
fig, axs = plt.subplots(1, 3, figsize=(15, 5.2))
cases = [("A. Klein glide, Pin+ spinors", [(1, "once round: +1 or −1\n(a sign, and the hand flips)")], 1, OR),
         ("B. Klein glide, Pin− spinors", [(1j, "once round: +i or −i\n(a quarter turn, and the hand flips)")], -1, OR),
         ("C. Half-turn screw (orientable)", [(1j, "once round: +i or −i\n(a quarter turn, hand kept)")], -1, BL)]
for ax, (title, once, twice, col) in zip(axs, cases):
    th = np.linspace(0, 2 * np.pi, 300); ax.plot(np.cos(th), np.sin(th), color=GY, lw=1)
    ax.plot(1, 0, "ko", ms=8); ax.text(1.08, 0.12, "start: 1", fontsize=9)
    z = once[0][0]
    for zz in (z, -z):
        ax.annotate("", (zz.real if isinstance(zz, complex) else zz, zz.imag if isinstance(zz, complex) else 0), (0, 0),
                    arrowprops=dict(arrowstyle="-|>", color=col, lw=2))
    ax.text(0, -1.45, once[0][1], ha="center", fontsize=9, color=col)
    ax.plot(twice, 0, "s", ms=12, mfc="none", mec="k", mew=2)
    ax.text(twice * 0.95 + (-0.55 if twice < 0 else -0.6), 0.18 if twice < 0 else -0.42, f"twice round: {'+1' if twice > 0 else '−1'}\n(forced, no choice)", fontsize=9)
    ax.set_xlim(-1.7, 1.7); ax.set_ylim(-1.8, 1.4); ax.set_aspect("equal"); ax.axis("off"); ax.set_title(title, fontsize=11)
fig.suptitle("ku5 · What going round a twist does to a spinor's phase (U7): Pin− gives a quarter turn per trip and a forced half turn per two trips.\n"
             "The glide moves only space, so cosmic time stays one single-valued clock.", fontsize=12)
draft(fig); fig.tight_layout(rect=(0, 0.02, 1, 0.9)); fig.savefig("figures/ku5_phases.png", dpi=140); plt.close(fig)
print("wrote figures/ku5_phases.png")
