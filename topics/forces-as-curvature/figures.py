"""DRAFT pictures for forces-as-curvature. Run after checks.py. Orange = known, blue = ours."""
import json, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle

OR, BL, GR, RD = "#e07b00", "#1f6fc2", "#777777", "#c0392b"
R = json.load(open("results.json"))
plt.rcParams.update({"font.size": 11, "figure.dpi": 130})


def save(fig, name):
    fig.tight_layout()
    fig.savefig(f"figures/{name}.png", bbox_inches="tight")
    plt.close(fig)

# f1: one recipe ------------------------------------------------------------
fig, ax = plt.subplots(1, 2, figsize=(12, 4.8), gridspec_kw={"width_ratios": [1, 1.6]})
a = ax[0]; a.set_xlim(-0.3, 1.5); a.set_ylim(-0.4, 1.4); a.axis("off"); a.set_aspect("equal")
pts = [(0, 0), (1, 0), (1, 1), (0, 1), (0, 0)]
labs = ["a₁", "a₂", "ā₃", "ā₄"]
for i in range(4):
    (x0, y0), (x1, y1) = pts[i], pts[i + 1]
    a.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=18, color=BL, lw=2.5))
    a.text((x0 + x1) / 2 + (0.1 if x0 == x1 else 0) * (1 if x0 > 0.5 else -1.6), (y0 + y1) / 2 + (0.08 if y0 == y1 and y0 > 0.5 else -0.12 if y0 == y1 else 0),
           labs[i], color=BL, fontsize=14, ha="center")
a.add_patch(Rectangle((0.15, 0.15), 0.7, 0.7, color=OR, alpha=0.15))
a.text(0.5, 0.5, "leftover\n= curvature\n× area", ha="center", va="center", color=OR, fontsize=12)
a.text(0.5, -0.33, "a₁ ∘ₖ a₂ ∘ₖ ā₃ ∘ₖ ā₄", ha="center", color=BL, fontsize=13)
a.set_title("One recipe: chain the link values round a small loop")
b = ax[1]; b.axis("off")
rows = [["Force", "What the link values are", "k·ab term", "Leftover read as"],
        ["Electromagnetism", "phases (1×1, U(1)), k = i·g", "cancels (they commute)", "F = dA"],
        ["Weak", "2×2 matrices, SU(2)", "kept: [A, A]", "W field strength"],
        ["Strong", "3×3 matrices, SU(3)", "kept: [A, A], size C_A = 3", "gluon field strength"],
        ["Gravity (curved)", "4×4 Lorentz turns", "kept: ω∧ω", "Riemann curvature"],
        ["Gravity (twisted)", "5×5 turn + shift", "turn part zero", "torsion T = de"]]
tb = b.table(cellText=rows[1:], colLabels=rows[0], loc="center", cellLoc="left", colWidths=[0.2, 0.33, 0.27, 0.22])
tb.auto_set_font_size(False); tb.set_fontsize(10); tb.scale(1, 1.9)
for (i, j), c in tb.get_celld().items():
    if i == 0:
        c.set_text_props(weight="bold")
    if j == 2 and i > 0:
        c.set_text_props(color=BL)
b.set_title("Same chain, different ring (known: Weyl, Yang–Mills, Utiyama, Kibble, Cartan)")
save(fig, "f1_one_recipe")

# f2: plaquettes -------------------------------------------------------------
fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
cols = {"U(1) electromagnetism": GR, "SU(2) weak": BL, "SU(3) strong": "#2e8b57"}
for name, d in R["gauge"].items():
    rows = np.array(d["rows"])
    ax[0].loglog(rows[:, 0], rows[:, 1], "o-", color=cols[name], label=name)
ax[0].loglog([0.01, 0.08], [2e-4 * 1, 2e-4 * 64], "k:", lw=1)
ax[0].text(0.03, 3e-3, "slope 2", fontsize=9)
ax[0].set_xlabel("loop size ε"); ax[0].set_ylabel("|leftover/ε² − (dA − dA + k[A,A])|")
ax[0].set_title("A. With the k·ab term: the loop reads the field exactly"); ax[0].legend(fontsize=9)
names_ = list(R["gauge"].keys())
miss = [R["gauge"][n]["rows"][-1][2] for n in names_]
Fs = [R["gauge"][n]["F_size"] for n in names_]
x = np.arange(3)
ax[1].bar(x - 0.18, Fs, 0.36, color=OR, label="size of the field")
ax[1].bar(x + 0.18, miss, 0.36, color=BL, label="miss if the k·ab term is dropped")
ax[1].set_xticks(x); ax[1].set_xticklabels(["U(1)\nlight", "SU(2)\nweak", "SU(3)\nstrong"])
ax[1].set_title("B. Drop k·ab: light is untouched, weak and strong are wrong")
ax[1].legend(fontsize=9)
save(fig, "f2_plaquettes")

# f3: why counting stopped at 1+1 ---------------------------------------------
fig, ax = plt.subplots(figsize=(7.5, 4.4))
dims = ["2D", "3D", "4D"]
riem = [1, 6, 20]; ric = [1, 6, 10]; wey = [0, 0, 10]
x = np.arange(3)
ax.bar(x, ric, 0.5, color=OR, label="Ricci part (what a volume count reads)")
ax.bar(x, wey, 0.5, bottom=ric, color=BL, label="Weyl part (only loops read it)")
for i in range(3):
    ax.text(i, riem[i] + 0.4, f"{riem[i]} numbers", ha="center")
ax.set_xticks(x); ax.set_xticklabels(dims); ax.set_ylabel("independent curvature numbers at a point")
ax.set_title("Why counting stopped: outside a mass Ricci = 0,\nso the whole pull sits in the Weyl part")
ax.legend(fontsize=9, loc="upper left"); ax.set_ylim(0, 24)
save(fig, "f3_why_counting_stopped")

# f4: Schwarzschild from loops ------------------------------------------------
fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
sec = R["schwarzschild_sectional"]
rr = np.linspace(2.4, 21, 200)
style = {"t-r": (RD, 2), "theta-phi": ("#8e44ad", 2), "t-theta": (BL, -1), "t-phi": ("#2e8b57", -1), "r-theta": (OR, -1), "r-phi": (GR, -1)}
for nm, (c, m) in style.items():
    ax[0].plot(rr, m / rr**3, color=c, lw=1.2, alpha=0.7)
    d = np.array(sec[nm])
    ax[0].plot(d[:, 0], d[:, 1], "o", color=c, ms=6, label=f"{nm}: {'+2' if m > 0 else '−1'}·M/r³")
ax[0].set_xscale("log"); ax[0].set_yscale("symlog", linthresh=1e-4)
ax[0].set_xlabel("r (units of M)"); ax[0].set_ylabel("curvature of the plane K")
ax[0].set_title("A. Six planes through a point: dots from loops,\nlines exact (Schwarzschild 1916)"); ax[0].legend(fontsize=8)
ax[0].text(2.6, 2e-4, "two planes on the upper line,\nfour on the lower", fontsize=8)
for planes, c, lab in ((("t-r", "t-theta", "t-phi"), BL, "planes through time"), (("t-r", "r-theta", "r-phi"), OR, "planes through r")):
    tot = [sum(np.array(sec[p])[i, 1] for p in planes) for i in range(len(sec["t-r"]))]
    ax[1].plot(np.array(sec["t-r"])[:, 0], np.abs(tot) + 1e-16, "o-", color=c, label=lab)
ax[1].set_yscale("log"); ax[1].set_xlabel("r (units of M)"); ax[1].set_ylabel("|sum of the three curvatures|")
ax[1].set_title("B. They cancel in threes: Ricci = 0 (vacuum),\nso a volume count is blind to M here"); ax[1].legend(fontsize=9)
save(fig, "f4_schwarzschild_from_loops")

# f5: Gravity Probe B ----------------------------------------------------------
fig, ax = plt.subplots(1, 2, figsize=(12, 4.4))
conv = np.array(R["gpb"]["convergence"])
ax[0].loglog(conv[:, 0], np.abs(conv[:, 1] - R["gpb"]["rate"]), "o-", color=BL)
ax[0].set_xlabel("steps n in the chain round one orbit"); ax[0].set_ylabel("error (mas/yr)")
ax[0].set_title("A. 2^m identical Lorentz steps: error falls like 1/n²")
m, e = R["gpb"]["measured"]
ax[1].errorbar([0], [m], yerr=[e], fmt="s", color=OR, capsize=6, ms=8, label="Gravity Probe B measured (Everitt 2011)")
ax[1].plot([1], [R["gpb"]["rate"]], "o", color=BL, ms=9, label=f"one 4D chain: {R['gpb']['rate']:.1f}")
ax[1].plot([2], [R["gpb"]["GR_with_oblateness"]], "^", color=GR, ms=8, label="GR with Earth's flattening: 6606.1")
ax[1].set_xlim(-0.6, 2.6); ax[1].set_xticks([]); ax[1].set_ylabel("geodetic drift (milliarcsec / yr)")
ax[1].set_title("B. Gyroscope turned by the chain vs measurement"); ax[1].legend(fontsize=9, loc="lower right")
ax[1].set_ylim(6570, 6625)
save(fig, "f5_gravity_probe_b")

# f6: colour factor ------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7.5, 4.3))
meas = R["casimir"]["measured"]
for i, (ex, v) in enumerate(meas.items()):
    ax.errorbar([i], [v[0]], yerr=[v[1]], fmt="s", color=OR, capsize=6, ms=8)
    ax.text(i + 0.08, v[0] + 0.05, ex, color=OR)
ax.axhline(3, color=BL, lw=2, label="SU(3) commutators: C_A = 3")
ax.axhline(0, color=GR, lw=2, ls="--", label="no k·ab term (abelian): C_A = 0")
ax.set_xlim(-0.6, 1.8); ax.set_xticks([]); ax.set_ylabel("gluon self-charge C_A")
ax.set_title("The k·ab term's size, measured in four-jet events at LEP")
ax.legend(fontsize=9, loc="center right")
save(fig, "f6_colour_factor")

# f7: curved vs twisted ---------------------------------------------------------
fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
tw = R["twist"]; planes = list(tw.keys()); x = np.arange(len(planes))
ax[0].bar(x - 0.2, [tw[p]["LC_turn"] for p in planes], 0.4, color=OR, label="turn left over (curvature)")
ax[0].bar(x + 0.2, [tw[p]["LC_shift"] + 1e-12 for p in planes], 0.4, color=BL, label="shift left over (twist)")
ax[0].set_title("A. GR's connection: loops come back turned, never shifted")
ax[1].bar(x - 0.2, [tw[p]["W_turn"] + 1e-12 for p in planes], 0.4, color=OR, label="turn left over (curvature)")
ax[1].bar(x + 0.2, [tw[p]["W_shift"] for p in planes], 0.4, color=BL, label="shift left over (twist)")
ax[1].set_title("B. Teleparallel connection, same spacetime:\nnever turned, only shifted")
for a in ax:
    a.set_yscale("log"); a.set_ylim(1e-12, 3); a.set_xticks(x); a.set_xticklabels(planes, rotation=30, fontsize=9)
    a.set_ylabel("size of leftover / area"); a.legend(fontsize=8, loc="lower left")
    a.axhline(1e-6, color="k", ls=":", lw=1); a.text(5.45, 2e-6, "numerical noise below", fontsize=8, ha="right")
save(fig, "f7_curved_or_twisted")

# f8: Aharonov-Bohm -------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7.5, 4.2))
fl = np.linspace(0, 2.2, 300)
ax.plot(fl, np.mod(2 * np.pi * fl + np.pi, 2 * np.pi) - np.pi, color=BL, label="U(1) chain round the solenoid (any path)")
ax.plot([0.5], [np.pi], "s", color=OR, ms=9, label="Tonomura et al. 1986: π at flux h/2e")
ax.plot([0, 1], [0, 0], "s", color=OR, ms=9, mfc="white")
ax.set_xlabel("enclosed flux (units of h/e)"); ax.set_ylabel("electron phase left over (rad)")
ax.set_title("Field zero on the path, phase still left over: the loop sees curvature inside")
ax.legend(fontsize=9, loc="lower right")
save(fig, "f8_aharonov_bohm")
print("figures written")

# f9: frame dragging as the twist of the time leg ---------------------------------
fd = R["frame_drag_twist"]
prof = np.array(fd["orbit_profile"])
fig, ax = plt.subplots(1, 2, figsize=(12, 4.4))
ax[0].plot(np.degrees(prof[:, 0]), prof[:, 2], color=BL, label="drift about the spin axis")
ax[0].plot(np.degrees(prof[:, 0]), prof[:, 1], color=GR, label="sideways drift (averages to 0)")
ax[0].axhline(fd["mas_per_yr"], color=BL, ls=":", label=f"orbit average {fd['mas_per_yr']:.1f}")
ax[0].set_xlabel("position round the polar orbit (degrees from the north pole)"); ax[0].set_ylabel("mas / yr")
ax[0].set_title("A. Gyroscope drift from the twist of the time leg, de⁰"); ax[0].legend(fontsize=8)
m, e = fd["measured"]
ax[1].errorbar([0], [m], yerr=[e], fmt="s", color=OR, capsize=6, ms=8, label="Gravity Probe B measured (Everitt 2011)")
ax[1].plot([1], [fd["mas_per_yr"]], "o", color=BL, ms=9, label=f"time-leg twist, polar average: {fd['mas_per_yr']:.1f}")
ax[1].plot([2], [fd["GPB_GR"]], "^", color=GR, ms=8, label="GP-B's own GR figure: 39.2")
ax[1].set_xlim(-0.6, 2.6); ax[1].set_xticks([]); ax[1].set_ylabel("frame dragging (mas / yr)")
ax[1].set_title("B. Frame dragging: twist vs measurement"); ax[1].legend(fontsize=8, loc="lower right"); ax[1].set_ylim(20, 50)
save(fig, "f9_frame_drag_twist")
print("f9 written")

# f10: Kaluza-Klein, one 5D chain ---------------------------------------------------
kk = R["kk"]
fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
pl = list(kk["planes"].keys()); vals = [kk["planes"][p] for p in pl]
cols_ = [BL if "w" in p else GR for p in pl]
ax[0].bar(range(len(pl)), vals, color=cols_)
ax[0].set_xticks(range(len(pl))); ax[0].set_xticklabels(pl, rotation=40, fontsize=9)
ax[0].set_ylabel("size of loop leftover / area")
ax[0].set_title("A. Flat spacetime + a magnetic field, as one 5D shape:\nloops through w (blue) carry the field; x–y (grey) bends by 3B²/4")
tr = np.array(kk["traj"])
ax[1].plot(tr[:, 0], tr[:, 1], color=BL, lw=2, label="straight line (geodesic) in 5D, seen in x–y")
th_ = np.linspace(0, 2 * np.pi, 200); cx = (tr[:, 0].max() + tr[:, 0].min()) / 2; cy = (tr[:, 1].max() + tr[:, 1].min()) / 2
ax[1].plot(cx + kk["radius"] * np.cos(th_), cy + kk["radius"] * np.sin(th_), ":", color=OR, label="Lorentz-force circle, radius v/(qB/m)")
ax[1].set_aspect("equal"); ax[1].legend(fontsize=8, loc="lower center")
ax[1].set_title(f"B. 5D curvature scalar from loops = −B²/2 ({kk['R5_loops']:.4f});\nfree fall in 5D = a charge circling in 4D")
save(fig, "f10_kaluza_klein")
print("f10 written")

# f11: SU(2) from a hidden 3-sphere --------------------------------------------------
nk = R["nonabelian_kk"]
fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
sw = np.array(nk["sweep"])
ff = np.linspace(0, sw[:, 0].max() * 1.05, 100)
ax[0].plot(ff, -(1.3**2 / 8) * ff**2, color=OR, label="Kaluza–Klein formula −(ρ²/16) FᵃFᵃ")
ax[0].plot(sw[:, 0], sw[:, 1], "o", color=BL, ms=8, label="hidden 3-sphere: from loops")
ax[0].plot(sw[:, 0], sw[:, 2], "s", color=GR, ms=7, label="hidden flat torus: from loops (always 0)")
ax[0].set_xlabel("c₁·c₂ (a constant field: no slopes, only the k·ab part)")
ax[0].set_ylabel("7D curvature minus the sphere's own")
ax[0].set_title("A. A field made only of the commutator term\nexists only if the hidden space's symmetries twist")
ax[0].legend(fontsize=8)
ps_ = nk["planes_sphere"]; pt_ = nk["planes_torus"]
lab7 = ["t", "x", "y", "z", "φ", "θ", "ψ"]
keys = list(ps_.keys())
labs_ = [lab7[int(k.split("-")[0])] + lab7[int(k.split("-")[1])] for k in keys]
xx = np.arange(len(keys))
ax[1].bar(xx - 0.2, [ps_[k] for k in keys], 0.4, color=BL, label="hidden 3-sphere")
ax[1].bar(xx + 0.2, [pt_[k] + 1e-12 for k in keys], 0.4, color=GR, label="hidden torus")
ax[1].set_yscale("log"); ax[1].set_ylim(1e-12, 10)
ax[1].set_xticks(xx); ax[1].set_xticklabels(labs_, rotation=90, fontsize=8)
ax[1].set_ylabel("loop leftover / area"); ax[1].legend(fontsize=8)
ax[1].set_title("B. Same constant field, 21 loop planes of one 7×7 chain")
save(fig, "f11_su2_hidden_sphere")
print("f11 written")
