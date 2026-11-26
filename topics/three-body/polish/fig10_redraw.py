"""Redraw Figure 10 (the unlisted loop types) from the polished orbits and the re-read words.
Run: python3 fig10_redraw.py <dir with tbpy.py + libtb.so>   (matplotlib)"""
import sys, os, json
sys.path.insert(0, sys.argv[1] if len(sys.argv) > 1 else "..")
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from tbpy import suvakov, traj, shape, special_points
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(os.path.dirname(HERE), "figures")
KNOWN, GREY, RED, GREEN = "#d9822b", "#888888", "#c53030", "#2f855a"
BODY = ["#2b6cb0", "#c53030", "#222222"]
plt.rcParams.update({"font.size": 11, "axes.titlesize": 12})
SP = special_points()
P = {r["id"]: r for r in json.load(open(os.path.join(HERE, "polished.json")))}
A = {o["id"]: o for o in json.load(open(os.path.join(HERE, "words_all.json")))}

def sphere(ax, alpha=0.10, view=(22, -60)):
    u, v = np.mgrid[0:2 * np.pi:36j, 0:np.pi:18j]
    ax.plot_wireframe(np.cos(u) * np.sin(v), np.sin(u) * np.sin(v), np.cos(v), color=GREY, alpha=alpha, lw=0.5)
    th = np.linspace(0, 2 * np.pi, 200); ax.plot(np.cos(th), np.sin(th), 0 * th, color=GREY, lw=0.8)
    for k, p in SP.items():
        if k.startswith("collision"): ax.scatter(*p, color=RED, marker="x", s=50, depthshade=False)
        elif k.startswith("Lagrange"): ax.scatter(*p, color=KNOWN, s=40, depthshade=False)
        else: ax.scatter(*p, color=GREEN, s=25, depthshade=False)
    ax.set_box_aspect((1, 1, 1)); ax.set_axis_off(); ax.view_init(*view)

show = sorted([i for i in A if not A[i]["new_match"]], key=lambda i: float(P[i]["T"]))
rows = (len(show) + 1) // 2
fig = plt.figure(figsize=(16, 4.2 * rows + 1))
for k, i in enumerate(show):
    p1, p2, T = float(P[i]["p1"]), float(P[i]["p2"]), float(P[i]["T"])
    t, Xs = traj(suvakov(p1, p2), T, int(T * 600))
    ax = fig.add_subplot(rows, 4, 2 * k + 1)
    for b in range(3): ax.plot(Xs[:, 2 * b], Xs[:, 2 * b + 1], color=BODY[b], lw=0.8)
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(f"T = {T:.6f}, closes to {P[i]['digits']} digits\np = ({p1:.9f}, {p2:.9f})", fontsize=9.5)
    ax = fig.add_subplot(rows, 4, 2 * k + 2, projection="3d"); sphere(ax)
    ax.plot(*shape(Xs).T, color="#0891b2", lw=0.9)
    ax.set_title(f"loop type, {len(A[i]['word'])} letters (re-read)", fontsize=9.5)
fig.suptitle(f"The {len(show)} loop types not in Li–Liao's 695-family list, from the polished orbits "
             "(each a real periodic orbit; real space left, shape sphere right)", y=1.0)
fig.savefig(os.path.join(OUT, "fig10-unlisted-orbits.png"), dpi=130, bbox_inches="tight")
print("wrote fig10-unlisted-orbits.png with", len(show), "orbits")

# ---- Figure 9 rings: same landscape, rings coloured by the re-read loop types
sys.path.insert(1, os.path.dirname(HERE))
import checks_lib as C
RES = json.load(open(os.path.join(os.path.dirname(HERE), "results.json")))
L = np.load(os.path.join(os.path.dirname(HERE), "landscape.npz")); E = RES["E"]
match = {i: bool(A[i]["new_match"]) for i in A}
fig, axs = plt.subplots(1, 2, figsize=(17, 7.6), layout="constrained")
for ax, key, title in ((axs[0], "dabs", "A. How close each start comes back (dark = closer)"), (axs[1], "tabs", "B. When it comes closest (time)")):
    F = L[key]
    im = ax.imshow(np.log10(F) if key == "dabs" else F, origin="lower", extent=(0, 0.8, 0, 0.8), cmap="magma" if key == "dabs" else "viridis",
                   vmin=-2 if key == "dabs" else 0, vmax=0.3 if key == "dabs" else 40)
    plt.colorbar(im, ax=ax, shrink=0.85, label="log₁₀ closest return" if key == "dabs" else "time of closest return")
    for idx, f in enumerate(E):
        if f["kind"] == "absolute" or abs(f["alpha"]) < 1e-6:
            m = match.get(idx, f["match"])
            ax.plot(f["p1"], f["p2"], "o", mfc="none", mec=KNOWN if m else "#00e5ff", mew=1.6, ms=7 if m else 11)
    for name in ("figure-eight", "butterfly I", "moth I", "yin-yang Ia", "goggles", "dragonfly", "bumblebee"):
        a, b, _ = C.SUV[name]; ax.annotate(name, (a, b), xytext=(5, 5), textcoords="offset points", color="white", fontsize=8.5)
    ax.set_xlabel("p₁ (sideways push)"); ax.set_ylabel("p₂ (upward push)"); ax.set_title(title)
axs[0].text(0.01, 0.79, "orange rings: orbits whose loop type is in Li–Liao's 2017 list\ncyan rings: periodic orbits whose loop type is not in it (re-read, polished)", color="white", fontsize=9, va="top")
fig.suptitle("The minimum search. Each pixel starts the three bodies in a line (−1, 0, +1) with pushes (p₁, p₂), (p₁, p₂), (−2p₁, −2p₂); 449,773 starts, each run to time 40.")
fig.savefig(os.path.join(OUT, "fig9-return-landscape.png"), dpi=130, bbox_inches="tight")
print("wrote fig9-return-landscape.png;", sum(1 for i in match if not match[i]), "cyan")
