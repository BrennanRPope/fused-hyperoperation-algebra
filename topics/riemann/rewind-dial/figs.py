"""Rewind-dial figures. DRAFT."""
import numpy as np, json, sys, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, "../knauf-field"); from knauf import chain
from rotate import setup, F, start_zeros
from oddpath import setup_odd
plt.rcParams.update({"font.size": 10})
mv = json.load(open("movie_n12.json")); od = {**json.load(open("oddpath_8_10.json")), **json.load(open("oddpath_12.json"))}

# ---------- Fig A: the rewind dial ----------
fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
a = ax[0]; t = np.linspace(0, np.pi/2, 100)
a.plot(np.cos(t), np.sin(t), color="tab:purple", lw=3)
a.plot([0, 1.0], [0, 0], color="0.5"); a.plot([0, 0], [0, 1.0], color="0.5")
a.plot(1, 0, "o", color="k", ms=8); a.plot(0, 1, "o", color="k", ms=8)
a.text(1.05, 0.02, "0°: magnetic field\n(fixed chain passes Lee–Yang:\nevery zero on the line)", fontsize=8.5, va="bottom")
a.text(0.08, 1.04, "90°: height field = p's Wick turn (the original temperature view;\nzeros on the line here ⇔ RH for the chain)", fontsize=8.5, va="bottom")
t2 = np.linspace(0.25, 1.3, 50); a.plot(0.62*np.cos(t2), 0.62*np.sin(t2), color="tab:purple", lw=1.2)
a.annotate("", (0.62*np.cos(1.32), 0.62*np.sin(1.32)), (0.62*np.cos(1.25), 0.62*np.sin(1.25)), arrowprops=dict(arrowstyle="->", color="tab:purple"))
a.text(0.05, 0.12, "rewind angle θ", color="tab:purple", fontsize=10)
a.text(0.5, -0.08, "push on the spins (magnetisation)", fontsize=8.5, ha="center", va="top")
a.text(-0.06, 0.5, "push on the energy (height ln φ)", fontsize=8.5, ha="right", va="center", rotation=90)
a.set_xlim(-0.2, 2.1); a.set_ylim(-0.25, 1.35); a.axis("off")
a.set_title("A1. The rewind: turn the field from the spins to the energy")
b = ax[1]; ns = [8, 10, 12, 14, 16]; odd, even = [], []
for n in ns:
    p = chain(n).astype(float); H = np.log(p); Ho = 0.5*(H-H[::-1]); He = 0.5*(H+H[::-1])
    odd.append(np.var(Ho)/(np.var(Ho)+np.var(He))); even.append(1-odd[-1])
b.bar(ns, odd, color="tab:green", label="odd half (the chain's built-in field):\nrewinds with every zero on the line")
b.bar(ns, even, bottom=odd, color="tab:red", label="even half (links that pair spins):\nrewinding it pushes zeros off the line")
for n, o in zip(ns, odd): b.text(n, o/2, "%.0f%%" % (100*o), ha="center", color="w", fontsize=9)
b.set_xlabel("chain length n"); b.set_ylabel("share of the energy's spread"); b.set_ylim(0, 1.45)
b.legend(loc="upper center", fontsize=8.5); b.set_title("A2. Flipping every spin splits the energy in two")
fig.tight_layout(); fig.savefig("figures/figA_rewind_dial.png", dpi=150); plt.close(fig)

# ---------- Fig B: movie of zeros along the full rewind ----------
fig, ax = plt.subplots(1, 2, figsize=(12, 4.8))
degs = sorted(int(d) for d in mv["frames"])
for d in degs:
    z = np.array(mv["frames"][str(d)])
    if z.size: sc = ax[0].scatter(np.full(len(z), d), z[:, 0], c=z[:, 1], cmap="viridis", vmin=0, vmax=36, s=12)
ax[0].axhline(0, color="k", lw=0.8); cb = fig.colorbar(sc, ax=ax[0]); cb.set_label("height of the zero (Im z)")
ax[0].set_xlabel("rewind angle θ (degrees)"); ax[0].set_ylabel("distance from the line (Re z)")
ax[0].set_title("B1. Full rewind, n = 12: zeros leave the line at once")
for d, c in ((0, "k"), (5, "tab:blue"), (45, "tab:orange"), (90, "tab:red")):
    z = np.array(mv["frames"][str(d)]); ax[1].scatter(z[:, 0], z[:, 1], s=14, color=c, label="θ = %d°" % d)
ax[1].axvline(0, color="k", lw=0.8); ax[1].set_xlabel("Re z   (at 90°: Re z = 3/2 − Re β)"); ax[1].set_ylabel("Im z")
ax[1].set_title("B2. Snapshots: on the line at 0°, spread out by 5°"); ax[1].legend(loc="upper right", fontsize=8.5)
fig.tight_layout(); fig.savefig("figures/figB_full_rewind_movie.png", dpi=150); plt.close(fig)

# ---------- Fig C: first-order speed off the line ----------
fig, ax = plt.subplots(figsize=(8, 4.6))
for n, c in ((8, "tab:blue"), (10, "tab:green"), (12, "tab:orange"), (14, "tab:red")):
    S = setup(n); z0 = start_zeros(S, 36.0); G = (S["w"]*S["H"]*np.exp(np.multiply.outer(z0, S["a"]*S["M"]))).sum(-1)
    _, Fz = F(S, z0, 0.0); zp = -(z0*G)/Fz; ax.scatter(z0.imag, np.abs(zp.real), s=14, color=c, label="n = %d" % n)
So = setup_odd(10); z0 = start_zeros(So, 36.0); G = (So["w"]*So["H"]*np.exp(np.multiply.outer(z0, So["a"]*So["M"]))).sum(-1)
_, Fz = F(So, z0, 0.0); zp = -(z0*G)/Fz; ax.scatter(z0.imag, np.abs(zp.real) + 1e-3, marker="x", color="k", label="odd half only, n = 10 (exactly 0, drawn at the floor)")
ax.set_yscale("log"); ax.set_ylim(5e-4, 100); ax.set_xlabel("height of the zero at θ = 0 (Im z)")
ax.set_ylabel("speed off the line as θ leaves 0   |Re dz/dθ|")
ax.set_title("C. The full rewind pushes every zero off at first order: faster higher up (in waves),\nand the same for every chain length")
ax.legend(fontsize=8.5, loc="lower right")
fig.tight_layout(); fig.savefig("figures/figC_first_order_speed.png", dpi=150); plt.close(fig)

# ---------- Fig D: odd half vs full ----------
fig, ax = plt.subplots(1, 2, figsize=(12, 4.6), sharey=True)
for d in sorted(int(x) for x in od["12"]):
    z = np.array(od["12"][str(d)])
    if z.size: ax[0].scatter(np.full(len(z), d), z[:, 0], c=z[:, 1], cmap="viridis", vmin=0, vmax=36, s=12)
ax[0].axhline(0, color="k", lw=0.8); ax[0].set_ylim(-1.05, 1.05)
ax[0].text(45, 0.45, "15 to 22 zeros at each angle,\nall within 1e-8 of the line", ha="center", fontsize=9.5)
ax[0].set_title("D1. Rewinding the odd half only (n = 12):\nevery zero stays on the line, 0° to 90°")
ax[0].set_xlabel("rewind angle θ (degrees)"); ax[0].set_ylabel("distance from the line (Re z)")
for d in degs:
    z = np.array(mv["frames"][str(d)])
    if z.size: ax[1].scatter(np.full(len(z), d), z[:, 0], c=z[:, 1], cmap="viridis", vmin=0, vmax=36, s=12)
ax[1].axhline(0, color="k", lw=0.8); ax[1].set_xlabel("rewind angle θ (degrees)")
ax[1].set_title("D2. Rewinding the whole energy (n = 12):\nthe even half throws zeros off")
fig.tight_layout(); fig.savefig("figures/figD_odd_vs_full.png", dpi=150); plt.close(fig)
print("figs done")
