import numpy as np, json, glob, mpmath as mp, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 10})
def m_of_x(x):
    q = mp.exp(-mp.pi*x); return (mp.jtheta(2, 0, q)/mp.jtheta(3, 0, q))**4
def amp(y):
    x = mp.sqrt(y); return float(2*mp.asin(mp.sqrt(m_of_x(x)))*180/mp.pi)
# ---- Fig A: zeta as a weight over pendulum swings
d = json.load(open("swingdensity.json")); th = np.array(d["theta"]); w = np.array(d["w"]); u = np.array(d["u"])
fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
ax[0].plot(th, u, color="k"); ax[0].axhline(0, color="gray", lw=0.8); ax[0].axvline(90, color="tab:blue", ls="--")
ax[0].set_xlabel("swing amplitude θ₀ (degrees)"); ax[0].set_ylabel("Riemann's height u = ½ ln(K′/K)")
ax[0].set_title("A1. Riemann's height variable is the log of the\npendulum's period ratio (imaginary over real)")
ax[0].text(93, 0.35, "90° swing: u = 0\n(square torus, self-dual)", color="tab:blue", fontsize=9)
ax[0].text(8, 0.85, "rest side: K′ ≫ K, u > 0", fontsize=8); ax[0].text(100, -0.9, "separatrix side: K ≫ K′, u < 0", fontsize=8)
ax[0].set_ylim(-1.1, 1.1)
ax[1].fill_between(th, w, color="tab:red", alpha=0.35); ax[1].plot(th, w, color="tab:red")
ax[1].axvline(90, color="tab:blue", ls="--")
half = th[w > w.max()/2]
ax[1].annotate("", xy=(half.max(), w.max()/2), xytext=(half.min(), w.max()/2), arrowprops=dict(arrowstyle="<->"))
ax[1].text(half.min()+3, w.max()/2*1.07, f"above half its peak: {half.min():.0f}° to {half.max():.0f}°", fontsize=9)
ax[1].set_xlabel("swing amplitude θ₀ (degrees)"); ax[1].set_ylabel("Riemann's Φ, per degree of swing")
ax[1].set_title("A2. ζ's spin is a pendulum: Φ spreads over swings,\npeaked at 90°, mirror-symmetric (θ₀ ↔ 180° − θ₀)")
ax[2].axis("off")
ax[2].text(0.0, 0.95, "Dictionary (checked to 15 digits)", fontsize=11, weight="bold", transform=ax[2].transAxes)
lines = ["Jacobi theta  θ(x) = Σ e^(−πn²x)  =  √(2K/π)   when  x = K′/K",
         "Riemann's height  u = ½ ln x = ½ ln(K′/K)",
         "Functional equation s ↔ 1 − s  =  swap K and K′ (m ↔ 1 − m)",
         "Critical line  =  mirror symmetry about the 90° swing",
         "Field z (the RH variable) couples to u: Ξ(z) = ∫ Φ(u) cos(zu) du",
         "RH  ⇔  this pendulum 'spin' is Lee–Yang in the field z",
         "",
         "Dial-gas piece: separatrix (m = 1) = phase transition.",
         "Here: the 90° swing (m = ½) = ζ's centre."]
for i, L in enumerate(lines): ax[2].text(0.0, 0.84 - 0.09*i, L, fontsize=9, transform=ax[2].transAxes)
fig.tight_layout(); fig.savefig("figures/figA_zeta_as_pendulum.png", dpi=150); plt.close(fig)
# ---- Fig B: one swing = one lattice shape; zeros leave the line except at Euler-product swings
rows = []
for f in glob.glob("zeros_y*.json"):
    z = json.load(open(f)); rows.append((z["y"], z["offline"], z["total"]))
rows.sort()
fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
ys = np.array([r[0] for r in rows]); off = np.array([r[1] for r in rows]); amps = np.array([amp(y) for y in ys])
cn1 = {1: "ζ·β", 2: "ζ·L₋₈", 3: "ζ·L₋₃", 4: "ζ·β", 7: "ζ·L₋₇"}
for yy, o, a in zip(ys, off, amps):
    col = "tab:green" if yy in cn1 else ("tab:orange" if o == 0 else "tab:red")
    h = o if o else 0.25
    ax[0].bar(np.log(yy), h, width=0.11, color=col)
    lab = "" if yy == 51 else ("y=49, 51" if yy == 49 else f"y={yy:g}" + (f"\n{a:.1f}°" if yy in cn1 else ""))
    ax[0].text(np.log(yy), h+0.12, lab, ha="center", fontsize=7)
ax[0].set_xlabel("ln y, torus shape (K′/K)² = y; the swing amplitude falls from 90° at y = 1")
ax[0].set_ylabel("zeros off the line, height ≤ 30")
ax[0].set_title("B1. Turning the swing: zeros leave the line\nexcept at the Euler-product swings (green)")
ax[0].set_ylim(0, 8.2)
ax[0].text(0.0, 7.0, "green: Euler-product swings (class number one): none off the line\nred: off-line zeros found below height 30\norange: none below 30 (long tori; expected higher up)", fontsize=8)
# real-zero threshold
yy = np.logspace(0, 2.2, 200)
ax[1].semilogx(yy, [amp(v) for v in yy], color="k")
ystar = json.load(open("pend.json"))["ystar"]
ax[1].axvline(ystar, color="tab:red", ls="--"); ax[1].text(ystar*1.05, 60, f"ratio √y = {np.sqrt(ystar):.4f}:\na real zero enters (½, 1)\n(Bateman–Grosswald 7.0556)\nswing {amp(ystar):.3f}°", fontsize=8, color="tab:red")
for c in (1, 2, 3, 4, 7): ax[1].plot(c, amp(c), "o", color="tab:green")
ax[1].set_xlabel("torus shape y = (K′/K)²  (y ↔ 1/y gives the mirror swing)"); ax[1].set_ylabel("swing amplitude θ₀ (degrees)")
ax[1].set_title("B2. The swing dial and its far wall: very long tori\n(tiny or near-separatrix swings) get real off-line zeros")
fig.tight_layout(); fig.savefig("figures/figB_swing_dial_epstein.png", dpi=150); plt.close(fig)
print("ok", rows)
