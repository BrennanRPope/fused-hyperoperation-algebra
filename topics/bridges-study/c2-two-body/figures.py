# DRAFT. Figure for the C2 follow-up. Colour key: green = Euclidean Kerr (known GR), magenta = ours/the law's reading,
# orange = the model under test (two half-masses), grey = reference.
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
GR, MA, OR, GY = "#2e8b57", "#c2187a", "#dd8a1e", "#8a8a8a"
fig = plt.figure(figsize=(15, 4.8))
ax = fig.add_subplot(1, 3, 1); al = 0.6; l = np.arange(6); w = 0.27
ax.bar(l-w, (-al)**l, w, color=GR, label="Euclidean Kerr: M(−α)ˡ")
ax.bar(l, 0.5*(al**l+(-al)**l), w, color=OR, label="two half-masses at ±α (the river's two-body rule)")
ax.bar(l+w, (-al)**l, w, color=MA, alpha=0.55, label="one mass moved to −α")
ax.axhline(0, color="black", lw=0.7); ax.set_xticks(l); ax.set_xlabel("multipole order l (0 = mass, 1 = dipole, …)")
ax.set_ylabel("moment (M = 1, α = 0.6)"); ax.legend(fontsize=7.6, loc="upper right")
ax.set_title("c2d-a. The two half-masses fail at the dipole:\nKerr's moments are one mass moved by α (G3, G4)", fontsize=9.5)
ax = fig.add_subplot(1, 3, 2); x = np.linspace(0, 4, 300)
ax.plot(x, x/np.sqrt(1+x**2), color=MA, lw=2.4, label="turn per period / 2π = α/√(M²+α²) = tanh u")
ax.axhline(1, color=GY, ls=":"); ax.set_xlabel("α / M  (= sinh u)"); ax.set_ylabel("fraction of a full turn")
ax.set_title("c2d-b. Each trip round imaginary time turns the\nframe by 2π·tanh u, never a full turn (G5)", fontsize=9.5); ax.legend(fontsize=8, loc="lower right")
ax = fig.add_subplot(1, 3, 3, projection="3d"); frac = 0.6/np.sqrt(1+0.36)
s = np.linspace(0, 3, 600); ang = 2*np.pi*s
ax.plot(np.cos(ang), np.sin(ang), s*frac, color=MA, lw=2.4)
for k in range(4): ax.plot(np.cos(np.linspace(0, 2*np.pi, 100)), np.sin(np.linspace(0, 2*np.pi, 100)), k*frac, color=GY, lw=0.6)
ax.set_xlabel("imaginary time τ (round)"); ax.set_zlabel("frame turn φ / 2π"); ax.set_xticks([]); ax.set_yticks([])
ax.set_title("c2d-c. Spin in imaginary time is a helix:\nτ goes round once, φ climbs by tanh u", fontsize=9.5)
fig.suptitle("c2d. The river's two-body rule doesn't rebuild Kerr; spin lives in a helical twist of imaginary time", fontsize=11.5, weight="bold")
fig.tight_layout(); fig.savefig("figures/c2d_two_body_test.png", dpi=150, bbox_inches="tight"); print("ok")
