# DRAFT. Figures for candidates C1-C3. Colour key: green = proved/known physics, magenta = the law / ours,
# orange = a model under test, grey = reference.
import json, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
GR, MA, OR, GY = "#2e8b57", "#c2187a", "#dd8a1e", "#8a8a8a"
R = json.load(open("results.json"))

# c1
fig, axs = plt.subplots(1, 2, figsize=(12.5, 4.6))
ax = axs[0]; ratio = np.logspace(-2, 2, 200)
for A1, col in [(1, OR), (10, "#b5651d"), (100, "#7a3e00")]:
    A2 = A1*ratio; S = lambda A: A**1.5
    ax.loglog(ratio, (S(A1+A2)-S(A1)-S(A2))/(S(A1)*S(A2)), color=col, lw=2, label=f"A^1.5 (Tsallis-Cirto), first patch area {A1}")
ax.axhline(1e-3, color="white", lw=0)
ax.set_xlabel("second patch area / first patch area"); ax.set_ylabel("k that the split would need")
ax.set_title("c1a. A power of area is not the law at any fixed k:\nthe 'k' changes with every split and size (check D1)", fontsize=10)
ax.text(0.012, 2e-4, "area law (δ = 1): k = 0 for every split", color=GR, fontsize=9); ax.legend(fontsize=8)
ax = axs[1]; lam = np.linspace(-0.8, 0.8, 400)
ax.axvspan(0, 0.8, color="#fde7f1"); ax.axvspan(-0.8, 0, color="#eef7ee")
ax.plot(lam, -lam/(1-lam), color=MA, lw=2.6)
ax.axvline(0, color="black", lw=0.8); ax.axhline(0, color="black", lw=0.8)
ax.text(0.05, -2.6, "λ > 0: negative tension\n(excess angle, negative energy)", fontsize=8.4)
ax.text(-0.78, 0.25, "λ < 0: positive tension\n(cone, deficit angle)", fontsize=8.4)
ax.text(0.03, -0.75, "λ = 0: area law,\nno brane; the only\npoint linear mixing allows", color="#5a0b38", fontsize=8.4)
ax.set_xlabel("dial λ put on the black-hole area (k = λ, Rényi n = 1 − λ)"); ax.set_ylabel("4G × brane tension")
ax.set_title("c1b. Through Dong's dictionary, any λ ≠ 0 on area needs a brane (check D5)", fontsize=10); ax.set_ylim(-4.2, 0.9)
fig.suptitle("c1. Non-additive horizon entropies sorted onto the law's dial", fontsize=11.5, weight="bold")
fig.tight_layout(); fig.savefig("figures/c1_horizon_entropies.png", dpi=150, bbox_inches="tight"); plt.close(fig)

# c2
fig, axs = plt.subplots(1, 3, figsize=(15, 4.7))
x = np.linspace(-3, 3, 401); X, Z = np.meshgrid(x, x); rho = np.abs(X)
a = 1.0; b = rho**2+Z**2-a**2; r = np.sqrt((b+np.sqrt(b*b+4*a*a*Z*Z))/2); HL = r**3/(r**4+a*a*Z*Z)
al = 1.0; HE = 0.5*(1/np.sqrt(rho**2+(Z+al)**2+1e-9)+1/np.sqrt(rho**2+(Z-al)**2+1e-9))
lv = np.linspace(0.15, 2.5, 18)
ax = axs[0]; ax.contour(X, Z, np.clip(HL, 0, 3), levels=lv, cmap="viridis"); ax.plot([-a, a], [0, 0], color=MA, lw=4)
ax.set_title("c2a. Lorentzian Kerr: the scalar is the field of a\nsource moved by i·a (a ring of radius a, magenta)", fontsize=9.5); ax.set_aspect("equal"); ax.set_xlabel("x"); ax.set_ylabel("z (spin axis)")
ax = axs[1]; ax.contour(X, Z, np.clip(HE, 0, 3), levels=lv, cmap="viridis"); ax.plot([0, 0], [-al, al], "o", color=MA, ms=10)
ax.set_title("c2b. Turn time imaginary (a = iα): the two i's give −1,\nand the scalar is two real half-masses at z = ±α (check E2)", fontsize=9.5); ax.set_aspect("equal"); ax.set_xlabel("x")
ax = axs[2]; aa = np.linspace(0, 0.999, 300); rp = 1+np.sqrt(1-aa**2); rm = 1-np.sqrt(1-aa**2)
ax.plot(aa, (rp-rm)/(4*np.pi*(rp**2+aa**2)), color=GR, lw=2.4, label="T from smoothness (κ/2π), check E3")
ax2 = ax.twinx(); ax2.plot(aa, np.pi*(rp**2+aa**2), color=MA, lw=2, ls="--", label="S = A/4 from the action, check E4")
ax.set_xlabel("spin a / M"); ax.set_ylabel("temperature (M = 1)", color=GR); ax2.set_ylabel("entropy", color=MA)
ax.set_title("c2c. Continued back to real time: Kerr's\ntemperature and entropy, from the shifted source", fontsize=9.5)
h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels(); ax.legend(h1+h2, l1+l2, fontsize=8, loc="lower left")
fig.suptitle("c2. The river's shift by i, in imaginary time, gives Kerr's thermodynamics", fontsize=11.5, weight="bold")
fig.tight_layout(); fig.savefig("figures/c2_euclidean_kerr.png", dpi=150, bbox_inches="tight"); plt.close(fig)

# c3
F = R["F3_track"]; N = F["N"]; st = np.array(F["steps"])
fig, axs = plt.subplots(1, 2, figsize=(12.5, 4.5))
ax = axs[0]
ax.plot(st, F["hawking"], color=GY, lw=2.4, label="Hawking's saddle only: no extra term, ever")
ax.plot(st, F["full"], color=MA, lw=2.4, marker="o", ms=4, label="with replica wormholes")
ax.axvline(N/2, color=GR, ls=":", lw=1.6); ax.text(N/2+0.3, -0.85, "Page time", color=GR, fontsize=9)
ax.set_xlabel(f"qubits radiated (of {N})"); ax.set_ylabel("extra term on the law's k  (q = 2)")
ax.set_title("c3a. Early vs late radiation: the law's extra term appears\nonly after the Page time, and only from wormholes (check F3)", fontsize=9.5); ax.legend(fontsize=8)
ax = axs[1]; rr = np.arange(0, N+1)
pur = lambda n: (2.0**n+2.0**(N-n))/(2.0**N+1)
ax.plot(rr, -np.log2(1/2.0**rr), color=GY, lw=2, label="Hawking: radiation keeps getting more mixed")
ax.plot(rr, -np.log2([pur(n) for n in rr]), color=MA, lw=2.4, label="with wormholes: the Page curve")
ax.axvline(N/2, color=GR, ls=":", lw=1.6)
ax.set_xlabel(f"qubits radiated (of {N})"); ax.set_ylabel("order-2 entropy of the radiation (bits)")
ax.set_title("c3b. The same two saddles give the Page curve\n(Lubkin 1978 average; PSSY 2019 replica wormholes)", fontsize=9.5); ax.legend(fontsize=8)
fig.suptitle("c3. Replica wormholes are where the law's extra term comes from", fontsize=11.5, weight="bold")
fig.tight_layout(); fig.savefig("figures/c3_page_cross_term.png", dpi=150, bbox_inches="tight"); plt.close(fig)

# c1c: the roundness dial p, where the horizon entropies land
fig, ax = plt.subplots(figsize=(5.6, 5.4)); t = np.linspace(0, 2*np.pi, 800)
def ball(p): c, s_ = np.cos(t), np.sin(t); return np.sign(c)*np.abs(c)**(2/p), np.sign(s_)*np.abs(s_)**(2/p)
for p, col, lab in [(2/3, OR, "p = 2/3: Tsallis-Cirto (δ = 3/2) and Barrow at Δ = 1"), (0.8, "#e9b26b", "p = 0.8: Barrow at Δ = 0.5"),
                    (1, MA, "p = 1: area law (Bekenstein-Hawking), plain addition"), (2, GY, "p = 2: the ball (reference)")]:
    x_, y_ = ball(p); ax.plot(x_, y_, color=col, lw=2.4 if p == 1 else 1.8, label=lab)
ax.set_aspect("equal"); ax.axhline(0, color="black", lw=0.5); ax.axvline(0, color="black", lw=0.5)
ax.set_title("c1c. Power-of-area entropies sit on the existing\nroundness dial p = 1/δ, on its non-convex side (check D7)", fontsize=9.5)
ax.legend(fontsize=7.6, loc="lower center", bbox_to_anchor=(0.5, -0.36)); fig.tight_layout()
fig.savefig("figures/c1c_roundness_dial.png", dpi=150, bbox_inches="tight"); plt.close(fig)
print("ok")
