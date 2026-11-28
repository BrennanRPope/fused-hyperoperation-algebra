import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig = plt.figure(figsize=(16, 5.2))
# a. phase helices around clocks' lines, far-seen
A = fig.add_subplot(1, 3, 1, projection="3d")
M = 1.0; T = np.linspace(0, 6, 600)
for x0, rr, col in [(0, np.inf, "#555"), (3, 5.0, "#2f6db5"), (6, 2.6, "#8e44ad"), (9, 2.05, "#c0392b")]:
    N = 1.0 if np.isinf(rr) else np.sqrt(1 - 2*M/rr)
    ph = 2*np.pi*1.2*N*T
    A.plot(x0 + 0.6*np.cos(ph), 0.6*np.sin(ph), T, color=col, lw=1.2)
    A.plot([x0, x0], [0, 0], [0, 6], color=col, lw=0.6, ls=":")
    A.text(x0 - 0.6, 0, 6.4, "far away" if np.isinf(rr) else f"r = {rr}M\nN = {N:.2f}", fontsize=7, color=col)
A.set_xticks([]); A.set_yticks([]); A.set_zlabel("far-away time")
A.set_title("a. Same clock (same E, winds at E/ħ), seen from far:\ncloser to the mass, the helix winds slower; at r = 2M it stops", fontsize=9)
A.view_init(elev=12, azim=-70)
# b. complex plane: (1+iv) and its square
B = fig.add_subplot(1, 3, 2)
th = np.linspace(0, np.pi/2, 100)
for vv, col, lab in [(0.45, "#2f6db5", "v = 0.45")]:
    w = 1 + 1j*vv; w2 = w*w
    B.annotate("", (w.real, w.imag), (0, 0), arrowprops=dict(arrowstyle="->", color=col, lw=1.5))
    B.annotate("", (w2.real, w2.imag), (0, 0), arrowprops=dict(arrowstyle="->", color="#c0392b", lw=1.8))
    B.text(w.real+0.03, w.imag-0.06, "1 + x, x = i·v\n(the rim, or the river)", fontsize=8, color=col)
    B.text(w2.real+0.03, w2.imag+0.02, "1 + x∘x = (1+iv)²", fontsize=8, color="#c0392b")
    B.plot([w2.real, w2.real], [0, w2.imag], ls="--", color="gray", lw=0.8)
    B.annotate("real part 1 − v²:\nthe slowdown (squared)", (w2.real/2, 0), xytext=(0.05, -0.35), fontsize=8,
               arrowprops=dict(arrowstyle="->", lw=0.7))
    B.annotate("imaginary ÷ real = 2v/(1−v²):\nthe helix's turns per loop", (w2.real, w2.imag/2), xytext=(0.95, 0.15), fontsize=8,
               arrowprops=dict(arrowstyle="->", lw=0.7))
B.axhline(0, color="gray", lw=0.5); B.axvline(0, color="gray", lw=0.5)
B.set_xlim(-0.1, 1.9); B.set_ylim(-0.5, 1.1); B.set_aspect("equal")
B.set_title("b. One composition, two answers: the law (k = 1) on x = i·v,\ncomposed with itself", fontsize=9)
# c. phase per loop = 2S
C = fig.add_subplot(1, 3, 3)
Ms = np.linspace(0.1, 2, 100)
C.plot(Ms, 8*np.pi*Ms**2, color="#2f6db5", label="no spin: M·β = 8πM² = 2S")
for s in [0.3, 0.6, 0.9]:
    Mk = np.linspace(0.3, 2, 6); a = s*Mk; rp = Mk + np.sqrt(Mk**2 - a**2)
    kap = (rp - Mk)/(rp**2 + a**2); Om = a/(rp**2 + a**2); beta = 2*np.pi/kap
    C.plot(Mk, beta*(Mk - 2*Om*Mk*a), "o", ms=4, label=f"spin a/M = {s}: β(M − 2ΩJ)")
    C.plot(Mk, 2*np.pi*(rp**2 + a**2), "x", color="k", ms=4)
C.set_xlabel("mass M (Planck units)"); C.set_ylabel("phase per imaginary-time loop (radians)")
C.legend(fontsize=7, frameon=False); C.text(0.15, 75, "× = 2S from the horizon area", fontsize=7)
C.set_title("c. Bigger mass: winds faster and its loop is longer;\nphase per loop = 2 × entropy (known: Smarr)", fontsize=9)
fig.suptitle("DRAFT · Phase winding at E/ħ and gravity's clock slowdown", fontsize=11)
fig.tight_layout(); fig.savefig("figures/c2g_winding_slowdown.png", dpi=150); print("saved")
