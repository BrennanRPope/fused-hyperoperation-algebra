import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 3, figsize=(15, 4.8))
# a. dial plane
A = ax[0]
kk = np.linspace(-2.5, 2.5, 300)
A.fill_between(kk, kk**2/4, 3, color="#f3e3cf", alpha=0.7)
A.fill_between(kk, -3, kk**2/4, color="#d7e5f5", alpha=0.7)
A.plot(kk, kk**2/4, color="gray", lw=0.8)
A.axhline(0, color="#333", lw=2)
A.text(0.9, -0.55, "d = 0: the law as it is\n(a+b+kab)", fontsize=8)
A.plot([0], [1], "o", color="#c0392b"); A.plot([0], [-1], "o", color="#2f6db5")
A.annotate("real-time Kerr helix\n(k=0, d=+1: rotation)", (0, 1), xytext=(-2.4, 1.9), fontsize=8, arrowprops=dict(arrowstyle="->", lw=0.7))
A.annotate("imaginary-time Kerr helix\n(k=0, d=-1: boost)", (0, -1), xytext=(-2.4, -2.3), fontsize=8, arrowprops=dict(arrowstyle="->", lw=0.7))
A.annotate("", (0.08, -0.9), (0.08, 0.9), arrowprops=dict(arrowstyle="<->", color="purple"))
A.text(0.18, -0.15, "× i flips d", color="purple", fontsize=8)
A.text(1.3, 2.4, "rotation-like", fontsize=8, color="#8a5a1f"); A.text(1.3, -2.6, "boost-like", fontsize=8, color="#2f6db5")
A.set_xlim(-2.5, 2.5); A.set_ylim(-3, 3); A.set_xlabel("k (first dial)"); A.set_ylabel("d (second dial, already in the repo)")
A.set_title("a. The two-dial law: a∘b = (a+b+kab)/(1 − d·ab)\n(Coleman–McGuinness 1991); d = 0 is the law as it is", fontsize=10)
# b. helix = rim composed with itself
B = ax[1]
vv = np.linspace(0, 0.95, 200)
B.plot(vv, 2*vv/(1-vv**2), color="#c0392b", label="real time: v∘v at d=+1  (= Ω_H/κ)")
B.plot(vv, 2*vv/(1+vv**2), color="#2f6db5", label="imaginary time: t∘t at d=−1  (= tanh u)")
M = 1.0
for s in np.linspace(0.1, 0.95, 7):
    a = s*M; rp = M+np.sqrt(M**2-a**2); B.plot(a/rp, a/np.sqrt(M**2-a**2), "o", color="#c0392b", ms=4)
    al = 3*s*M; rpE = M+np.sqrt(M**2+al**2); B.plot(al/rpE, al/np.sqrt(M**2+al**2), "s", color="#2f6db5", ms=4)
B.set_ylim(0, 4); B.set_xlabel("horizon rim speed (v = a/r₊, or t = α/r₊)"); B.set_ylabel("helix: turns per loop ÷ 2π")
B.text(0.05, 3.3, "dots = actual Kerr horizons", fontsize=8, color="#555")
B.legend(fontsize=8, frameon=False, loc="center left")
B.set_title("b. The helix is the horizon rim composed with itself", fontsize=10)
# c. the multiplicative coordinate read off the moved mass
C = ax[2]
th = np.linspace(0, 2*np.pi, 300); C.plot(np.cos(th), np.sin(th), color="#c0392b", lw=0.8)
chi = 0.8; C.plot([0, np.cos(chi)], [0, np.sin(chi)], color="#c0392b"); C.plot(np.cos(chi), np.sin(chi), "o", color="#c0392b")
C.text(np.cos(chi)+0.05, np.sin(chi)+0.05, "$(r_+ - M + ia)/M = e^{i\\chi}$\nreal time: on the circle", fontsize=8, color="#c0392b")
uu = 0.9; C.plot([0, np.exp(uu)], [0, 0], color="#2f6db5"); C.plot(np.exp(uu), 0, "s", color="#2f6db5")
C.text(np.exp(uu)-1.5, -0.5, "$(r_+ + \\alpha - M)/M = e^{u}$\nimaginary time: on the line", fontsize=8, color="#2f6db5")
C.axhline(0, color="gray", lw=0.4); C.axvline(0, color="gray", lw=0.4)
C.set_aspect("equal"); C.set_xlim(-1.3, 3.0); C.set_ylim(-1.3, 1.5)
C.set_title("c. The law's multiplier at the horizon pole:\nhelix = tan χ (circle) or tanh u (line)", fontsize=10)
fig.suptitle("DRAFT · The Kerr helix on the law's d dial", fontsize=11)
fig.tight_layout(); fig.savefig("figures/c2f_helix_from_law.png", dpi=150); print("saved")
