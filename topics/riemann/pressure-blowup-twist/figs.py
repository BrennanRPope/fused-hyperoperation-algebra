import json, math
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
O = json.load(open("pbt.json"))
C = np.load("collapse.npz"); V = np.load("velocity.npz"); T = np.load("twist_traj.npz")

# ---------------- Fig A: blow-up ---------------------------------------------------------------------
fig, ax = plt.subplots(1, 2, figsize=(14, 5.4))
a = ax[0]
left, gaps = C["left"], C["gaps"]
a.loglog(left, gaps, "o", color="C3", ms=4, label=f"our flow, rewound: gap of a colliding pair (slope {O['A1_collapse']['exponent_gap_vs_time_left']:.3f})")
x = np.array([1e-7, 3e-2])
ref = gaps[len(gaps) // 2] / left[len(left) // 2]**0.5
for p, lab, c, ls in [(0.5, "Keller–Segel particles (overdamped log pull): ½  ✓ same law", "C3", "-"),
                      (0.4, "bubble collapse, Rayleigh 1917 (with inertia): 2/5", "C0", "--"),
                      (2 / 3, "gravitational free fall (with inertia): 2/3", "C2", ":")]:
    a.loglog(x, ref * x**p * (left[len(left) // 2]**(0.5 - p)), ls, color=c, lw=1, label=lab)
a.set_xlabel("heat time left before the collision"); a.set_ylabel("gap between the two zeros")
a.set_title("A1. A backward collision is a finite-time collapse with gap ∝ √(time left),\n"
            "gap² = 8 × time left exactly: the law of overdamped particles pulling as 1/r", fontsize=10)
a.legend(fontsize=7.5, loc="upper left")
a = ax[1]
r = O["A2_threshold"]
xs = np.array([q["S_d0sq"] for q in r]); ts = np.array([q["t"] for q in r]); own = np.array([q["own_pair"] for q in r])
tf = np.array([q["t_formula_frozen"] if q["t_formula_frozen"] is not None else np.nan for q in r])
a.semilogy(xs[own], ts[own], "o", color="C3", label="the chosen pair collapses")
a.semilogy(xs[~own], ts[~own], "s", color="C1", label="a neighbouring pair collapses instead")
a.semilogy(xs, tf, "k-", lw=1, label="frozen-crowd formula (collapses only if S·d₀² < 2)")
a.axvline(2, color="0.5", lw=0.8); a.text(2.05, 1.6, "threshold S·d₀² = 2\n(the 'critical mass')", fontsize=8)
a.set_xlabel("crowd squeeze × gap²  (S·d₀²) of the chosen pair"); a.set_ylabel("rewind time to the first collision")
a.set_title("A2. Like Keller–Segel, each pair has a threshold. Unlike Keller–Segel,\n"
            "staying under it never saves the system: some other pair always collapses", fontsize=10)
a.legend(fontsize=8, loc="lower right")
plt.tight_layout(); plt.savefig("figures/figA_blowup.png", dpi=130); plt.close()

# ---------------- Fig B: pressure --------------------------------------------------------------------
fig, ax = plt.subplots(1, 3, figsize=(18, 5.3))
a = ax[0]
from numpy.polynomial.hermite import hermroots
h = np.sort(hermroots([0] * 12 + [1]))
from scipy.integrate import solve_ivp
def f(t, z):
    D = z[:, None] - z[None, :]; np.fill_diagonal(D, np.inf); return 2 * (1 / D).sum(1) - 1.0 * z
z0 = np.sort(np.random.default_rng(3).normal(0, 2, 12))
zeq = np.sort(solve_ivp(f, [0, 40], z0, rtol=1e-12, atol=1e-14, method="DOP853").y[:, -1])
a.plot(z0, np.zeros(12) - 0.25, ".", color="0.5", ms=7, label="12 random starting zeros")
a.plot(zeq, np.zeros(12), "o", color="C2", ms=9, label="the same zeros after the pressured flow (a = 1) settles")
a.plot(h * math.sqrt(2), np.zeros(12) + 0.15, "x", color="k", ms=8, label="zeros of the Hermite polynomial H₁₂ (×√2)")
a.set_ylim(-0.4, 0.6); a.set_yticks([])
a.set_xlabel("position")
a.set_title(f"B1. With an inward pressure −a·z the zeros DO stop: at the Hermite\nzeros (Stieltjes 1885). Match to {O['B1_hermite']['max_err']:.0e}", fontsize=10)
a.legend(fontsize=8, loc="upper left")
a = ax[1]
s = np.linspace(-3, 1.2, 400)
for aa, c in [(1.0, "C0"), (2.0, "C2")]:
    a.plot(s, (np.exp(2 * aa * s) - 1) / (2 * aa), color=c, label=f"pressure a = {aa:g}: real heat time t(s)")
    a.axhline(-1 / (2 * aa), color=c, ls=":", lw=0.8)
a.plot(s, s, "k--", lw=0.8, label="no pressure: t = s")
Tc = O["A1_collapse"]["rewind_time_to_collision"] * 20
a.set_ylim(-1.2, 2); a.axhline(0, color="k", lw=0.4)
a.set_xlabel("time on the pressured clock s"); a.set_ylabel("real heat time t reached")
a.set_title("B2. Pressure is a zoom plus a slower clock (exact, error {:.0e}).\nRun backward forever and t only reaches −1/(2a): it hides collisions, never removes them".format(O["B2_zoom_clock"]["max_err"]), fontsize=10)
a.text(-2.9, -0.62, "t = −1/(2a): the backward wall", fontsize=8, color="C0")
a.legend(fontsize=8, loc="upper left")
a = ax[2]
zz, res = V["z"], V["res"]
m = slice(0, 120)
a.bar(zz[m] / 2, -res[m], width=0.6, color=np.where(-res[m] > 0, "C2", "C3"))
a.axhline(0, color="k", lw=0.5)
a.set_xlabel("height γ (zeros 1000–1119 shown)"); a.set_ylabel("pressure needed to hold each zero still at t = 0")
bv = O["B4_velocity"]
a.set_title("B3. The pressure that would freeze ζ's zeros at t = 0 is not one number:\n"
            f"it flips sign at {bv['frac_sign_changes']*100:.0f}% of neighbours (rms {bv['rms_residual_velocity']:.2f}, 1.7× the smooth drift)", fontsize=10)
plt.tight_layout(); plt.savefig("figures/figB_pressure.png", dpi=130); plt.close()

# ---------------- Fig C: twist -----------------------------------------------------------------------
fig, ax = plt.subplots(1, 2, figsize=(14, 5.3))
a = ax[0]
for k, lab, c in [("0.0000", "θ = 0 (plain heat)", "k"), ("0.3927", "θ = π/8", "C0"), ("0.7854", "θ = π/4", "C2"),
                  ("1.5708", "θ = π/2 (Schrödinger-type)", "C3")]:
    d = O["C_twist"][k]
    a.plot(d["tau"], d["max_im"], "o-" if k != "0.0000" else "-", ms=3, color=c, label=lab)
a.set_xlabel("twisted heat time |t| (along t = |t|·e^{iθ})"); a.set_ylabel("furthest any zero is off the line")
a.set_title("C1. Twist the heat time by any angle and the zeros leave the line\nat once, at speed ∝ sin θ (300 of our zeros, γ ≈ 2700–2850)", fontsize=10)
a.legend(fontsize=8)
a = ax[1]
Z = T["Z"]
for j in range(0, Z.shape[1], 1):
    a.plot(Z[:, j].real / 2, Z[:, j].imag / 2, color="C3", lw=0.6, alpha=0.7)
a.plot(Z[0].real / 2, Z[0].imag / 2, "k.", ms=3, label="start: on the line")
a.plot(Z[-1].real / 2, Z[-1].imag / 2, "o", color="C3", ms=2.5, label="after |t| = 0.05 at θ = π/2")
a.set_xlim(0, 60); a.set_xlabel("position along the line (γ, offset)"); a.set_ylabel("off the line (γ units)")
a.set_title("C2. Under the Schrödinger-type flow each zero moves sideways\nby its own heat-flow speed, turned 90°. Untwisting is just running it back", fontsize=10)
a.legend(fontsize=8, loc="lower left")
plt.tight_layout(); plt.savefig("figures/figC_twist.png", dpi=130); plt.close()
print("ok")
