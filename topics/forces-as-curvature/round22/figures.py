"""DRAFT picture for round 22. Run after checks.py. Blue = light, orange = matter, grey = reference."""
import json, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OR, BL, GR = "#e07b00", "#1f6fc2", "#777777"
R = json.load(open("results.json"))
gd = lambda t: 2 * np.arctan(np.tanh(np.asarray(t) / 2))
plt.rcParams.update({"font.size": 10.5, "figure.dpi": 130})
fig, ax = plt.subplots(2, 2, figsize=(16, 13)); ax = ax.ravel()

# A. conformal diagram: angle on the slice (chi) across, conformal time (eta = gd(tau)) up
a = ax[0]
a.fill_between([np.pi / 2, np.pi], 0, np.pi / 2, color=GR, alpha=0.12)
a.text(3 * np.pi / 4, 0.25, "never reached by light\nsent from the pole after the bounce\n(beyond the horizon)", ha="center", color=GR, fontsize=9)
for tau in (0, 0.5, 1, 2, 3):
    a.axhline(gd(tau), color=GR, lw=0.6, ls=":")
    a.text(np.pi + 0.03, gd(tau), f"τ = {tau}, slice radius a = cosh τ = {math.cosh(tau):.2f}", va="center", fontsize=8, color=GR)
a.plot([0, np.pi / 2], [0, np.pi / 2], color=BL, lw=2.6, label="the law's light ray from the bounce:\nχ = gd(τ) → π/2 (M3)")
for i, r in enumerate(R["horizon"]["later"]):
    e0 = gd(r["tau0"]); a.plot([0, np.pi / 2 - e0], [e0, np.pi / 2], color=BL, lw=1.4, ls="--", label="light sent later, τ0 = 0.5, 1, 2:\nreaches π/2 − gd(τ0)" if i == 0 else None)
for i, r in enumerate(R["matter"]):
    b_ = r["beta"]; tau = np.linspace(0, 12, 600)
    chi = np.arctan2(np.sinh(b_) * np.tanh(tau), 1.0)            # comoving angle of the particle (from its placed line)
    u = np.linspace(0, 0.9999999, 4000); D = 1 - u**2
    W, T, X = (1 + u**2) / D, 2 * u * np.cosh(b_) / D, 2 * u * np.sinh(b_) / D
    a.plot(np.arctan2(X, W), gd(np.arcsinh(T)), color=OR, lw=1.8, label="free particles leaving the pole,\nrapidity β = 0.5, 1, 2: slow to rest\nat χ = gd(β) (M5)" if i == 0 else None)
    a.text(r["gd_beta"] - 0.03, np.pi / 2 - 0.12, f"β = {b_}", color=OR, fontsize=8.5, ha="right")
a.set_xlim(0, np.pi); a.set_ylim(0, np.pi / 2 + 0.05)
a.set_xticks([0, np.pi / 4, np.pi / 2, 3 * np.pi / 4, np.pi]); a.set_xticklabels(["0 (pole)", "π/4", "π/2", "3π/4", "π (far pole)"])
a.set_yticks([0, np.pi / 4, np.pi / 2]); a.set_yticklabels(["0 (bounce)", "π/4", "π/2 (τ → ∞)"])
a.set_xlabel("angle χ on the slice's 3-sphere"); a.set_ylabel("conformal time η = gd(τ)")
a.legend(fontsize=8.5, loc="center right", bbox_to_anchor=(1.0, 0.62)); a.set_title("A. Light and matter in the law's de Sitter (future half, conformal picture)")

# B. a photon mode on the growing slices: the wave stretches and freezes, its energy falls as 1/a, the photon count stays
b = ax[1]
t = np.linspace(0, 5, 2000); k = 12; a_ = np.cosh(t)
b.plot(t, 2.3 + 0.8 * np.cos(k * gd(t)), color=BL, lw=1.3, label=f"mode amplitude cos(k η), k = {k}: its period stretches with a, then it freezes")
b.plot(t, 1 / a_, color=BL, lw=2.4, label="energy / starting energy = 1/a (the redshift, M2)")
b.plot(t, np.ones_like(t), color=GR, lw=1.5, ls="--", label="energy / local frequency (photon count) / start = 1")
rk = [r for r in R["redshift_modes"] if r["k"] == 12][0]
b.text(5.0, 0.2, f"integrated to t = 6: E(6)/E(0) = {rk['E_end_over_E0']:.5f},\n1/a = {1 / rk['a_end']:.5f}; E·a drift {rk['Ea_drift']:.0e}", ha="right", fontsize=8.5, color=GR)
b.set_xlabel("time τ since the bounce"); b.set_yticks([0, 0.5, 1]); b.set_ylim(0, 3.9)
b.legend(fontsize=8.5, loc="upper right"); b.set_title("B. A photon mode on the law's growing 3-spheres")

# C. the same light on the law's mesh
c = ax[2]
L = R["mesh_light"]; tt = np.array(L["proj_t"]); pr = np.array(L["proj"]); E2 = np.array(L["E2"])
ts = np.linspace(0, tt[-1], 400)
c.plot(ts, np.cos(math.sqrt(L["law_field_freq2"]) * gd(ts)), color=GR, lw=1.2, label=f"cos(√λ gd(τ)), λ = {L['law_field_freq2']:.3f} (the field's frequency² on the mesh)")
c.plot(tt, pr, "o", color=BL, ms=4, label=f"the law's own field (round 18, curl {L['law_field_curl']:+.2f}) left at rest at the bounce")
c.plot(tt, E2 / E2[0], "s", color=BL, ms=3.5, mfc="none", label="a generic field's energy / start")
c.plot(ts, 1 / np.cosh(ts), color=GR, lw=1.2, ls="--", label="1/a")
c.axhline(0, color=GR, lw=0.5)
c.text(0.97, 0.30, f"level 1 mesh, {L['edges']} edges\ncosmic-time run vs frozen-sphere run in η: {L['cosmic_vs_conformal']:.0e}\nE·a drift {L['E_times_a_drift']:.0e}; Gauss {L['gauss']:.0e}; gradient kept {L['gauge_kept']:.0e}",
       transform=c.transAxes, ha="right", fontsize=8.5, color=GR)
c.set_xlabel("time τ since the bounce"); c.set_ylabel("amplitude, or energy / start")
c.legend(fontsize=8.5, loc="upper right"); c.set_title("C. Discrete Maxwell on the law's mesh, slices growing (M4)")

# D. matter: momentum falls as 1/a
d = ax[3]
for r in R["matter"]:
    tau = np.array(r["tau"]); p = np.array(r["p"])
    d.semilogy(tau, p, "o", color=OR, ms=4)
    tg = np.linspace(0.05, 6, 300); d.semilogy(tg, np.sinh(r["beta"]) / np.cosh(tg), color=GR, lw=1)
    d.text(0.15, math.sinh(r["beta"]) * 1.12, f"β = {r['beta']}", color=OR, fontsize=9)
d.plot([], [], "o", color=OR, label="momentum γv of a particle on the law's tilted time line (rest frame of the slices)")
d.plot([], [], color=GR, label="sinh β / cosh τ: momentum × a fixed (textbook)")
d.set_xlabel("time τ since the bounce"); d.set_ylabel("momentum per unit mass, γv (log)")
d.legend(fontsize=8.5, loc="lower left"); d.set_title("D. Matter: the same 1/a as light (M5)")

fig.suptitle("DRAFT round 22: light, then matter, in the law's de Sitter spacetime", fontsize=13)
fig.tight_layout(rect=(0, 0, 1, 0.97))
fig.savefig("figures/f29_light_and_matter_in_de_sitter.png", bbox_inches="tight")
