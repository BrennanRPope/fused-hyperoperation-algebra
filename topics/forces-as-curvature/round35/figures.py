"""DRAFT round 35 picture: f42_wave_energy.png, from results.json."""
import json, os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
R = json.load(open(os.path.join(HERE, "results.json")))
GREY, BLUE, ORANGE, RED, PURPLE, GREEN = "#777777", "#1f6fc5", "#e08a1e", "#c8372d", "#7b4fa8", "#2e9e4f"
fig = plt.figure(figsize=(20, 13.5))
gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.18)
fig.suptitle("Round 35 (DRAFT): the energy gravitational waves carry, and how it bends spacetime", fontsize=16)

# A: speed of crests and of energy on the pieces
ax = fig.add_subplot(gs[0, 0])
def cone(ks): return 2 * np.arcsinh(np.sqrt(np.sum(np.sin(np.asarray(ks) / 2)**2)))
kk = np.linspace(0.01, 1.05, 200)
for nm, v, c in (("along an axis", [1, 0, 0], BLUE), ("along a body diagonal", [1, 1, 1], PURPLE)):
    d = np.array(v, float) / np.linalg.norm(v)
    ph = [cone(x * d) / x for x in kk]; gr = [(cone((x + 1e-6) * d) - cone((x - 1e-6) * d)) / 2e-6 for x in kk]
    ax.plot(kk, ph, "--", color=c, lw=1.5, label=f"wave crests, {nm}")
    ax.plot(kk, gr, "-", color=c, lw=2.5, label=f"energy, {nm}")
ax.axhline(1, color=GREY, ls=":", lw=1); ax.text(0.02, 1.006, "light speed", color=GREY, fontsize=11)
worst = max(max(abs(r["plus"] - 1), abs(r["cross"] - 1)) for r in R["E1"]["rows"] if r["k"] == 0.005)
ax.text(0.03, 0.835, f"+ and × carry energy alike:\nthe pieces' stiffness for each is k²/4 to {worst:.0e}", color=ORANGE, fontsize=11.5)
ax.set_ylim(0.74, 1.02); ax.set_xlabel("wavenumber × size of a piece, ka"); ax.set_ylabel("speed (light = 1)")
ax.set_title("A. On the pieces, a wave's energy travels at light speed at long wavelengths (E1)\nnear one piece across it lags 3 times as much as the crests do")
ax.legend(fontsize=10, loc="lower left")

# B: where a wave's energy goes on a slice
ax = fig.add_subplot(gs[0, 1])
ax.bar([0], [0.25], color=ORANGE, width=0.55, label="the wave's motion  ⟨K_ij K^ij⟩")
ax.bar([0], [0.25], bottom=[0.25], color=BLUE, width=0.55, label="its ripples, which bend the slice less on average  −⟨³R⁽²⁾⟩")
ax.bar([1], [0.5], color=GREEN, width=0.55, label="the smooth bending they demand underneath  16π ρ_gw")
ax.text(0.5, 0.25, "=", fontsize=40, ha="center", va="center", color=GREY)
ax.set_xticks([0, 1]); ax.set_xticklabels(["the wave, averaged\nover a wavelength", "the slice's\nbackground bending"], fontsize=12)
ax.set_ylabel("in units of ω² (A² + B²)  (A, B = + and × amplitudes)"); ax.set_ylim(0, 0.95); ax.set_xlim(-0.6, 1.6)
ax.text(0.68, 0.58, "law's local dial rises by 8π ρ_gw / 3:\nthe same rule as matter (rounds 23, 25, 30)", color=PURPLE, fontsize=11.5)
ax.text(-0.55, 0.58, "in spacetime: curvature of a beam\nof light-like matter, energy = flux,\nno sideways pressure (Isaacson)", color=GREY, fontsize=11)
ax.set_title("B. A wave's energy bends spacetime the way matter does (E2, exact)")
ax.legend(fontsize=10, loc="upper left")

# C: orbiting pairs, power against eccentricity
ax = fig.add_subplot(gs[1, 0])
ee = np.linspace(0, 0.85, 300); f = (1 + 73 / 24 * ee**2 + 37 / 96 * ee**4) / (1 - ee**2)**3.5
ax.semilogy(ee, f, color=BLUE, lw=2.5, label="Peters–Mathews formula")
kep = R["E3"]["kepler"]; f0 = kep[0]["peters_mathews"]
ax.semilogy([r["e"] for r in kep], [r["numeric"] / f0 for r in kep], "o", color=ORANGE, ms=11, label="energy flux summed over a sphere, orbit worked out numerically")
ax.axvline(0.617134, color=RED, ls=":", lw=1.5); ax.text(0.625, 2.0, "Hulse–Taylor\npulsar, e = 0.617\n(about 12 ×)", color=RED, fontsize=11)
ax.set_xlabel("orbit eccentricity e"); ax.set_ylabel("power radiated, relative to a circular orbit")
ax.set_title("C. The energy an orbiting pair sends off (E3)\nsumming E1's flux over a sphere gives the quadrupole formula; stretched orbits radiate far more")
ax.legend(fontsize=10.5, loc="upper left")

# D: Hulse-Taylor
ax = fig.add_subplot(gs[1, 1])
yrs = 1975 + np.array(R["E4"]["years"]); sh = np.array(R["E4"]["shift_s"])
ax.plot(yrs, sh, color=BLUE, lw=2.5, label="predicted from E3's power and the measured orbit")
ax.fill_between(yrs, sh * (0.9983 - 0.0016), sh * (0.9983 + 0.0016), color=RED, alpha=0.35, label="measured: 0.9983 ± 0.0016 of the prediction")
ax.axhline(0, color=GREY, ls=":", lw=1); ax.text(1998, 2.5, "an orbit that did not radiate", color=GREY, fontsize=11); ax.set_ylim(-80, 9)
ax.set_xlabel("year"); ax.set_ylabel("periastron arrives earlier by (seconds)")
ax.set_title("D. The Hulse–Taylor binary pulsar's orbit shrinks as predicted (E4)\nperiod falls 2.40 × 10⁻¹² s per s; about 8 × 10²⁴ W goes out as gravitational waves")
ax.text(1976, -48, "measured timing points are not re-plotted here;\nthe red band (too thin to see here) is the\npublished measurement (Weisberg and Huang 2016)", color=GREY, fontsize=10.5)
ax.legend(fontsize=10.5, loc="lower left")
fig.savefig(os.path.join(HERE, "figures", "f42_wave_energy.png"), dpi=110, bbox_inches="tight")
