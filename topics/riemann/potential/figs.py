import json, math
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
O = json.load(open("potential.json"))
R = np.load("ring_traj.npz"); T = np.load("toy_traj.npz")
snaps = np.load("../known-tools/rewind_snaps.npz")

# ---------------- Fig 1: the snag -------------------------------------------------------------------
fig, ax = plt.subplots(1, 3, figsize=(17, 5.2))
a = ax[0]
ex = O["excess"]; ts = [0, 0.01, 0.05, 0.2]
a.plot(ts, [ex["ours"][str(float(t))] for t in ts], "o-", color="C3", label="our zeros (1000–4500) under heat flow")
a.axhline(ex["GUE"], color="C0", ls="--", lw=1); a.text(0.12, ex["GUE"] + 0.0006, "random-matrix (GUE) spacing", color="C0", fontsize=8)
a.axhline(0, color="k", lw=1); a.text(0.12, 0.0006, "evenly spaced crystal = the lowest energy", fontsize=8)
a.set_ylim(-0.002, 0.022)
a.set_xlabel("heat time t"); a.set_ylabel("excess log energy over the crystal (per zero, per neighbour)")
a.set_title("1a. The snag, measured: the zeros' potential energy only falls\nas t goes forward. Its bottom is the crystal at t → ∞, not t = 0", fontsize=10)
a.annotate("t = 0: nothing special,\njust part way down", xy=(0, ex["ours"]["0.0"]), xytext=(0.05, 0.017),
           arrowprops=dict(arrowstyle="->", lw=0.8), fontsize=8)
a.text(0.075, 0.0035, f"(Poisson, no repulsion: {ex['Poisson']:.3f}, off the top)", fontsize=7.5, color="0.4")
a.legend(fontsize=8, loc="lower left")

a = ax[1]
a.plot(R["ttr"], R["Wrr"], color="C1", lw=2, label="off-line pair being pulled in")
a.plot(R["ttr2"], R["Wrr2"], color="C2", lw=2, label="after capture: a real pair that repels")
a.axhline(O["ring_toy"]["W_crystal"], color="k", ls=":", lw=1)
a.text(0.34, O["ring_toy"]["W_crystal"] - 0.75, "crystal (lowest energy\nfor zeros on the line)", fontsize=8)
a.axvline(O["ring_toy"]["t_capture"], color="0.5", lw=0.8)
a.text(O["ring_toy"]["t_capture"] + 0.008, O["ring_toy"]["W_at_capture"] - 0.4, "capture:\nenergy spikes", fontsize=8, color="0.3")
a.set_xlabel("heat time t"); a.set_ylabel("log potential energy W of 39 zeros on a ring")
a.set_title("1b. Off the line the flow climbs the energy. The off-line start sits\nBELOW the crystal; capture is uphill, then it rolls back down", fontsize=10)
a.legend(fontsize=8, loc="upper right")

a = ax[2]
z0, zc, ze = R["z0"], R["zcap"], R["zend"]
a.plot(z0.real, z0.imag, "o", color="C1", ms=5, label="t = 0: a pair at height ±0.6 off the line")
a.plot(zc.real, zc.imag - 1.2, "o", color="C3", ms=5, label=f"t = {O['ring_toy']['t_capture']:.3f}: the pair lands (collision)")
a.plot(ze.real, ze.imag - 2.4, "o", color="C2", ms=5, label="t = 0.62: all real, relaxing to the crystal")
for yy in (0, -1.2, -2.4): a.axhline(yy, color="0.8", lw=0.6)
a.set_xlim(12, 28); a.set_ylim(-3.1, 1.2); a.set_yticks([])
a.set_xlabel("position on the ring (zeros 13–27 of 39 shown)")
a.set_title("1c. The ring toy behind 1b (rows offset for clarity)", fontsize=10)
a.legend(fontsize=8, loc="upper left")
plt.tight_layout(); plt.savefig("figures/fig1_snag.png", dpi=130); plt.close()

# ---------------- Fig 2: where the edge shows -----------------------------------------------------------
fig, ax = plt.subplots(1, 3, figsize=(17, 5.2))
a = ax[0]
L = O["lehmer_pair"]; S, d0, tc = L["S"], L["d0"], L["t_c"]
tt = np.linspace(tc - 3e-3, tc + 3e-3, 2001)
u = 2 / S + (d0**2 - 2 / S) * np.exp(-4 * S * tt)
a.plot(tt, -0.5 * np.log(np.abs(u)), color="C3")
a.axvline(tc, color="C3", ls=":", lw=0.8); a.axvline(0, color="k", lw=1)
a.text(tc - 2.9e-3, -0.5 * np.log(abs(u[0])) + 1.6, "rewound past here the pair is\noff the line (complex)", fontsize=8)
a.text(1e-4, -0.5 * np.log(abs(u[-1])) + 1.6, "t = 0: the pair\nis real", fontsize=8)
a.set_xlabel("heat time t"); a.set_ylabel("pair energy  −log|gap|")
a.set_title(f"2a. The closest pair we have (γ ≈ {L['gamma_mid']:.0f}): its energy spikes\nat the collision, t_c = {tc:.2e}, just before 0", fontsize=10)

a = ax[1]
late = np.array(O["latest_collision_by_height"])
a.semilogy(late[:, 0], -late[:, 1], "o-", color="C3", label="latest collision (energy spike) in each height band, t < 0")
a.set_xlabel("height γ"); a.set_ylabel("how far before t = 0 the last spike is  (−t)")
a.set_title("2b. Λ is the time of the LAST energy spike, not a minimum.\nOurs all sit before 0; the closest gets nearer 0 higher up (Rodgers–Tao)", fontsize=10)
a.text(400, 2e-4, f"An off-line zero (β = 1) above 3·10¹² would put its spike AFTER 0,\nat t ≈ {O['offline_world_tcap_3e12']:.3f} (mean-field estimate, 'squeeze' page)", fontsize=8, color="C1")
a.legend(fontsize=8, loc="upper right"); a.set_ylim(1e-4, 0.3)

a = ax[2]
lt = np.array(O["li_toy"]["t"]); lm = np.array(O["li_toy"]["min_lambda"])
la = np.array(O["li_toy"]["t_after"]); lma = np.array(O["li_toy"]["min_lambda_after"])
a.semilogy(lt, -lm, "o", color="C1", label="before capture: the smallest Li number is NEGATIVE (size shown)")
a.semilogy(la, lma, "s", color="C2", label="after capture: every Li number is positive (smallest shown)")
tcap = O["offline_toy"]["t_capture"]
a.axvline(tcap, color="0.5", lw=0.8); a.text(tcap + 0.01, 1e8, f"capture\nt = {tcap:.3f}", fontsize=8)
a.set_xlabel("heat time t (toy: 76 zeros with one off-line pair, mirrored)")
a.set_ylabel("|smallest Li number λ_n, n ≤ 3000|")
a.set_title("2c. The 'trend against wiggle' balance (Li numbers) flips exactly\nat the edge, wherever the edge is. Putting it at t = 0 is RH", fontsize=10)
a.legend(fontsize=8, loc="upper right")
plt.tight_layout(); plt.savefig("figures/fig2_edge.png", dpi=130); plt.close()

# ---------------- Fig 3: what an argument would need ------------------------------------------------------
fig, a = plt.subplots(figsize=(15, 4.6)); a.axis("off")
boxes = [
 (0.01, "1. An energy defined for zeros\nanywhere in the plane", "Have it: log energy W,\nor de Bruijn's off-line height²", "C2"),
 (0.26, "2. Its lowest state =\nall zeros on the line", "W fails: off-line is cheaper.\nHeight² works (zero iff all real)", "C1"),
 (0.51, "3. The flow lowers it", "On the line, W: yes (theorem).\nHeight²: yes (de Bruijn 1950)", "C2"),
 (0.76, "4. ζ itself, at t = 0, is\nalready in the lowest state", "This IS the Riemann hypothesis.\nNo flow energy can supply it", "C3"),
]
for x, h, b, c in boxes:
    a.add_patch(FancyBboxPatch((x, 0.30), 0.22, 0.52, boxstyle="round,pad=0.01", fc="white", ec=c, lw=2.5))
    a.text(x + 0.11, 0.70, h, ha="center", va="center", fontsize=11, weight="bold")
    a.text(x + 0.11, 0.45, b, ha="center", va="center", fontsize=9.5, color=c)
for x in (0.235, 0.485, 0.735):
    a.annotate("", xy=(x + 0.02, 0.56), xytext=(x - 0.003, 0.56), arrowprops=dict(arrowstyle="-|>", lw=1.5))
a.text(0.5, 0.13, "Steps 1–3 only reach the lowest state at t = Λ, the edge. They say the edge exists, and Rodgers–Tao's version of the same\n"
       "energy argument says it is not before 0 (Λ ≥ 0). Pinning it AT 0 needs step 4, a fact about the primes at t = 0, which no flow energy contains.",
       ha="center", fontsize=10)
a.set_title("3. What 'least energy pins Λ = 0' would need, link by link (green: in hand · orange: needs the right energy · red: missing)", fontsize=11)
a.set_xlim(0, 1); a.set_ylim(0, 0.95)
plt.savefig("figures/fig3_what_it_needs.png", dpi=130, bbox_inches="tight"); plt.close()
print("ok")
