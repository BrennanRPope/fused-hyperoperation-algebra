"""DRAFT picture for round 15. Run after checks.py. Blue = the law, orange = Regge read-back, grey = reference."""
import json, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OR, BL, GR = "#e07b00", "#1f6fc2", "#777777"
R = json.load(open("results.json"))
plt.rcParams.update({"font.size": 10.5, "figure.dpi": 130})
fig, ax = plt.subplots(1, 3, figsize=(16.5, 5.2))

a = ax[0]
rows = R["boost_loops"]; lv = [r["level"] for r in rows]
a.semilogy(lv, [r["eps_max"] for r in rows], "o-", color=OR, ms=8, label="largest boost deficit at a hinge")
a.semilogy(lv, [r["twist_max"] for r in rows], "s-", color=BL, ms=8, label="largest twist going round two hinges")
a.semilogy(lv, [max(r["err2"], r["err4"]) for r in rows], "^:", color=GR, ms=7, label="loop vs the law at q = −1 (mismatch)")
a.set_xticks(lv); a.set_xticklabels([f"{l}\n(time step {0.2 / 2**l:g})" for l in lv])
a.set_xlabel("refinement level of the slabs"); a.set_ylabel("radians (rapidity for the deficit)")
a.set_ylim(1e-17, 1); a.legend(fontsize=8, loc="center right")
a.text(1, 3e-15, "the law composes every loop\nto rounding (Ungar's identity)", ha="center", fontsize=8.5, color=GR)
a.set_title("A. Boost loops on the glued slabs:\nthe law composes them; the twist fades fastest")

b = ax[1]
r3 = R["shift_read"]["3"]; kap, F, G, t0, tau = r3["kappa"], r3["F"], r3["G"], r3["t0"], r3["tau"]
tt = np.linspace(-1.6, 1.6, 300); n = (tt - t0) / tau; Rr = math.exp(kap * tau)
b.plot(tt, F * Rr**n, "--", color=BL, label="forward run of one step: (1 + kτ)ⁿ")
b.plot(tt, G * Rr**-n, ":", color=BL, lw=2, label="run of its inverse: (1 + kτ)⁻ⁿ")
b.plot(tt, F * Rr**n + G * Rr**-n, color=GR, lw=1.5, label="their sum")
ta = t0 + tau * np.arange(len(r3["a"]))
b.plot(ta, r3["a"], "o", color=OR, ms=5, label="glued slabs (level 3), started at t = 0.5")
b.plot([r3["bounce_t"]], [r3["throat"]], "*", color="k", ms=13)
b.annotate(f"the runs balance: the bounce\nt = {r3['bounce_t']:+.3f}, size {r3['throat']:.4f} ≈ 1/k", (r3["bounce_t"], r3["throat"]), (-1.55, 0.62),
           fontsize=8.5, arrowprops=dict(arrowstyle="->", color=GR))
b.axhline(0, color="k", lw=0.6); b.text(1.55, 0.06, "size 0 = the absorbing centre −1/k", ha="right", fontsize=8, color=GR)
b.set_xlabel("proper time t (units of 1/k = √(3/Λ))"); b.set_ylabel("size of the universe a")
b.set_ylim(-0.1, 3); b.legend(fontsize=8, loc="upper right")
b.set_title("B. Time as the fused shift: the slabs are one\nstep's run forward plus its inverse's run")

c = ax[2]
for L, mk in (("1", "^"), ("2", "s"), ("3", "o")):
    d = R["rate_drift"][L]
    c.plot(d["piece"], d["rate"], mk + "-", color=OR, ms=4, mfc="white" if L != "3" else OR, label=f"level {L}")
c.axhline(1, color=GR, lw=1.5, label="the shift's rate √(Λ/3)")
c.axvline(1, color=BL, ls=":", lw=1.5); c.text(0.115, 2.6, "dotted line: pieces = 1/k,\nthe absorbing centre's size", color=BL, fontsize=8.5)
c.set_xscale("log"); c.set_xlabel("size of one space piece, in units of 1/k"); c.set_ylabel("rate read from 3 slices: arccosh(...)/τ")
c.legend(fontsize=8, loc="upper left")
c.set_title("C. The slabs keep the shift's rate until\ntheir pieces near 1/k, then leave it")
fig.tight_layout(); fig.savefig("figures/f22_boost_loops_and_time_as_the_shift.png", bbox_inches="tight")
print("f22 written")
