"""DRAFT pictures for round 6. Run after checks.py. Orange = known, blue = ours, grey = control."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OR, BL, GR = "#e07b00", "#1f6fc2", "#777777"
R = json.load(open("results.json"))
plt.rcParams.update({"font.size": 10.5, "figure.dpi": 130})
fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.8))
ym = R["yang_mills"]
# A: convergence of action and force
ar = np.array([r[0] for r in ym["action_rows"]]); ae = np.array([r[3] for r in ym["action_rows"]])
fr = np.array([r[0] for r in ym["force_rows"]]); fe = np.array(ym["force_errs"])
ax[0].loglog(ar, ae, "o-", color=BL, label="loop action vs Yang–Mills action")
ax[0].loglog(fr, fe, "s-", color=BL, alpha=0.6, label="force on a link vs Yang–Mills equation")
ax[0].loglog(ar, ae[0] * (ar / ar[0])**2, ":", color=OR, label="slope 2 (lattice spacing²)")
ax[0].set_xlabel("loop size a"); ax[0].set_ylabel("relative difference")
ax[0].set_title("A. Strong force: Wilson's loops → Yang–Mills\naction and field equations (SU(3), 4D)")
ax[0].legend(fontsize=8)
# B: the k.ab-only field is pushed by its own commutator
cf = np.array(ym["comm_field_loops"]); cd = ym["sign"] * np.array(ym["comm_field_DF"])
xx = np.arange(8)
ax[1].bar(xx - 0.2, cd, 0.4, color=OR, label="Yang–Mills: −g²[A,[A,A]]-type term")
ax[1].bar(xx + 0.2, cf, 0.4, color=BL, label="force read from the 6 loops at a link")
ax[1].axhline(0, color=GR, lw=0.8)
ax[1].set_xticks(xx); ax[1].set_xticklabels([f"{i+1}" for i in xx]); ax[1].set_xlabel("colour direction")
ax[1].set_ylabel("force (× g/a³)")
ax[1].set_title("B. Rounds 4–5's field made only of k·ab\nis not at rest: its own twist pushes it")
ax[1].legend(fontsize=8)
# C: Friedmann
fr_ = R["friedmann"]; sw = np.array(fr_["sweep"])
hh = np.linspace(0.55, 1.45, 100); Lam = fr_["Lambda"]; H0 = np.sqrt(Lam / 3)
tg = np.linspace(0, 2, 801); bump = np.exp(-((tg - 1.0) / 0.25)**2)
curve = [np.trapezoid(np.exp(3 * f * H0 * tg) * (6 * (f * H0)**2 - 2 * Lam) * bump, tg) for f in hh]
ax[2].plot(hh, curve, "-", color=OR, label="Friedmann: 6aȧ² − 2Λa³ (bump-weighted)")
ax[2].plot(sw[:, 0], sw[:, 1], "o", color=BL, ms=8, label="varying the action read from loops")
ax[2].axhline(0, color=GR, lw=0.8); ax[2].axvline(1, color=GR, lw=0.8, ls=":")
ax[2].set_xlabel("expansion rate / √(Λ/3)"); ax[2].set_ylabel("δS / δ(clock rate)")
ax[2].set_title("C. Gravity: vary the loop curvature as an action →\nFriedmann; only H = √(Λ/3) is at rest")
ax[2].legend(fontsize=8)
fig.tight_layout(); fig.savefig("figures/f13_field_equations_from_loops.png", bbox_inches="tight")
print("f13 written")
