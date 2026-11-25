import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, os
HERE = os.path.dirname(os.path.abspath(__file__))
fig, axes = plt.subplots(2, 1, figsize=(10, 4.2))
for ax, (lab, done, todo, unit) in zip(axes, [("modules", 62, 57, "modules"), ("theorems (core library count)", 1750, 1269 + 10, "theorems")]):
    ax.barh([0], [done], color="#2e86c1"); ax.barh([0], [todo], left=[done], color="#f5b041")
    ax.text(done/2, 0, f"ported on the branch: {done}", ha="center", va="center", color="white", fontsize=10, weight="bold")
    ax.text(done + todo/2, 0, f"newer, not ported: {todo}", ha="center", va="center", fontsize=10, weight="bold")
    ax.set_yticks([]); ax.set_title(lab, loc="left", fontsize=10); ax.spines[["top", "right", "left"]].set_visible(False)
axes[1].set_xlabel("count (the 62 ported modules became 1,970 theorems in mathlib form)")
fig.suptitle("Mathlib edition: about half the library is ported; everything since the morning of Nov 25 is not", weight="bold")
fig.tight_layout(); fig.savefig(os.path.join(HERE, "fig1-port-coverage.png"), dpi=130)
