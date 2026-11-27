"""DRAFT figures for the prior-art audit. Run: python3 figures.py  (matplotlib)"""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from collections import Counter
from claims import C, AREAS

COL = {"new": "#2f6fd6", "partly": "#f08c28", "print": "#a3a3a3"}
NAME = {"new": "new as far as searched", "partly": "partly known", "print": "already in print"}
PRED = {"out": ("#d62728", "-", "prediction, ruled out"), "nolead": ("#7f7f7f", "-.", "worked, no lead (round 28)"),
        "math": ("#111111", "--", "math yes/no check"), "far": ("#8c4fb3", ":", "fixed, out of reach")}

# Figure 1: the claim map
fig, ax = plt.subplots(figsize=(19, 13)); ax.set_xlim(0, 5); ax.set_ylim(-2.6, 23.6); ax.axis("off")
for a, area in enumerate(AREAS):
    items = [c for c in C if c[1] == a]
    order = {"new": 0, "partly": 1, "print": 2}
    items.sort(key=lambda c: order[c[3]])
    ax.text(a + 0.5, 23.2, area, ha="center", va="center", fontsize=12.5, weight="bold")
    for i, (cid, _, label, st, conf, pred) in enumerate(items):
        y = 22.3 - i * 1.0
        ec, ls, lw = ("white", "-", 0.5)
        if pred: ec, ls, _ = PRED[pred]; lw = 3.2
        box = FancyBboxPatch((a + 0.04, y - 0.4), 0.92, 0.8, boxstyle="round,pad=0.01,rounding_size=0.08",
                             fc=COL[st], ec=ec, lw=lw, ls=ls)
        ax.add_patch(box)
        txt = f"{cid}  {label}" + ("  (?)" if conf == "low" else "")
        ax.text(a + 0.08, y, txt, va="center", fontsize=9.6, color="white" if st != "print" else "#222")
hs = [plt.Rectangle((0, 0), 1, 1, fc=COL[k]) for k in COL] + \
     [plt.Rectangle((0, 0), 1, 1, fc="white", ec=v[0], lw=3, ls=v[1]) for v in PRED.values()]
ax.legend(hs, [NAME[k] for k in COL] + [v[2] for v in PRED.values()], loc="lower right", ncol=2, fontsize=11,
          frameon=False, bbox_to_anchor=(1.0, 0.0))
ax.text(0.0, -1.2, "(?) = low confidence: thin search, or sources seen only as abstracts.", fontsize=10.5, color="#444")
ax.set_title("Figure 1. Every claim the repo marks as ours, coloured by what the literature search found (73 claims)",
             fontsize=14, loc="left")
plt.tight_layout(); plt.savefig("figures/fig1-claim-map.png", dpi=110); plt.close()

# Figure 2: counts by area
fig, ax = plt.subplots(figsize=(10, 4.6))
ys = range(len(AREAS))[::-1]
for y, a in zip(ys, range(len(AREAS))):
    cnt = Counter(c[3] for c in C if c[1] == a); left = 0
    for st in ["new", "partly", "print"]:
        n = cnt.get(st, 0)
        if n:
            ax.barh(y, n, left=left, color=COL[st]); ax.text(left + n/2, y, str(n), ha="center", va="center",
                                                             color="white" if st != "print" else "#222", fontsize=11)
        left += n
ax.set_yticks(list(ys)); ax.set_yticklabels(AREAS, fontsize=11); ax.set_xlabel("claims")
for s in ["top", "right"]: ax.spines[s].set_visible(False)
ax.legend([plt.Rectangle((0, 0), 1, 1, fc=COL[k]) for k in COL], [NAME[k] for k in COL], frameon=False, loc="lower right")
ax.set_title("Figure 2. Where the new pieces are: mostly pure mathematics; none in the one-law physics program",
             loc="left", fontsize=12)
plt.tight_layout(); plt.savefig("figures/fig2-by-area.png", dpi=120); plt.close()

# Figure 3: the closest prior programs
fig, ax = plt.subplots(figsize=(13, 7.4)); ax.set_xlim(0, 13); ax.set_ylim(0, 7.4); ax.axis("off")
def box(x, y, w, h, title, body, fc):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.05,rounding_size=0.15", fc=fc, ec="#555"))
    ax.text(x + 0.15, y + h - 0.3, title, weight="bold", fontsize=11, va="top")
    ax.text(x + 0.15, y + h - 0.75, body, fontsize=9.6, va="top", linespacing=1.35)
box(4.3, 2.85, 4.4, 1.85, "This project's physics rounds",
    "one law u⊕v = (u+v)(1 − q ū v)⁻¹ with dial q:\nS³ at q>0, boosts at q<0, de Sitter with a\nLorentzian dot, gauge fields as loops", "#dbe7fb")
nb = [(0.2, 5.0, "Nesterov & Sabinin 2000–2004", "same quaternion loop with a curvature\nconstant; S², de Sitter, FRW slices;\ncurvature = associator holonomy", "#eeeeee"),
      (8.8, 5.0, "κ-stereographic model (2020)", "Bachmann–Bécigneul–Ganea; Skopek et al.:\nthe identical 3D formula, κ = q,\nboth signs; no physics", "#eeeeee"),
      (0.2, 0.3, "Ungar's gyrovector spaces", "Einstein / Möbius addition (q<0),\nThomas–Wigner rotation as gyration,\nangle = defect", "#eeeeee"),
      (8.8, 0.3, "Relative locality / DSR (2011)", "a composition law of momenta makes\ncurvature and torsion of momentum space;\nmakes testable time-delay predictions", "#eeeeee")]
for x, y, t, b, fc in nb:
    box(x, y, 4.0, 2.0, t, b, fc)
    cx = x + 4.0 if x < 4 else x
    cy = y + 1.0
    tx = 4.3 if x < 4 else 8.7
    ty = 4.5 if y > 3 else 3.0
    ax.plot([cx, tx], [cy, ty], color="#888", lw=1.5)
ax.text(4.45, 2.45, "Only here (all rebuild known physics):\n• tie to the scalar law a+b+kab and its fused shift\n• Regge slabs on corners the law grows\n• the q=+1 twist read as the Hopf photon",
        fontsize=10, va="top", color="#1f4e9e")
ax.set_title("Figure 3. The one-law program's nearest neighbours: its skeleton is Nesterov–Sabinin's, its formula is κ-addition",
             loc="left", fontsize=12)
plt.tight_layout(); plt.savefig("figures/fig3-neighbours.png", dpi=120); plt.close()
print("figures written")
