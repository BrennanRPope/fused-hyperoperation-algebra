"""Figures for the econ-forecasting page. Run: python3 figures.py (writes figures/*.png)"""
import os, json, warnings
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import load
warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "figures")
INK, MUTE, GRID = "#15131a", "#5d5a66", "#e4e3ea"
BLUE, RED, GREEN, AMBER = "#1d2bd8", "#cf3519", "#1b8a5a", "#b07a00"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11, "axes.edgecolor": MUTE,
    "axes.labelcolor": INK, "xtick.color": MUTE, "ytick.color": MUTE, "axes.spines.top": False,
    "axes.spines.right": False, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
    "lines.linewidth": 2, "savefig.facecolor": "white", "figure.facecolor": "white"})
def save(fig, name): fig.savefig(os.path.join(OUT, name), dpi=170, bbox_inches="tight"); plt.close(fig)

# ---------- e1: average the corners, not the rates ----------
fig, (a, b) = plt.subplots(1, 2, figsize=(11, 4.3))
T = np.arange(0, 201)
rates = [0.005, 0.025, 0.045]
facs = [(1 + r) ** -T for r in rates]
avg = np.mean(facs, axis=0)
for r, f, c in zip(rates, facs, [GREEN, MUTE, AMBER]):
    a.plot(T, f, color=c, lw=1.6, label=f"{100*r:.1f}% for ever")
a.plot(T, avg, color=BLUE, lw=2.6, label="average of the three factors")
a.plot(T, (1 + np.mean(rates)) ** -T, color=RED, ls="--", lw=2, label="factor at the average rate (2.5%)")
a.set_yscale("log"); a.set_ylim(1e-4, 1.2); a.set_xlabel("years ahead"); a.set_ylabel("discount factor (log scale)")
a.set_title("Three equally likely rate worlds", fontsize=11.5); a.legend(fontsize=9, frameon=False, loc="lower left")
ce = avg[1:] ** (-1 / T[1:]) - 1
b.plot(T[1:], 100 * ce, color=BLUE, lw=2.6, label="rate implied by the averaged factor")
b.axhline(2.5, color=RED, ls="--", lw=2, label="average rate")
b.axhline(0.5, color=GREEN, lw=1.4, label="lowest rate")
b.set_ylim(0, 3); b.set_xlabel("years ahead"); b.set_ylabel("discount rate, %")
b.set_title("The long-run rate falls toward the lowest corner", fontsize=11.5); b.legend(fontsize=9, frameon=False, loc="center right", bbox_to_anchor=(1, 0.62))
save(fig, "e1_average_corners.png")

# ---------- e2: backtest ----------
J = json.load(open(os.path.join(HERE, "discount_backtest.json")))
bt = J["backtest"]; hs = sorted(bt, key=int)
names = [("constant mean rate", RED, "constant average rate (baseline)"), ("averaged factors, block 1", AMBER, "averaged factors, years shuffled"),
         ("averaged factors, block 10", GREEN, "averaged factors, 10-year blocks"), ("averaged, uncertain mean", BLUE, "averaged factors, uncertain long-run mean")]
fig, (a, b) = plt.subplots(1, 2, figsize=(11, 4.3))
w = 0.2; x = np.arange(len(hs))
for i, (k, c, lab) in enumerate(names):
    a.bar(x + (i - 1.5) * w, [bt[h][k]["mae_log"] for h in hs], w, color=c, label=lab)
    b.plot([int(h) for h in hs], [bt[h][k]["ratio"] for h in hs], "o-", color=c, label=lab)
a.set_xticks(x, [f"{h} yrs" for h in hs]); a.set_ylabel("mean |log error| (lower is better)")
a.set_title("Forecasting the discount factor history delivered", fontsize=11.5); a.legend(fontsize=8.5, frameon=False, loc="upper left")
b.axhline(1, color=INK, lw=1); b.set_xlabel("horizon, years"); b.set_ylabel("realized ÷ forecast (1 = unbiased)")
b.set_title("Bias: above 1 means the future was undervalued", fontsize=11.5)
save(fig, "e2_backtest.png")

# ---------- e3: the square cut from the corner, and the Gordon pole ----------
fig, (a, b) = plt.subplots(1, 2, figsize=(11, 4.6), gridspec_kw={"width_ratios": [1, 1.25]})
m, s = 0.35, 0.22
a.add_patch(Rectangle((0, 0), 1 + m, 1 + m, fill=False, ec=MUTE, lw=1.5, ls="--"))
a.add_patch(Rectangle((0, 0), 1 + m + s, 1 + m - s, color=BLUE, alpha=0.18, ec=BLUE, lw=2))
a.add_patch(Rectangle((1 + m, 1 + m - s), s, s, color=RED, alpha=0.35, ec=RED, lw=1.5))
a.add_patch(Rectangle((1 + m, 0), s, 1 + m - s, color=GREEN, alpha=0.25, ec=GREEN))
a.text((1 + m) / 2, (1 + m - s) / 2, "(1+m+s)(1+m−s)\nwhat you end with", ha="center", va="center", fontsize=10, color=INK)
a.text(1 + m + s / 2, 1 + m - s / 2, "s²\nshort", ha="center", va="center", fontsize=10, color=RED)
a.text((1 + m) / 2, 1 + m - s / 2, "this strip, (1+m) by s ...", ha="center", va="center", fontsize=9, color=MUTE)
a.text(1 + m + s / 2, (1 + m - s) / 2, "... moves here, but is s too long", ha="center", va="center", fontsize=8, color=GREEN, rotation=90)
a.text(0.02, 1 + m + 0.03, "the box you'd have with steady growth (1+m)²", fontsize=9, color=MUTE)
a.set_xlim(-0.05, 1.75); a.set_ylim(-0.05, 1.55); a.set_aspect("equal"); a.axis("off")
a.set_title("Volatility drag: up s, down s", fontsize=11.5)
r = 0.07; g = np.linspace(-0.02, 0.0685, 400)
P = (1 + g) / (r - g)
b.plot(100 * g, P, color=BLUE); b.axvline(100 * r, color=RED, ls="--", lw=1.6); b.text(100 * r - 0.1, 260, "pole at g = r", color=RED, ha="right")
b.set_ylim(0, 300); b.set_xlabel("dividend growth g, % (discount rate r = 7%)"); b.set_ylabel("price ÷ dividend")
mm = load.shiller_monthly()
for y, lab in (("1950-01", "1950"), ("1982-01", "1982"), ("2009-03", "2009"), ("2022-12", "2022"), ("2000-01", "2000")):
    row = mm.loc[y].iloc[0]; pd_ = row.P / row.D
    gg = (r - 1 / pd_) / (1 + 1 / pd_)      # g that makes D/(r (-) g) = P/D
    b.plot(100 * gg, pd_, "o", color=INK, ms=5); b.annotate(lab, (100 * gg, pd_), textcoords="offset points", xytext=(-26, 2), fontsize=9)
b.set_title("A growing stream is a perpetuity at r ⊖ g", fontsize=11.5)
save(fig, "e3_drag_and_pole.png")

# ---------- e4: effective dimension of the S&P 500 ----------
D = pd.read_csv(os.path.join(HERE, "dimension_SP500_const.csv"), index_col=0, parse_dates=True)
fig, a = plt.subplots(figsize=(11, 3.9))
a.plot(D.index, D.d_eff, color=BLUE, lw=1.6, label="independent axes (effective dimension)")
a.set_ylabel("effective dimension", color=BLUE)
a2 = a.twinx(); a2.plot(D.index, 100 * D.rv22, color=RED, lw=1, alpha=0.6); a2.set_ylabel("trailing 1-month volatility, %", color=RED)
a2.spines["right"].set_visible(True); a2.grid(False)
a.set_title(f"S&P 500 members: {int(D.n.median())} stocks behave like about {D.d_eff.median():.0f} independent ones", fontsize=11.5)
save(fig, "e4_dimension.png")
print("figures written")
