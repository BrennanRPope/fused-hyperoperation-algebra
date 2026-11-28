"""Pictures for the clock test (DRAFT). Run after checks.py:  python3 figures.py
Colour key: blue = the law's own reading; orange = textbook physics; grey = reference lines."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

BLUE, ORANGE, GREY, RED = "#1f77b4", "#ff7f0e", "#7f7f7f", "#c0392b"
R = json.load(open("results.json"))
plt.rcParams.update({"font.size": 10, "figure.dpi": 130})


def save(fig, name):
    fig.savefig(f"figures/{name}.png", bbox_inches="tight")
    plt.close(fig)


# k1 has it been said before?
rows = ["Page-Wootters 1983", "Connes-Rovelli thermal time 1994", "Stueckelberg / Horwitz-Piron tau",
        "de Broglie internal clock 1924", "Klein 1926 (fifth-dimension circle)", "Bohm 1952 + Valentini 1991",
        "'t Hooft cellular automaton 2016", "Hestenes zitterbewegung 1990", "Wootters ubit 2013",
        "Wrapped-clock idea (Brennan)"]
cols = ["one clock\neverywhere", "deterministic\nat the root", "wrapping\ngives the i", "wrapping gives\nthe odds", "twist couples\nto the clock"]
M = np.array([[2, 0, 1, 0, 0],
              [2, 0, 0, 1, 0],
              [2, 0, 0, 0, 0],
              [0, 1, 1, 0, 0],
              [0, 1, 2, 0, 1],
              [1, 2, 0, 2, 0],
              [1, 2, 0, 1, 0],
              [0, 1, 2, 0, 1],
              [2, 0, 2, 0, 0],
              [2, 2, 2, 2, 2]])
fig, ax = plt.subplots(figsize=(9.5, 5.6))
cmap = {0: "white", 1: "#fde0c5", 2: ORANGE}
for i in range(len(rows)):
    for j in range(len(cols)):
        col = cmap[M[i, j]] if i < len(rows) - 1 else ("#cfe2f3" if M[i, j] else "white")
        ax.add_patch(plt.Rectangle((j, len(rows) - 1 - i), 1, 1, facecolor=col, edgecolor=GREY))
        mark = {0: "", 1: "partly", 2: "yes"}[M[i, j]]
        ax.text(j + 0.5, len(rows) - 0.5 - i, mark if i < len(rows) - 1 else "claims", ha="center", va="center", fontsize=8)
ax.set_xlim(0, len(cols)); ax.set_ylim(0, len(rows))
ax.set_xticks(np.arange(len(cols)) + 0.5); ax.set_xticklabels(cols, fontsize=8)
ax.set_yticks(np.arange(len(rows)) + 0.5); ax.set_yticklabels(rows[::-1], fontsize=8)
ax.xaxis.tick_top()
ax.set_title("k1 · Has it been said before? Each part has; no one source says all five", pad=40)
save(fig, "k1_prior_work")

# k2 a clock in place of the energy zero
fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
t = np.linspace(0, 7, 400)
for e, lab in zip([1, 3, 6], ["level 1", "level 3", "level 6"]):
    axs[0].plot(t, np.angle(np.exp(-1j * e * 2 * np.pi / 8 * t)), color=ORANGE, alpha=0.4 + 0.2 * (e > 1))
axs[0].set_xlabel("clock reading t"); axs[0].set_ylabel("system phase (wrapped)")
axs[0].set_title("A · Frozen whole, moving part (K1)\nread at each clock tick, the system turns as e^(-iHt)")
names = ["shared turn\n(balanced)", "separate turns\n(mirror part)"]
axs[1].bar([0, 1], [0.5, 0.5], color=[ORANGE, GREY])
axs[1].text(0, 0.52, "obeys clock + system = 0\nruns forward (error 5e-15)", ha="center", fontsize=8)
axs[1].text(1, 0.52, "breaks it; obeys\nclock - system = whole turns:\nruns backward", ha="center", fontsize=8)
axs[1].set_xticks([0, 1]); axs[1].set_xticklabels(names); axs[1].set_ylim(0, 0.8)
axs[1].set_ylabel("share of a real-number history")
axs[1].set_title("B · Real numbers, separate quarter turns (K3)\nhalf the history is a mirror world")
save(fig, "k2_clock_replaces_energy_zero")

# k3 gravity as the clock, and the law's slice as the switch
fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
errs = R["K4"]["errs"]
Ms = np.array([float(k) for k in errs]); E = np.array(list(errs.values()))
axs[0].loglog(Ms, E, "o-", color=ORANGE, label="matter state vs e^(-iH tau)")
axs[0].loglog(Ms, E[0] * 10 / Ms, "--", color=GREY, label="1/M")
axs[0].set_xlabel("heaviness M of the gravitational clock"); axs[0].set_ylabel("error")
axs[0].legend(fontsize=8); axs[0].set_title("A · The i comes from gravity's phase direction (K4)")
eps = 0.02
for q, col in zip([0.25, 1.0, 4.0], [BLUE, "#5aa0d8", "#9cc6ea"]):
    at = 1 / np.sqrt(q)
    sol = solve_ivp(lambda a, y: [y[1], (a ** 2 - q * a ** 4) / eps ** 2 * y[0]], [0, 2.2 * at], [1, 0],
                    max_step=at / 4000, rtol=1e-10, atol=1e-14)
    y = sol.y[0]; a = sol.t
    yn = y / np.max(np.abs(y[a > at]))
    axs[1].plot(a * np.sqrt(q), np.clip(yn, -1.5, 1.5) * 0.4 + {0.25: 0, 1.0: 1.2, 4.0: 2.4}[q], color=col, lw=0.7,
                label=f"q = {q}")
axs[1].axvline(1, color=GREY, ls="--"); axs[1].text(1.03, 1.75, "a = 1/sqrt(q)\n(the law's slice)", fontsize=8, color=GREY)
axs[1].set_xlabel("size a, in units of 1/sqrt(q)"); axs[1].set_yticks([])
axs[1].set_title("B · Real de Sitter state (K5): no clock below the slice,\nan oscillating clock above it")
axs[1].legend(fontsize=8, loc="lower left")
save(fig, "k3_gravity_clock_and_slice")

# k4 Bell and the shared clock
fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
labels = ["carried\nby each", "shared,\nno signal", "shared +\n0.41 bit avg", "shared +\n1 bit (TB)", "quantum\nmechanics", "no-signal\nmaximum"]
vals = [2, 2, 2 * np.sqrt(2), R["D2"]["S_TB"], 2 * np.sqrt(2), 4]
cols = [GREY, GREY, ORANGE, ORANGE, ORANGE, GREY]
axs[0].bar(range(6), vals, color=cols)
axs[0].axhline(2, color=GREY, ls="--"); axs[0].axhline(2 * np.sqrt(2), color=ORANGE, ls=":")
axs[0].set_xticks(range(6)); axs[0].set_xticklabels(labels, fontsize=7); axs[0].set_ylabel("CHSH")
axs[0].set_title("A · What a shared clock needs to reach 2 sqrt 2 (D1-D3)")
ts = R["D4"]["ts"]; Sc = R["D4"]["S_curve"]
axs[1].plot(ts, Sc, color=ORANGE, label="one quantum clock, coupled turn")
axs[1].axhline(2, color=GREY, ls="--", label="separate turns: stays at 2")
axs[1].axhline(2 * np.sqrt(2), color=GREY, ls=":")
axs[1].set_xlabel("clock reading t"); axs[1].set_ylabel("best CHSH of the pair at that reading")
axs[1].legend(fontsize=8); axs[1].set_title("B · Shared quantum clock (D4): the coupling crosses the line,\nnot the sharing")
save(fig, "k4_shared_clock_bell")

# k5 helical time: the pitch scan
fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
th = np.linspace(0, 4 * np.pi, 800)
for m, col, lab in [(0.5, ORANGE, "spin 1/2"), (1.0, GREY, "spin 1 (m = 1)")]:
    z = np.exp(-1j * m * th)
    axs[0].plot(th / np.pi, z.real, color=col, label=f"{lab}: real part")
    axs[0].plot(th / np.pi, z.imag, color=col, ls="--", label=f"{lab}: i part")
for x, txt in [(0.5, "quarter"), (1, "half: spin 1/2 gets -i"), (2, "full: spin 1/2 gets -1"), (4, "two turns")]:
    axs[0].axvline(x, color=GREY, lw=0.5); axs[0].text(x - 0.12, 1.1, txt, fontsize=7, rotation=90, va="bottom")
axs[0].set_xlabel("twist per loop (units of pi)"); axs[0].set_ylim(-1.2, 3.3)
axs[0].legend(fontsize=7, loc="lower left"); axs[0].set_title("A · Phase a twist gives (T1)")
for n in range(-2, 3):
    axs[1].plot(th / np.pi, 2 * np.pi * n + 0.5 * th, color=ORANGE, lw=0.8)
    axs[1].plot(th / np.pi, 2 * np.pi * n + 1.0 * th, color=GREY, lw=0.5, ls="--")
axs[1].set_ylim(-8, 8); axs[1].set_xlabel("twist per loop (units of pi)"); axs[1].set_ylabel("level x period")
axs[1].set_title("B · Closed helical time (T3): spin 1/2 (orange)\nneeds two full turns to come back")
save(fig, "k5_pitch_scan")

# k6 the law's helix and how much helix nature allows
fig, axs = plt.subplots(1, 2, figsize=(11, 4.4))
tau = np.linspace(0, 3, 600)
for k, col, lab, ls in [(0.6, BLUE, "real k: growth (de Sitter)", "-"), (1j * 3, BLUE, "imaginary k: turn (quantum phase)", "--"),
                        (0.6 + 3j, RED, "complex k: spiral, barred by the sign table", ":")]:
    z = np.exp(k * tau)
    axs[0].plot(z.real, z.imag, color=col, ls=ls, label=lab)
axs[0].set_aspect("equal"); axs[0].legend(fontsize=7, loc="lower left"); axs[0].set_xlim(-7, 7); axs[0].set_ylim(-7, 7)
axs[0].set_title("A · The law's fused shift e^(k tau) (T5)")
H0 = R["T6"]["H0"]
names = ["allowed spin of\nthe cosmos (CMB)", "1 rad per\nHubble time", "electron\nCompton clock", "proton\nCompton clock"]
v = [R["T6"]["Omega_max"], H0, R["T6"]["pitch_e"], R["T6"]["pitch_p"]]
axs[1].bar(range(4), np.log10(v), color=[GREY, BLUE, ORANGE, ORANGE])
for i, x in enumerate(v):
    axs[1].text(i, np.log10(x) + 0.8, f"{x:.0e}", ha="center", fontsize=8)
axs[1].set_xticks(range(4)); axs[1].set_xticklabels(names, fontsize=8); axs[1].set_ylabel("log10 turn rate (rad/s)")
axs[1].set_title("B · Pitch of the helix (T6): space barely turns;\nthe fast helix is internal and set by mass")
save(fig, "k6_law_helix_and_pitch")

# k7 scorecard
items = [("A clock in place of the energy zero keeps the i", "holds", "textbook PW + why-i; tie ours"),
         ("Real numbers without a shared turn", "fails: mirror world", "ours (check K3)"),
         ("Gravity's phase direction sets the i's sign", "holds", "textbook (Banks; Kiefer-Singh)"),
         ("Law's slice 1/sqrt(q) = where the clock switches on", "holds", "HH textbook; reading ours"),
         ("Unknown clock phase makes the odds", "only inside QM", "Zurek; global phase does nothing"),
         ("Shared classical clock reaches 2 sqrt 2", "fails (2)", "Bell"),
         ("Shared clock + signal reaches 2 sqrt 2", "holds, >= 0.41 bit", "Toner-Bacon; Pironio"),
         ("Helical time is the source of the i", "fails: spin-dependent", "ours (T1, T2)"),
         ("Any pitch moves 2 sqrt 2 or the Born odds", "no", "T4"),
         ("A universal pitch singles itself out", "no (data: ~0)", "Saadeh et al. 2016")]
fig, ax = plt.subplots(figsize=(10, 4.6)); ax.axis("off")
tab = ax.table(cellText=[list(r) for r in items], colLabels=["claim", "verdict", "status"], loc="center",
               cellLoc="left", colWidths=[0.52, 0.2, 0.28])
tab.auto_set_font_size(False); tab.set_fontsize(8); tab.scale(1, 1.35)
for (r, c), cell in tab.get_celld().items():
    if r > 0 and c == 1:
        txt = cell.get_text().get_text()
        cell.set_facecolor("#d9ead3" if txt.startswith(("holds", "no")) and "only" not in txt else "#f4cccc" if txt.startswith("fails") else "#fff2cc")
ax.set_title("k7 · Scorecard (DRAFT)")
save(fig, "k7_scorecard")
print("figures written")
