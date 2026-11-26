import numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lib import *
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
S = lambda f: (plt.tight_layout(), plt.savefig("figures/"+f, dpi=130), plt.close())
out = {}

# Fig 1: why field matters and what pushes back
B = np.linspace(1, 30, 300)
fig, ax = plt.subplots(figsize=(8.5, 4.2))
ax.plot(B, (B/12.2)**4, lw=2.5, label="fusion power per volume ∝ B⁴ (same plasma pressure ratio)")
ax.plot(B, (B/12.2)**2, lw=2.5, label="magnetic pressure B²/2μ₀ ∝ B² (the stress the coil must hold)")
ax.plot(B, Ic(B, 20)/Ic(12.2, 20), lw=2.5, label="current one REBCO tape carries at 20 K (fit)")
ax.axvline(12.2, color="gray", ls=":"); ax.text(12.4, 4.6, "SPARC plasma\ncentre ≈12.2 T", color="gray")
ax.axvline(20, color="gray", ls="--"); ax.text(20.3, 3.6, "TFMC peak on\ncoil 20 T (2021)", color="gray")
ax.set_ylim(0, 8); ax.set_xlabel("field B (tesla)"); ax.set_ylabel("relative to the value at 12.2 T")
ax.set_title("Fig 1. Field buys power as B⁴, but stress grows as B² and the tape carries less current")
ax.legend(loc="upper left", fontsize=8.5); S("fig1_why_field_and_what_fights_back.png")

# Fig 2: critical surface
fig, ax = plt.subplots(figsize=(7.5, 4.2))
for T, c in [(4.2, "C0"), (20, "C1"), (30, "C2"), (50, "C3"), (77, "C4")]:
    ax.plot(B, Ic(B, T), lw=2, color=c, label=f"{T:g} K")
ax.scatter([20], [Ic(20, 20)], s=60, color="k", zorder=5); ax.annotate("TFMC-like point\n20 T, 20 K", (20, Ic(20, 20)), (22, 1100), arrowprops=dict(arrowstyle="->"))
ax.set_yscale("log"); ax.set_ylim(10, 5000); ax.set_xlabel("field perpendicular to tape (T)"); ax.set_ylabel("critical current of one tape (A)")
ax.set_title("Fig 2. REBCO critical surface: colder buys current, field takes it away\n(Succi–Ballarino fit form, Faraday Factory constants; a fit to measured data)", fontsize=9.5)
ax.legend(title="temperature"); S("fig2_rebco_critical_surface.png")

# Fig 3: derating composition and the corner test
fig, axs = plt.subplots(1, 2, figsize=(10.5, 4.3))
a = axs[0]; B1, B2, T1, T2 = 10, 20, 20, 30
pts = {(0,0):(B1,T1),(1,0):(B2,T1),(0,1):(B1,T2),(1,1):(B2,T2)}
for (i,j),(b,t) in pts.items():
    sgn = "+" if (i+j)%2==0 else "−"
    a.scatter(i, j, s=2600, color="C0" if sgn=="+" else "C3", alpha=.25)
    a.text(i, j, f"sign {sgn}\n{b} T, {t} K\nIc = {Ic(b,t):.0f} A", ha="center", va="center", fontsize=8.5)
a.plot([0,1,1,0,0],[0,0,1,1,0], color="gray", lw=1)
cs = corner_sum_log(Ic, B1, B2, T1, T2)
a.set_xlim(-.6, 1.6); a.set_ylim(-.6, 1.7); a.axis("off")
a.set_title(f"Corner sum of ln Ic (signs ±) = {cs:+.3f}\n0 would mean field and temperature deratings simply multiply", fontsize=9)
a = axs[1]
Bg, Tg = np.meshgrid(np.linspace(5, 25, 80), np.linspace(10, 50, 80))
dB = 1 - Ic(Bg, 20)/Ic(5, 20); dT = 1 - Ic(5, Tg)/Ic(5, 20)
gap = (1 - Ic(Bg, Tg)/Ic(5, 20)) - circ(dB, dT, -1)
im = a.contourf(Bg, Tg, 100*gap, 20, cmap="RdBu_r", vmin=-abs(100*gap).max(), vmax=abs(100*gap).max()); plt.colorbar(im, ax=a, label="true derating − ∘₋₁ prediction (percentage points)")
a.set_xlabel("field (T)"); a.set_ylabel("temperature (K)"); a.set_title("Error of composing deratings with ∘₋₁ (reference 5 T, 20 K)", fontsize=9)
fig.suptitle("Fig 3. Deratings compose by the union rule only when the critical surface separates; the corner sum measures the failure", fontsize=10)
S("fig3_derating_union_rule_and_corner_test.png"); out["max_gap_pp"] = float(abs(100*gap).max())

# Fig 4: weakest link and detection
fig, axs = plt.subplots(1, 2, figsize=(10, 4))
Lk = np.linspace(0, 20, 200)
for lam in [0.05, 0.2, 0.5]:
    axs[0].plot(Lk, 1-np.exp(-lam*Lk), lw=2, label=f"{lam} weak spots per km")
    ns = np.arange(1, 21, 4); axs[0].scatter(ns, [fold([lam*n/50]*50, -1) for n in ns], s=14, color="k")
axs[0].set_xlabel("km of tape in the coil"); axs[0].set_ylabel("chance the coil has ≥1 weak spot")
axs[0].set_title("∘₋₁ folded over segments (dots) = weakest-link law (lines)", fontsize=9); axs[0].legend(fontsize=8)
n = np.arange(1, 9)
for p in [0.5, 0.8, 0.95]:
    axs[1].plot(n, fold([p]*8, -1) if False else [fold([p]*m, -1) for m in n], "o-", label=f"each sensor catches it with p={p}")
axs[1].set_xlabel("independent detectors watching the same spot"); axs[1].set_ylabel("chance at least one catches it")
axs[1].set_title("Union rule on detectors = multiply the miss chances", fontsize=9); axs[1].legend(fontsize=8)
fig.suptitle("Fig 4. The union rule ∘₋₁ restates series reliability; it adds no new handle on quench", fontsize=10); S("fig4_union_rule_weakest_link_and_detection.png")

# Fig 5: strain composition
fig, ax = plt.subplots(figsize=(7.5, 4))
e = np.linspace(0, 0.006, 100)
ax.plot(100*e, 100*e*e, lw=2.5, label="cross term ε₁ε₂ that ∘₁ adds (ε₁ = ε₂ = ε)")
ax.axhspan(0, 0.0, color="w")
ax.fill_between(100*e, 0, 0*e+0.2, color="C1", alpha=.12, label="0.2–0.4% critical strains used in coil design")
ax.fill_between(100*e, 0.2, 0.4, color="C1", alpha=.25)
ax.set_yscale("symlog", linthresh=1e-4); ax.set_xlabel("each strain (%)"); ax.set_ylabel("strain (%)")
ax.set_title("Fig 5. Stacked strains compose by ∘₁ exactly, but the correction is\nhundreds of times smaller than the design limits", fontsize=10); ax.legend(fontsize=8.5); S("fig5_strain_composition_o1.png")

# Fig 6: semiring dial = KS / Lp
rng = np.random.default_rng(3); x = np.linspace(0, 1, 400)
g = 300 + 120*np.sin(9*x) + 90*np.exp(-((x-.72)/.03)**2) + rng.normal(0, 5, x.size)
fig, axs = plt.subplots(1, 2, figsize=(10.5, 4))
rhos = np.logspace(-3, 0, 60)
axs[0].semilogx(rhos, [KS(g, r) for r in rhos], lw=2.5, label="KS aggregate (log-sum-exp at dial 1/ρ)")
axs[0].axhline(g.max(), color="k", ls="--", label="true peak stress")
axs[0].fill_between(rhos, g.max(), g.max()+np.log(g.size)/rhos, color="C0", alpha=.12, label="guaranteed band  max + ln N / ρ")
axs[0].set_ylim(400, 1500); axs[0].set_xlabel("dial ρ (1/MPa)"); axs[0].set_ylabel("aggregate (MPa)"); axs[0].legend(fontsize=8)
axs[0].set_title("Turning the dial slides from 'sum' to 'max'", fontsize=9)
axs[1].plot(x, g, color="gray", lw=1, label="stress along a coil (toy)")
ax2 = axs[1].twinx()
for r, c in [(0.01, "C2"), (0.05, "C1"), (0.3, "C3")]:
    w = np.exp(r*(g-g.max())); ax2.plot(x, w/w.sum(), color=c, label=f"gradient weights, ρ={r}")
axs[1].set_xlabel("position along coil"); axs[1].set_ylabel("stress (MPa)"); ax2.set_ylabel("share of the gradient"); ax2.legend(fontsize=7.5, loc="upper left")
axs[1].set_title("The weights are a softmax: high ρ listens only to the worst spot", fontsize=9)
fig.suptitle("Fig 6. The semiring dial is the Kreisselmeier–Steinhauser aggregate (and Lp norms in log units): known and already used", fontsize=10)
S("fig6_semiring_dial_is_KS_aggregation.png")

# Fig 7: tape twist on a closed coil
def panel(X, ax3, axs, title):
    F = frames(X); U, hol = bishop(F['T'], F['N'][0]); A = spherical_area(F['T'])
    ax3.plot(*X.T, lw=2); 
    for i in range(0, len(X), 160):
        p = X[i]; ax3.plot(*np.stack([p-0.35*U[i], p+0.35*U[i]]).T, color="C1", lw=1.2)
    ax3.set_title(title, fontsize=9); ax3.set_box_aspect((np.ptp(X[:,0])+.1, np.ptp(X[:,1])+.1, np.ptp(X[:,2])+.1))
    u, v = np.mgrid[0:2*np.pi:40j, 0:np.pi:20j]
    axs.plot_wireframe(np.cos(u)*np.sin(v), np.sin(u)*np.sin(v), np.cos(v), color="lightgray", lw=.4)
    axs.plot(*F['T'].T, color="C3", lw=2.5); axs.set_box_aspect((1,1,1)); axs.set_axis_off()
    axs.set_title(f"tangent path on the sphere: encloses {abs(A):.3f} sr\nflat tape comes back turned {np.degrees(abs(wrap(hol))):.2f}°", fontsize=9)
    return F, hol, A
fig = plt.figure(figsize=(11, 8))
r1 = panel(modular_coil(), fig.add_subplot(2,2,1,projection="3d"), fig.add_subplot(2,2,2,projection="3d"), "non-planar modular coil\n(orange = tape face direction, carried along with no twist)")
r2 = panel(planar_coil(), fig.add_subplot(2,2,3,projection="3d"), fig.add_subplot(2,2,4,projection="3d"), "planar coil")
fig.suptitle("Fig 7. Swap rule at d = +1 on a coil: the untwisted tape's turn = 2π − area enclosed by the tangent path\n(Fuller 1971 / Berry phase; planar coils enclose exactly 2π and need no net twist)", fontsize=10)
S("fig7_tape_twist_is_enclosed_area.png")

# Fig 8: twist budget
rows = []
for name, X in [("gentle modular coil", modular_coil()), ("coil with a dent", dented_coil())]:
    F = frames(X); U, hol = bishop(F['T'], F['N'][0]); d = abs(wrap(hol)); L = F['L']
    sinang = np.abs(np.einsum('ij,ij->i', np.cross(U, F['N']), F['T']))
    for w in [0.004, 0.012]:
        rows.append((f"{name}\n{w*1000:.0f} mm tape", (d/L)**2*w**2/12, np.max(F['tau']**2*w**2/12), np.max(w*F['kappa']*sinang/2)))
fig, ax = plt.subplots(figsize=(10, 4.8))
xs = np.arange(len(rows)); wd = .27
for j, (lab, c) in enumerate([("floor set by closing the loop (swap rule)", "C2"), ("twist if tape follows the curve's own frame", "C0"), ("hard-way bend if tape is held untwisted", "C3")]):
    ax.bar(xs+(j-1)*wd, [r[j+1]*100 for r in rows], wd, color=c, label=lab)
ax.axhline(0.2, color="k", ls="--"); ax.text(xs[-1]+.45, .22, "0.2% design limit", ha="right")
ax.set_yscale("log"); ax.set_xticks(xs); ax.set_xticklabels([r[0] for r in rows], fontsize=8.5); ax.set_ylabel("peak strain (%)")
ax.set_title("Fig 8. The loop-closing floor is negligible; the real fight is local twist versus hard-way bending", fontsize=10); ax.legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=3)
S("fig8_twist_budget.png"); out["twist_rows"] = [(r[0].replace("\n"," "), r[1], r[2], r[3]) for r in rows]

# Fig 9: scorecard
limits = ["Critical current vs field, temperature", "Lorentz / hoop stress", "Strain on tape", "Quench detection & protection", "No-insulation coil behaviour", "Coil & stellarator path optimisation", "Tape twist on 3-D coils", "Radiation damage", "Cost"]
tools = ["union rule ∘₋₁", "corner sum", "∘₁ composition", "semiring dial", "swap rule d=+1"]
M = np.full((len(limits), len(tools)), 0)
# 0 nothing, 1 restates known, 2 small practical use
M[0,0]=1; M[0,1]=2; M[2,2]=1; M[3,0]=1; M[5,3]=1; M[1,3]=1; M[6,4]=1; M[6,3]=1; M[2,1]=1
cmap = matplotlib.colors.ListedColormap(["#eeeeee", "#9ecae1", "#fdae6b"])
fig, ax = plt.subplots(figsize=(9, 5.2)); ax.imshow(M, cmap=cmap, vmin=0, vmax=2, aspect="auto")
ax.set_xticks(range(len(tools))); ax.set_xticklabels(tools, rotation=20); ax.set_yticks(range(len(limits))); ax.set_yticklabels(limits)
lab = {0: "nothing", 1: "restates known", 2: "handy check (known method)"}
for i in range(len(limits)):
    for j in range(len(tools)): ax.text(j, i, lab[M[i,j]], ha="center", va="center", fontsize=7.5)
ax.set_title("Fig 9. Scorecard: what each tool of the theory does for each limit", fontsize=10); S("fig9_scorecard.png")
import json; json.dump(out, open("figures_output.json", "w"), indent=1, default=float); print(json.dumps(out, indent=1, default=float))
