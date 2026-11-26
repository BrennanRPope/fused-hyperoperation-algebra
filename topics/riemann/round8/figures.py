import json, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle
plt.rcParams.update({'font.size': 11, 'axes.spines.top': False, 'axes.spines.right': False})
C = {'blue': '#1f77b4', 'orange': '#ff7f0e', 'green': '#2ca02c', 'red': '#d62728', 'grey': '#7f7f7f', 'purple': '#9467bd'}
F = 'figures/'

# ---------- fig0 overview ----------
rows = [
 ("q-deformed prime gas", "order dial q", "Fisher zeros = a-points of ζ\n(ζ = −q) to first order", "known theory of a-points\n(Levinson): not a lead"),
 ("Calogero–Moser rewind", "rewind (backward heat flow)", "zeros at time t = eigenvalues\nof Z₀ + tL₀, L₀ nilpotent", "in print (Sokal and others):\nnot a lead"),
 ("Twist on a flat surface", "twist (crosscap), whole number", "Klein bottle zeta: 0 of 60\nzeros on the line", "twist breaks the functional\nequation: not a lead"),
 ("Twist on a hyperbolic surface", "twist on the transverse ladder", "glide: ladder sign (−1)ᵏ, cosh weight;\nZ⁻·Z⁺ = cover's factor", "Artin factorisation\nof the cover: not a lead"),
 ("Holes (genus)", "hole count g, crosscaps k", "−χ caps the zeros off the line\n(Otal–Rosas)", "ζ has no χ to cap:\nnot a lead"),
 ("Re-centre the blow-up at s = 1", "fused centre k, Möbius, order q", "k moves local holes, never the\nglobal pole; q = 1 moves it to 1.399", "coordinates change, the pole's\nnature does not: not a lead"),
]
fig, ax = plt.subplots(figsize=(15, 9.2)); ax.axis('off')
heads = ["What we tried", "Which dial", "What it gave", "Verdict"]
cols = [C['orange'], C['blue'], C['green'], C['red']]
W, H = 0.235, 0.115
for j, (hd, c) in enumerate(zip(heads, cols)):
    ax.text(0.02+j*0.245+W/2, 0.965, hd, ha='center', va='center', fontsize=15, color=c, weight='bold')
for i, r in enumerate(rows):
    y = 0.84 - i*0.14
    for j, (txt, c) in enumerate(zip(r, cols)):
        ax.add_patch(FancyBboxPatch((0.02+j*0.245, y-H/2), W, H, boxstyle='round,pad=0.004', fc='white', ec=c, lw=1.6))
        ax.text(0.02+j*0.245+W/2, y, txt, ha='center', va='center', fontsize=10.5)
ax.set_xlim(0, 1); ax.set_ylim(0, 1)
ax.set_title("Round 8: the two open ideas, holes and twists, and re-centring the blow-up", fontsize=16, pad=6)
fig.savefig(F+'fig0-round8-map.png', dpi=110, bbox_inches='tight'); plt.close(fig)

# ---------- fig1 q-gas vs a-points ----------
d = json.load(open('qgas.json'))
fig, axs = plt.subplots(1, 3, figsize=(15, 4.8), sharey=True)
for ax, n in zip(axs, ('1', '2', '3')):
    rows = d[n]['rows']
    for key, col, mk, lab in (('fisher', C['blue'], 'o', 'Fisher zero of the q-gas'), ('apoint', C['orange'], 's', 'a-point of ζ, ζ(β) = −q'), ('linear', C['grey'], '^', 'ρ − q/ζ′(ρ)')):
        pts = [(r['q'], r[key][0]) for r in rows if r[key] is not None and abs(r['q']) <= (0.9 if key != 'apoint' else 0.5)]
        pts.sort()
        ax.plot([p[0] for p in pts], [p[1] for p in pts], mk+'-', color=col, ms=5, lw=1.3, label=lab)
    ax.axhline(0.5, color='k', lw=0.8, ls=':'); ax.axvline(0, color='k', lw=0.6)
    ax.set_title(f"next to ρ_{n} = 1/2 + {d[n]['rho'][1]:.3f}i", fontsize=12)
    ax.set_xlabel("order dial q"); ax.set_ylim(-0.1, 1.5)
axs[0].set_ylabel("Re β")
axs[0].legend(loc='upper right', fontsize=9.5, frameon=False)
fig.suptitle("The q-gas's zeros start out as ζ's a-points: right of the line for q < 0, left for q > 0", fontsize=14)
fig.tight_layout(); fig.savefig(F+'fig1-qgas-apoints.png', dpi=110); plt.close(fig)

# ---------- fig2 Calogero-Moser rewind ----------
d = json.load(open('calogero.json'))
t = np.array(d['rewind_t']); re = np.array(d['traj_re']); im = np.array(d['traj_im'])
fig, axs = plt.subplots(1, 2, figsize=(15, 5.4), gridspec_kw={'width_ratios': [1.6, 1]})
ax = axs[0]
pos = re[0] > 0
for j in np.where(pos)[0]:
    cplx = np.abs(im[:, j]) > 1e-6
    ax.plot(re[:, j], t, color=C['blue'], lw=1)
    ax.plot(re[cplx, j], t[cplx], color=C['red'], lw=2)
ax.set_xlabel("Re z (zero height)"); ax.set_ylabel("rewind time t (backward heat flow, t < 0)")
ax.set_title("Rewinding 30 zeros of ζ (blue: still real; red: the pair has gone complex)", fontsize=12)
ax.set_xlim(5, 110)
ax = axs[1]; ax.axis('off')
txt = ("Zeros under the heat flow obey\n\n"
       "   first order:  ż_k = 2 Σ 1/(z_k − z_j)\n"
       "   second order: z̈_k = −8 Σ 1/(z_k − z_j)³\n\n"
       "which is the rational Calogero–Moser system\n"
       "(coupling g² = −4). Its Lax matrix L is nilpotent,\n"
       "so the zeros at time t are the eigenvalues of\n\n"
       "            Z₀ + t L₀\n\n"
       f"checked to {max(e[1] for e in d['eig_formula_errors']):.0e} at five times.\n"
       f"First pair leaves the line at t = {d['first_complex_t']}:\n"
       f"the closest pair (γ{d['closest_index']}, γ{d['closest_index']+1}), gap {d['closest_gap']:.3f}.")
ax.text(0.0, 0.95, txt, va='top', fontsize=12, family='DejaVu Sans')
fig.tight_layout(); fig.savefig(F+'fig2-calogero-rewind.png', dpi=110); plt.close(fig)

# ---------- fig3 Klein bottle ----------
d = json.load(open('klein.json'))
fig, axs = plt.subplots(1, 2, figsize=(15, 6))
ax = axs[0]
for key, col, mk, lab in (('torus', C['grey'], 'o', 'torus (no twist): all on the line'), ('klein_even', C['red'], 's', 'Klein bottle (twist, even sector)'), ('klein_odd', C['blue'], '^', 'odd sector (torus minus Klein)')):
    z = np.array(d[key]); ax.plot(z[:, 0], z[:, 1], mk, color=col, ms=5, label=lab, alpha=0.85)
ax.axvline(0.5, color='k', lw=0.8, ls=':'); ax.set_xlabel("Re s"); ax.set_ylabel("Im s")
ax.set_title("Zeros of the spectral zeta, square torus vs Klein bottle", fontsize=12)
ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.12), ncol=2, frameon=False, fontsize=10)
ax = axs[1]
for sgn, col, lab in (('1', C['red'], 'twist strength 0 → 1 (Klein)'), ('-1', C['blue'], '0 → −1 (odd sector)')):
    for k, p in enumerate(d['traj'][sgn]):
        p = np.array(p); ax.plot(p[:, 0], p[:, 1], '-', color=col, lw=1.2, label=lab if k == 0 else None)
        ax.plot(p[-1, 0], p[-1, 1], 'o', color=col, ms=4)
    for p in d['traj'][sgn]: ax.plot(p[0][0], p[0][1], 'o', color='k', ms=3)
ax.axvline(0.5, color='k', lw=0.8, ls=':'); ax.set_xlabel("Re s"); ax.set_ylabel("Im s"); ax.set_ylim(0, 41)
ax.set_title("Turning the twist on: every zero leaves the line", fontsize=12)
ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.12), ncol=2, frameon=False, fontsize=10)
fig.tight_layout(); fig.savefig(F+'fig3-klein-twist.png', dpi=110); plt.close(fig)

# ---------- fig4 twist on the ladder ----------
fig, axs = plt.subplots(1, 2, figsize=(15, 5.4))
ax = axs[0]; ax.axis('off'); ax.set_xlim(0, 10); ax.set_ylim(0, 7)
for x0, title, signs, col in ((0.5, "orientation-preserving orbit", [+1]*5, C['blue']), (5.5, "glide reflection (twist)", [+1, -1, +1, -1, +1], C['red'])):
    ax.text(x0+1.9, 6.6, title, ha='center', fontsize=12, color=col, weight='bold')
    for k, sg in enumerate(signs):
        y = 0.8 + k*1.1
        ax.plot([x0, x0+3.8], [y, y], color=col, lw=2)
        ax.text(x0+1.9, y+0.2, f"1 {'−' if sg > 0 else '+'} e^-(s+{k})ℓ", ha='center', fontsize=11)
    ax.plot([x0, x0], [0.8, 0.8+4*1.1], color='k', lw=1); ax.plot([x0+3.8, x0+3.8], [0.8, 0.8+4*1.1], color='k', lw=1)
ax.text(5, 0.05, "rung by rung: (1 − x)(1 + x) = 1 − x²  ⇒  Z⁻_ℓ · Z⁺_ℓ = Z_2ℓ  (the cover sees one orbit of length 2ℓ)", ha='center', fontsize=11)
ax = axs[1]
l = np.linspace(0.1, 6, 300)
ax.plot(l, 1/(2*np.sinh(l/2)), color=C['blue'], lw=2, label="1/(2 sinh(ℓ/2)): Bose ladder, Σ e^-kℓ")
ax.plot(l, 1/(2*np.cosh(l/2)), color=C['red'], lw=2, label="1/(2 cosh(ℓ/2)): alternating ladder, Σ (−1)^k e^-kℓ")
ax.set_xlabel("orbit length ℓ"); ax.set_ylabel("weight in the trace formula"); ax.set_ylim(0, 2.5)
ax.legend(frameon=False, fontsize=10.5); ax.set_title("The twist turns the transverse ladder from Bose to alternating", fontsize=12)
fig.tight_layout(); fig.savefig(F+'fig4-twist-ladder.png', dpi=110); plt.close(fig)

# ---------- fig5 holes and crosscaps cap the off-line zeros ----------
fig, ax = plt.subplots(figsize=(13, 5.2))
g = np.arange(2, 7); k = np.arange(3, 11)
ax.bar(-(2-2*g)-0.18, -(2-2*g), width=0.34, color=C['blue'], label="orientable, genus g: at most −χ = 2g − 2 eigenvalues below 1/4 (Otal–Rosas)")
ax.bar(-(2-k)+0.18, -2*(2-k), width=0.34, color=C['red'], alpha=0.8, label="k crosscaps: at most −2χ via the double cover (−χ if Ballmann–Matthiesen–Mondal applies)")
for gi in g: ax.text(-(2-2*gi)-0.18, -(2-2*gi)+0.3, f"g={gi}", ha='center', fontsize=9, color=C['blue'])
for ki in k: ax.text(-(2-ki)+0.18, -2*(2-ki)+0.3, f"k={ki}", ha='center', fontsize=9, color=C['red'])
ax.set_xlabel("−χ (Gauss–Bonnet: area = −2πχ)"); ax.set_ylabel("most zeros of Selberg's zeta off the line\n(real, in (1/2, 1], counting s = 1)")
ax.set_title("Holes and twists are whole numbers, and −χ caps how many zeros can leave the line; ζ has no χ", fontsize=12)
ax.legend(loc='upper left', frameon=False, fontsize=9.5); ax.set_ylim(0, 19)
fig.tight_layout(); fig.savefig(F+'fig5-genus-cap.png', dpi=110); plt.close(fig)

# ---------- fig6 quantum gravity: GUE vs GOE ----------
d = json.load(open('topology.json'))
sg = np.array(d['sgrid']); sp = np.array(d['zeta_spacings'])
fig, ax = plt.subplots(figsize=(12, 5.2))
ax.hist(sp, bins=50, range=(0, 3.5), density=True, color='#cccccc', label=f"ζ zeros at heights 10000–11600 ({len(sp)} spacings)")
ax.plot(sg, d['p_gue'], color=C['blue'], lw=2.2, label=f"GUE: JT gravity summed over orientable surfaces (KS {d['ks_gue']:.3f})")
ax.plot(sg, d['p_goe'], color=C['red'], lw=2.2, label=f"GOE: crosscaps included, Stanford–Witten (KS {d['ks_goe']:.3f})")
ax.set_xlim(0, 3.5); ax.set_xlabel("nearest-neighbour spacing (unfolded)"); ax.set_ylabel("density")
ax.set_title("On the quantum-gravity side, ζ's zeros look like a sum with no twists", fontsize=12)
ax.legend(frameon=False, fontsize=10)
fig.tight_layout(); fig.savefig(F+'fig6-gue-not-goe.png', dpi=110); plt.close(fig)

# ---------- fig7 re-centring ----------
d = json.load(open('recentre.json'))
fig, axs = plt.subplots(1, 3, figsize=(16, 5.4))
ax = axs[0]
th = np.linspace(-0.12, 0.12, 400); ax.plot(-1+np.cos(th), np.sin(th), color=C['grey'], lw=1.2, label="image of Re s = 1/2: the circle |w + 1| = 1")
w = np.array(d['mobius_w']); ax.plot(w[:, 0], w[:, 1], 'o', color=C['blue'], ms=3.5, label="first 100 zeros, w = 1/(s − 1)")
ax.plot(w[:, 0], -w[:, 1], 'o', color=C['blue'], ms=3.5, alpha=0.35)
ax.plot(0, 0, '*', color=C['red'], ms=12, label="w = 0 is s = ∞; the pole s = 1 is now w = ∞")
ax.set_xlim(-0.008, 0.002); ax.set_ylim(-0.085, 0.085)
ax.set_xlabel("Re w (zoomed near w = 0)"); ax.set_ylabel("Im w")
ax.set_title("Möbius re-centring: the line becomes a circle", fontsize=12)
ax.legend(loc='upper left', frameon=False, fontsize=9)
ax = axs[1]
ks = sorted(float(k) for k in d['Ck']); cs = [d['Ck'][str(k) if str(k) in d['Ck'] else k] for k in ks]
ax.plot(ks, cs, 'o-', color=C['green'], lw=2)
ax.set_xlabel("fused centre dial k"); ax.set_ylabel("residue of ζ_k at s = 1")
ax.set_title("Fusion keeps the pole at s = 1; only its size moves", fontsize=12)
ax.annotate("k = −1: 6/π² (Fermi)", xy=(-1, cs[0]), xytext=(-0.9, 1.9), arrowprops=dict(arrowstyle='->'), fontsize=10)
ax.annotate("k → 2: the local hole of p = 2\nlands on s = 1; order 1 → 3/2", xy=(1.9, cs[-1]), xytext=(-0.6, 2.6), arrowprops=dict(arrowstyle='->'), fontsize=10)
ax = axs[2]
b = np.array(d['beta']); ax.plot(b, d['zeta_real'], color=C['blue'], lw=2, label="ζ(β): Bose, q = 0")
Zb = np.array(d['Z_boltz']); Zb[np.abs(Zb) > 12] = np.nan
ax.plot(b, Zb, color=C['red'], lw=2, label="1/(1 − P(β)): Boltzmann, q = 1")
ax.axvline(1, color=C['blue'], ls=':', lw=1); ax.axvline(d['wall'], color=C['red'], ls=':', lw=1)
ax.set_ylim(-6, 12); ax.set_xlabel("β (real)"); ax.set_title(f"The order dial does move it: to β = {d['wall']:.4f}", fontsize=12)
ax.legend(frameon=False, fontsize=10, loc='upper right')
fig.tight_layout(); fig.savefig(F+'fig7-recentre.png', dpi=110); plt.close(fig)

# ---------- fig8 subtracting vs multiplying ----------
fig, ax = plt.subplots(figsize=(12, 5))
st = np.array(d['stieltjes_zeros'])
import mpmath as mp
rz = [float(mp.im(mp.zetazero(n))) for n in range(1, 11)]
ax.plot([0.5]*10, rz, 'o', color=C['grey'], ms=8, label="zeros of ζ, (s − 1)ζ and ξ: the same, on the line")
ax.plot(st[:, 0], st[:, 1], 's', color=C['red'], ms=7, label="zeros of ζ(s) − 1/(s − 1): pole subtracted, zeros off the line")
for (x, y), y0 in zip(st, rz): ax.plot([0.5, x], [y0, y], color=C['red'], lw=0.8)
ax.axvline(0.5, color='k', ls=':', lw=0.8); ax.set_xlim(0.46, 0.54); ax.set_xlabel("Re s"); ax.set_ylabel("Im s")
ax.set_title("Multiplying the pole away keeps the zeros; subtracting it moves them", fontsize=12)
ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.13), frameon=False, ncol=1, fontsize=10)
fig.tight_layout(); fig.savefig(F+'fig8-remove-pole.png', dpi=110); plt.close(fig)
print('figures done')
