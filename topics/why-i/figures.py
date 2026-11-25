"""Figures for topics/why-i. Run: python3 figures.py  -> figures/*.png"""
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, FancyArrowPatch, Ellipse, Circle, Arc, FancyBboxPatch
OUT = 'figures/'
BLUE, ORANGE, GREEN, PURPLE, RED, GREY = '#1f77b4', '#ff7f0e', '#2ca02c', '#9467bd', '#d62728', '#777777'
plt.rcParams.update({'font.size': 11, 'figure.dpi': 130})
def arrow(ax, p, q, col='k', lw=2, ms=16, **kw):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle='-|>', mutation_scale=ms, color=col, lw=lw, **kw))

# ---------- fig1: the twist ----------
fig, axs = plt.subplots(1, 2, figsize=(11, 4.8))
for ax, (u, v, col, lab) in zip(axs, [((2, .6), (.7, 1.8), BLUE, 'shift-stretch A, then B'), ((.7, 1.8), (2, .6), ORANGE, 'B, then A')]):
    P = np.array([(0, 0), u, (u[0]+v[0], u[1]+v[1]), v])
    ax.add_patch(Polygon(P, closed=True, fc=col, alpha=.25, ec=col, lw=2))
    arrow(ax, (0, 0), u, col); arrow(ax, (0, 0), v, GREY)
    th = np.linspace(0, 2*np.pi*0.8, 50)
    cx, cy = (u[0]+v[0])/2, (u[1]+v[1])/2
    sgn = 1 if col == BLUE else -1
    ax.plot(cx + .35*np.cos(sgn*th), cy + .35*np.sin(sgn*th), color=col, lw=2)
    ax.annotate('', xy=(cx + .35*np.cos(sgn*th[-1]), cy + .35*np.sin(sgn*th[-1])),
                xytext=(cx + .35*np.cos(sgn*th[-3]), cy + .35*np.sin(sgn*th[-3])), arrowprops=dict(arrowstyle='-|>', color=col, lw=2))
    ax.text(cx, cy - .75, ('+ area' if sgn > 0 else '− area') + '\n(counter-clockwise)' * (sgn > 0) + '\n(clockwise)' * (sgn < 0),
            ha='center', fontsize=11, color=col, weight='bold')
    ax.set_xlim(-.4, 3.2); ax.set_ylim(-.6, 2.8); ax.set_aspect('equal'); ax.set_title(lab)
    ax.set_xlabel('stretch part'); ax.set_ylabel('shift part')
fig.suptitle('1. The twist: two stretch-and-shift maps in either order differ by an oriented area a·d − b·c.\n'
             'Swap the order and the area keeps its size but flips its sign. That sign is the "twist".', fontsize=12)
fig.tight_layout(); fig.savefig(OUT + 'fig1-the-twist.png', bbox_inches='tight'); plt.close(fig)

# ---------- fig2: area alone does not pick an i ----------
fig, axs = plt.subplots(1, 2, figsize=(11.5, 5))
ax = axs[0]
for ang, ratio, col in [(0, 1, GREEN), (20, 2.2, BLUE), (70, 1.6, ORANGE), (-35, 3, PURPLE), (110, 2.6, RED)]:
    w, h = 2*np.sqrt(ratio), 2/np.sqrt(ratio)
    ax.add_patch(Ellipse((0, 0), w, h, angle=ang, fc='none', ec=col, lw=2.2 if col == GREEN else 1.6))
ax.set_xlim(-2.2, 2.2); ax.set_ylim(-2.2, 2.2); ax.set_aspect('equal'); ax.axhline(0, c='#ddd'); ax.axvline(0, c='#ddd')
ax.set_title('Twist only: every one of these "unit circles"\nencloses the same area with the same orientation')
ax.text(0, -1.9, 'each ellipse is a different candidate i (a different J with J² = −1)', ha='center', fontsize=9.5)
ax = axs[1]
ax.add_patch(Circle((0, 0), 1, fc='none', ec=GREEN, lw=2.5))
u = np.array([1., 0]); Ju = np.array([0, 1.]); JJu = -u
arrow(ax, (0, 0), u, BLUE); arrow(ax, (0, 0), Ju, ORANGE); arrow(ax, (0, 0), JJu, RED)
ax.text(1.05, -.15, 'u', color=BLUE, fontsize=13, weight='bold'); ax.text(.08, 1.05, 'J u', color=ORANGE, fontsize=13, weight='bold')
ax.text(-1.45, -.2, 'J J u = −u', color=RED, fontsize=13, weight='bold')
ax.add_patch(Arc((0, 0), .6, .6, theta1=0, theta2=90, color=ORANGE, lw=2)); ax.add_patch(Arc((0, 0), .9, .9, theta1=90, theta2=180, color=RED, lw=2))
ax.set_xlim(-1.8, 1.8); ax.set_ylim(-1.5, 1.6); ax.set_aspect('equal'); ax.set_axis_off()
ax.set_title('Twist + a length: the length picks the round circle,\nand J = "quarter turn the positive way". Twice is −1.')
fig.suptitle('2. Known (compatible triples): an oriented area alone allows a whole half-plane of i\'s. Add a length and exactly one is left.', fontsize=12)
fig.tight_layout(); fig.savefig(OUT + 'fig2-area-alone-vs-area-plus-length.png', bbox_inches='tight'); plt.close(fig)

# ---------- fig3: flip the twist, flip i ----------
fig, axs = plt.subplots(1, 2, figsize=(10.5, 4.6))
for ax, sgn, title in [(axs[0], 1, 'Orientation as drawn:\nquarter turn counter-clockwise = i'), (axs[1], -1, 'Mirror (orientation flipped):\nquarter turn clockwise = −i')]:
    ax.add_patch(Circle((0, 0), 1, fc='none', ec=GREY, lw=1))
    arrow(ax, (0, 0), (1, 0), BLUE); arrow(ax, (0, 0), (0, sgn), ORANGE)
    ax.add_patch(Arc((0, 0), .7, .7, theta1=0 if sgn > 0 else -90, theta2=90 if sgn > 0 else 0, color=ORANGE, lw=2.5))
    ax.text(1.05, .05, 'u', color=BLUE, fontsize=13); ax.text(.1, sgn*1.08, 'J u' if sgn > 0 else '(−J) u', color=ORANGE, fontsize=13)
    ax.set_xlim(-1.4, 1.6); ax.set_ylim(-1.4, 1.4); ax.set_aspect('equal'); ax.set_axis_off(); ax.set_title(title)
fig.suptitle('3. The sign of the twist is the sign of i. Complex conjugation, mirror reflection and time reversal (Wigner)\n'
             'are the same move: i ↔ −i. Physics never decides which of the two square roots of −1 is "i".', fontsize=11.5)
fig.tight_layout(); fig.savefig(OUT + 'fig3-flip-the-twist-flip-i.png', bbox_inches='tight'); plt.close(fig)

# ---------- fig4: the order dial's group is the upper half-plane ----------
fig, ax = plt.subplots(figsize=(11.5, 6.6))
ax.axhline(0, color='k', lw=2)
ax.fill_between([-6, 6], 0, 4.2, color='#f3f6fb')
for cc in np.linspace(-6, 6, 13):
    for rr in [1, 2.5]:
        t = np.linspace(0, np.pi, 100); ax.plot(cc + rr*np.cos(t), rr*np.sin(t), color='#dde3ee', lw=.8)
dials = [(0, GREEN, 'k = 0: plain addition\n(horizontal: a horocycle)'), (1, BLUE, 'k = 1: mult-add'), (.4, PURPLE, 'k = 0.4'),
         (-1, ORANGE, 'k = −1: probability join'), (-.4, RED, 'k = −0.4')]
for kk, col, lab in dials:
    if kk == 0:
        ax.plot([-6, 6], [1, 1], color=col, lw=2.5); ax.text(3.6, 1.08, lab, color=col, fontsize=9.5)
        continue
    c0 = -1/kk; s = np.linspace(-6, 6, 600)
    lam = 1 + kk*s; m = (lam >= 0) & (lam <= 4.2)
    ax.plot(s[m], lam[m], color=col, lw=2.5)
    ax.plot([c0], [0], 'o', color=col, ms=9)
    ax.text(c0, -.32, f'centre −1/k = {c0:g}', color=col, ha='center', fontsize=9)
    sx = 1.6 if kk > 0 else -1.6; ax.text(sx, 1 + kk*sx + .12, lab, color=col, fontsize=9.5, ha='left' if kk > 0 else 'right')
ax.plot([0, 0], [0, 4.2], color='k', lw=2.5, ls='--'); ax.text(.08, 3.7, 'pure multiplication\n(k → ∞: a true geodesic)', fontsize=9.5)
ax.plot([0], [1], 'o', color='k', ms=12, zorder=5); ax.text(.15, .45, 'i = the identity map\n(do nothing)', fontsize=11, weight='bold', bbox=dict(fc='white', ec='none', alpha=.8))
ax.annotate('edge: stretch → 0', xy=(-5.5, 0), xytext=(-5.8, .5), fontsize=9.5, arrowprops=dict(arrowstyle='->'))
ax.set_xlim(-6, 6); ax.set_ylim(-.6, 4.2); ax.set_aspect('equal')
ax.set_xlabel('shift  s'); ax.set_ylabel('stretch  λ  (= 1 + k·a)')
ax.set_title('4. The order dial\'s home, drawn: every stretch-and-shift map x ↦ λx + s is a point s + iλ of the upper half-plane.\n'
             'Each dial k is a straight line through i (doing nothing) that ends at its centre −1/k on the edge. Faint arcs: the half-plane\'s own straight lines.',
             fontsize=11)
fig.tight_layout(); fig.savefig(OUT + 'fig4-dial-rays-in-the-half-plane.png', bbox_inches='tight'); plt.close(fig)

# ---------- fig5: why a phase ----------
fig, axs = plt.subplots(1, 2, figsize=(11.5, 4.8))
z = np.linspace(0, 6, 400)
ax = axs[0]
ax.plot(z, np.exp(.5*z), color=RED, lw=2.5, label='real: e^(λz) (grows)')
ax.plot(z, np.exp(-.5*z), color=ORANGE, lw=2.5, label='real: e^(−λz) (fades)')
ax.plot(z, np.abs(np.exp(.5j*z)), color=GREEN, lw=3, label='turn: |e^(iλz)| = 1 always')
ax.set_xlabel('total twist z (adds up as you compose)'); ax.set_ylabel('size'); ax.set_ylim(0, 6); ax.legend(loc='upper left')
ax.set_title('Twists add. If each twist multiplies by a number,\nreal numbers change the size; only a turn keeps it')
ax = axs[1]
t = np.linspace(0, 2*np.pi, 200); ax.plot(np.cos(t), np.sin(t), color=GREY)
for zz, col in [(0, 'k'), (.9, BLUE), (1.8, PURPLE), (2.7, ORANGE)]:
    ax.plot([np.cos(zz)], [np.sin(zz)], 'o', color=col, ms=10); ax.text(1.12*np.cos(zz), 1.12*np.sin(zz), f'z = {zz:g}', color=col, ha='center')
ax.text(0, 0, 'composing twists\n= walking round\nthe circle', ha='center', va='center')
ax.set_aspect('equal'); ax.set_axis_off(); ax.set_xlim(-1.5, 1.5); ax.set_ylim(-1.4, 1.4)
ax.set_title('With real numbers alone, a size-keeping twist must be a\nrotation of a plane; its quarter turn J has J² = −1')
fig.suptitle('5. Known (unitarity; Stone–von Neumann): the twist becomes "i" the moment you demand that composing never changes total size.', fontsize=12)
fig.tight_layout(); fig.savefig(OUT + 'fig5-why-a-phase.png', bbox_inches='tight'); plt.close(fig)

# ---------- fig6: commutator eigenvalues ----------
rng = np.random.default_rng(3)
fig, ax = plt.subplots(figsize=(7.5, 6))
for n, col in [(2, BLUE), (3, GREEN), (5, ORANGE), (8, PURPLE)]:
    ev = []
    for _ in range(60):
        M1, M2 = rng.normal(size=(n, n)), rng.normal(size=(n, n)); A1, A2 = M1 + M1.T, M2 + M2.T
        ev += list(np.linalg.eigvals(A1 @ A2 - A2 @ A1))
    ev = np.array(ev); ax.plot(ev.real, ev.imag, 'o', ms=4, alpha=.6, color=col, label=f'{n}×{n} real symmetric pairs')
ax.axhline(0, c='#ccc'); ax.axvline(0, c='k', lw=1)
ax.set_xlim(-75, 75); ax.set_xlabel('real part'); ax.set_ylabel('imaginary part'); ax.legend(loc='upper right', fontsize=9)
ax.set_title('6. Stueckelberg\'s point, as a picture: take two real "observables" (symmetric matrices).\n'
             'Their order cost AB − BA has eigenvalues only on the imaginary axis, ± pairs.\n'
             'Even an all-real theory meets i the moment order matters.', fontsize=10.5)
fig.tight_layout(); fig.savefig(OUT + 'fig6-twist-eigenvalues-are-imaginary.png', bbox_inches='tight'); plt.close(fig)

# ---------- fig7: the wall flips the sign of i ----------
N, Lb = 4096, 400.; xs = np.linspace(-Lb/2, Lb/2, N, endpoint=False); dx = xs[1]-xs[0]; xi = 2*np.pi*np.fft.fftfreq(N, d=dx)
seed = np.exp(-xs**2/30.)*np.exp(1j*1.5*xs); F = np.fft.fft(seed); psi = np.fft.ifft(np.where(xi > 0, F, 0))
def U(lam, s, p):
    return (np.interp((xs-s)/lam, xs, p.real, left=0, right=0) + 1j*np.interp((xs-s)/lam, xs, p.imag, left=0, right=0))/np.sqrt(abs(lam))
cases = [('starting wave (only + frequencies)', psi, BLUE), ('after mult-add element a = −0.5, k = 1\n(1 + ka = 0.5 > 0, before the wall)', U(.5, -.5, psi), GREEN),
         ('after a = −2, k = 1\n(1 + ka = −1 < 0, past the wall)', U(-1, -2, psi), RED), ('complex conjugate of the start\n(i → −i)', np.conj(psi), PURPLE)]
fig, axs = plt.subplots(1, 4, figsize=(14, 3.8), sharey=True)
order = np.argsort(xi)
for ax, (lab, p, col) in zip(axs, cases):
    P = np.abs(np.fft.fft(p))**2; P /= P.max()
    ax.fill_between(xi[order], P[order], color=col, alpha=.5); ax.axvline(0, color='k', lw=1)
    ax.set_xlim(-6, 6); ax.set_title(lab, fontsize=9.5); ax.set_xlabel('frequency (− | +)')
axs[0].set_ylabel('strength')
fig.suptitle('7. Known group theory (Gelfand–Naimark 1947), new reading: stretch-and-shift maps that keep orientation keep the sign of every frequency.\n'
             'Mult-add maps past the wall −1/k flip orientation, and they flip the frequency sign exactly as swapping i for −i does.', fontsize=11)
fig.tight_layout(); fig.savefig(OUT + 'fig7-the-wall-flips-i.png', bbox_inches='tight'); plt.close(fig)

# ---------- fig8: holonomy ----------
fig, axs = plt.subplots(1, 2, figsize=(11, 4.6))
ax = axs[0]
sq = np.array([(0, 0), (1.2, 0), (1.2, .8), (0, .8), (0, 0)])
ax.fill(sq[:, 0], sq[:, 1], color=BLUE, alpha=.2); ax.plot(sq[:, 0], sq[:, 1], color=BLUE, lw=2)
for p, q in zip(sq[:-1], sq[1:]): arrow(ax, p, p + .55*(q - p), BLUE)
ax.text(.6, .4, 'area A', ha='center', fontsize=13, color=BLUE); ax.set_aspect('equal'); ax.set_xlim(-.3, 1.5); ax.set_ylim(-.3, 1.1)
ax.set_title('Go round a loop (shift, boost, unshift, unboost)'); ax.set_axis_off()
ax = axs[1]
t = np.linspace(0, 2*np.pi, 200); ax.plot(np.cos(t), np.sin(t), color=GREY)
A = .96
ax.plot([1], [0], 'ko', ms=9); ax.text(1.08, -.1, '1 (start)')
ax.plot([np.cos(A)], [np.sin(A)], 'o', color=BLUE, ms=11); ax.text(1.1*np.cos(A), 1.1*np.sin(A), 'e^(iA): counter-clockwise loop', color=BLUE)
ax.plot([np.cos(A)], [-np.sin(A)], 'o', color=ORANGE, ms=11); ax.text(1.1*np.cos(A), -1.15*np.sin(A), 'e^(−iA): same loop, reversed', color=ORANGE)
ax.add_patch(Arc((0, 0), .5, .5, theta1=0, theta2=np.degrees(A), color=BLUE, lw=2)); ax.add_patch(Arc((0, 0), .7, .7, theta1=-np.degrees(A), theta2=0, color=ORANGE, lw=2))
ax.set_aspect('equal'); ax.set_axis_off(); ax.set_xlim(-1.3, 2.6); ax.set_ylim(-1.3, 1.3)
ax.set_title('The loop comes back turned by its area')
fig.suptitle('8. Known (Aharonov–Bohm 1959, Kostant–Souriau 1970, Berry 1984): carrying a state round an oriented area turns it by that area.\n'
             'The twist sets the angle and its orientation sets the direction of the turn.', fontsize=11.5)
fig.tight_layout(); fig.savefig(OUT + 'fig8-loop-turns-by-its-area.png', bbox_inches='tight'); plt.close(fig)

# ---------- fig9: map of routes ----------
fig, ax = plt.subplots(figsize=(13, 7.6)); ax.set_axis_off(); ax.set_xlim(0, 13); ax.set_ylim(0, 8.4)
def box(x, y, w, h, txt, col, fs=9.5, bold=False):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.08', fc=col, ec='k', alpha=.9))
    ax.text(x + w/2, y + h/2, txt, ha='center', va='center', fontsize=fs, weight='bold' if bold else 'normal', wrap=True)
box(4.9, 3.25, 3.2, 1.1, 'THE TWIST\noriented area (order dial: a·d − b·c)', '#ffe9b3', 11, True)
routes = [
 (0.2, 6.0, '+ a LENGTH\n→ J with J² = −1\n(compatible triples, Kähler)', '#d8ecd8', 'fig 2, 4: the order dial\'s group is the\nhyperbolic half-plane; its own length gives J'),
 (4.9, 6.15, '+ KEEP SIZE when composing\n→ twist acts by e^(iθ)\n(unitarity, Stone–von Neumann)', '#d8ecd8', 'fig 5'),
 (9.6, 6.0, '+ CARRY ROUND A LOOP\n→ holonomy e^(iA)\n(Kostant–Souriau, Aharonov–Bohm, Berry)', '#d8ecd8', 'fig 8'),
 (0.2, 0.55, 'CLIFFORD ALGEBRA\nunit oriented plane e₁e₂ squares to −1\n(Hestenes: i = spin plane)', '#dfe6f7', 'project Clifford page:\nx✳y − y✳x = 2 x∧y'),
 (4.9, 0.4, 'REAL OBSERVABLES\nAB − BA is antisymmetric, needs a J\n(Stueckelberg 1960)', '#dfe6f7', 'fig 6'),
 (9.6, 0.55, 'WHY C, NOT R OR H?\ncomposite systems decide\n(Hardy, Wootters, Renou 2021;\nstill argued in 2026)', '#f6dede', 'not settled by the twist'),
]
for x, y, txt, col, note in routes:
    box(x, y, 3.2, 1.25, txt, col)
    cx, cy = x + 1.6, y + (0 if y > 3 else 1.25)
    arrow(ax, (6.5, 4.35 if y > 3 else 3.25), (cx, cy + (.05 if y > 3 else .1)), GREY, lw=1.5)
    ax.text(x + 1.6, y + (1.55 if y > 3 else -.2), note, ha='center', fontsize=8, color='#333', style='italic', va='bottom' if y > 3 else 'top')
ax.text(6.5, 8.25, '9. Six known roads from the twist to i. Each needs ONE extra ingredient besides the twist (green); none gives ħ.',
        ha='center', fontsize=12, weight='bold')
fig.savefig(OUT + 'fig9-roads-from-twist-to-i.png', bbox_inches='tight'); plt.close(fig)
print('ok')

# ================= path integral section =================
from scipy.linalg import expm
# ---------- fig10: a chain of steps is a sum over paths ----------
fig, axs = plt.subplots(1, 2, figsize=(13, 5.4), gridspec_kw={'width_ratios': [1, 1.25]})
ax = axs[0]; ax.set_axis_off(); ax.set_xlim(0, 10); ax.set_ylim(0, 6)
for i, lab in enumerate(['a₁', 'a₂', 'a₃', 'a₄']):
    ax.add_patch(FancyBboxPatch((0.4 + 2.3*i, 3.2), 1.6, 1.0, boxstyle='round,pad=0.05', fc='#dfe6f7', ec=BLUE))
    ax.text(1.2 + 2.3*i, 3.7, f'1 + k·{lab}', ha='center', va='center', fontsize=12)
    if i < 3: arrow(ax, (2.05 + 2.3*i, 3.7), (2.65 + 2.3*i, 3.7), GREY)
ax.text(5, 4.9, 'time →   each short step is one stretch-and-shift', ha='center', fontsize=11)
ax.text(5, 2.3, 'a₄ ∘ₖ a₃ ∘ₖ a₂ ∘ₖ a₁  =  ( (1+ka₄)(1+ka₃)(1+ka₂)(1+ka₁) − 1 ) / k', ha='center', fontsize=11.5, color=BLUE, weight='bold')
ax.text(5, 1.35, 'Multiply out: one term for every time-ordered choice of steps\n(the skeleton of Dyson\'s series). Infinitely many tiny steps:\nVolterra\'s product integral (1887), the time-ordered exponential.', ha='center', fontsize=10)
ax.set_title('Whole-chain composition, nothing lost', fontsize=12)
ax = axs[1]
rng2 = np.random.default_rng(7); Lr, Nt = 6, 8
for x_ in range(Lr):
    for t_ in range(Nt + 1): ax.plot(t_, x_, 'o', color='#ccc', ms=5)
cols = [BLUE, ORANGE, GREEN, PURPLE, RED]
for c_ in cols:
    xs_ = [0]
    for t_ in range(Nt - 1): xs_.append(min(Lr - 1, max(0, xs_[-1] + rng2.choice([-1, 0, 1]))))
    xs_.append(3)
    if abs(xs_[-2] - 3) > 1: xs_[-2] = 2 if xs_[-2] < 3 else 4
    ax.plot(range(Nt + 1), xs_, '-', color=c_, lw=2.2, alpha=.85)
ax.plot([0], [0], 'ko', ms=10); ax.plot([Nt], [3], 'ks', ms=10)
ax.text(-.3, -.6, 'start', fontsize=10); ax.text(Nt - .3, 3.35, 'end', fontsize=10)
ax.set_xlabel('time step'); ax.set_ylabel('position on a lattice'); ax.set_xticks(range(Nt + 1)); ax.set_yticks(range(Lr))
ax.set_title('Pick positions between the steps and the same product becomes\na SUM OVER PATHS: each path weighted by its step factors (check P3)', fontsize=11)
fig.suptitle('10. A chain of weighted-circle steps, expanded, is literally a sum over paths. This part is a theorem (matrix multiplication), and it is known.', fontsize=12)
fig.tight_layout(); fig.savefig(OUT + 'fig10-chain-is-a-sum-over-paths.png', bbox_inches='tight'); plt.close(fig)

# ---------- fig11: the dial k as Wick rotation ----------
fig, axs = plt.subplots(1, 2, figsize=(12.5, 5.2))
ax = axs[0]
ax.axhline(0, color='k', lw=1); ax.axvline(0, color='k', lw=1)
t = np.linspace(np.pi/2, np.pi, 50); ax.plot(1.2*np.cos(t), -1.2*np.sin(t) * -1 * -1, color=GREY, lw=1.5, ls='--')
ax.plot([-1.2], [0], 'o', color=RED, ms=13); ax.text(-1.25, .15, 'k = −1/ħ (real)\nheat kernel, Feynman–Kac,\nstatistical mechanics\n(rigorous)', color=RED, ha='center', va='bottom', fontsize=10)
ax.plot([0], [-1.2], 'o', color=GREEN, ms=13); ax.text(.12, -1.25, 'k = −i/ħ (imaginary)\nFeynman\'s quantum propagator', color=GREEN, va='center', fontsize=10)
ax.plot([.6], [-.6], 'o', color=ORANGE, ms=10); ax.text(.68, -.5, 'any other k:\nsize grows or shrinks', color=ORANGE, fontsize=9.5)
ax.annotate('Wick rotation =\nturning the dial', xy=(-.85, -.85), xytext=(-1.75, -1.55), fontsize=10, arrowprops=dict(arrowstyle='->', color=GREY))
ax.set_xlim(-2, 1.8); ax.set_ylim(-1.8, 1.2); ax.set_aspect('equal'); ax.set_xlabel('real part of k'); ax.set_ylabel('imaginary part of k')
ax.set_title('Where the dial k points decides the physics')
ax = axs[1]
Ns = np.unique(np.logspace(0, 3, 40).astype(int)); h_ = 1.0; T_ = 2.0
for kk, col, lab in [(-1, RED, 'k = −1: shrinks to e^(−2)'), (-1j, GREEN, 'k = −i: size → 1 (kept)'), (np.exp(-1j*np.pi/3), ORANGE, 'k = e^(−iπ/3) (real part > 0): grows to e^1')]:
    ax.plot(Ns, [abs((1 + kk*h_*T_/N)**N) for N in Ns], 'o-', ms=3, color=col, label=lab)
ax.set_xscale('log'); ax.axhline(1, color='k', lw=.8, ls=':'); ax.set_xlabel('number of time slices N'); ax.set_ylabel('size of the chain (1 + k h T/N)^N')
ax.legend(); ax.set_title('Only an imaginary dial keeps size as the slices get fine (P6)')
fig.suptitle('11. Known: real-time and imaginary-time path integrals differ by Wick rotation. New reading: that rotation is a turn of the dial k,\n'
             'and the size-keeping rule is what points it straight down the imaginary axis. Its length 1/ħ is NOT fixed by anything here.', fontsize=11)
fig.tight_layout(); fig.savefig(OUT + 'fig11-dial-as-wick-rotation.png', bbox_inches='tight'); plt.close(fig)

# ---------- fig12: Levy area ----------
rng3 = np.random.default_rng(11)
vs = rng3.normal(size=(9, 2)) * [1, 1.1] + [0.35, 0.15]
pts = np.vstack([[0, 0], np.cumsum(vs, axis=0)])
fig, axs = plt.subplots(1, 2, figsize=(12.5, 5.2))
ax = axs[0]
poly = np.vstack([pts, pts[:1]])
ax.fill(poly[:, 0], poly[:, 1], color=BLUE, alpha=.18)
ax.plot(pts[:, 0], pts[:, 1], '-o', color=BLUE, lw=2, ms=4); ax.plot([pts[0, 0], pts[-1, 0]], [pts[0, 1], pts[-1, 1]], '--', color=GREY)
area = 0.5*np.sum(poly[:-1, 0]*poly[1:, 1] - poly[1:, 0]*poly[:-1, 1])
ax.text(pts[0, 0], pts[0, 1] - .4, 'start', ha='center'); ax.text(pts[-1, 0], pts[-1, 1] + .3, 'end', ha='center')
ax.set_title(f'A path of 9 steps. Signed area between path and chord\n(Lévy area) = {area:.2f} = half the sum of every pair\'s twist (P8)')
ax.set_aspect('equal')
ax = axs[1]; ax.set_axis_off(); ax.set_xlim(0, 10); ax.set_ylim(0, 6)
rows = [('Steps v₁, v₂, …, v_N', 'the path', '#eee'),
        ('Compose them, keeping only the total step\nand the total twist (Heisenberg group)', 'Chen\'s iterated integrals, level 2\n(Chen 1957; Lyons\' signatures)', '#dfe6f7'),
        ('Represent with size kept', 'phase e^(i·area/ħ)\n= the p·dq part of Feynman\'s action', '#d8ecd8')]
for i, (l, r, c) in enumerate(rows):
    y = 4.7 - 1.85*i
    ax.add_patch(FancyBboxPatch((0.2, y - .55), 4.6, 1.1, boxstyle='round,pad=0.05', fc=c, ec='k'))
    ax.add_patch(FancyBboxPatch((5.2, y - .55), 4.6, 1.1, boxstyle='round,pad=0.05', fc=c, ec='k'))
    ax.text(2.5, y, l, ha='center', va='center', fontsize=10); ax.text(7.5, y, r, ha='center', va='center', fontsize=10)
    if i < 2: arrow(ax, (2.5, y - .6), (2.5, y - 1.25), GREY)
ax.set_title('The twist summed along a path, and where it lands in physics')
fig.suptitle('12. Known (Lévy 1940s, Chen 1957, Weyl 1927): a path\'s accumulated twist is its signed area, and with size kept it is a quantum phase.', fontsize=12)
fig.tight_layout(); fig.savefig(OUT + 'fig12-twist-along-a-path.png', bbox_inches='tight'); plt.close(fig)

# ---------- fig13: convergence ----------
rng4 = np.random.default_rng(5)
Hr = rng4.normal(size=(4, 4)); Hr = (Hr + Hr.T)/2
A1 = rng4.normal(size=(3, 3)); A2 = rng4.normal(size=(3, 3)); A1, A2 = (A1 - A1.T)/2, (A2 - A2.T)/2
Ns = np.unique(np.logspace(0.5, 3.3, 25).astype(int))
eq = [np.linalg.norm(np.linalg.matrix_power(np.eye(4) - 1j*Hr/N, N) - expm(-1j*Hr)) for N in Ns]
eh = [np.linalg.norm(np.linalg.matrix_power(np.eye(4) - Hr/N, N) - expm(-Hr)) for N in Ns]
et = [np.linalg.norm(np.linalg.matrix_power(expm(A1/N) @ expm(A2/N), N) - expm(A1 + A2)) for N in Ns]
ec = [np.linalg.norm(A1 @ A2 - A2 @ A1)/(2*N) for N in Ns]
fig, ax = plt.subplots(figsize=(9, 5.4))
ax.loglog(Ns, eq, 'o-', color=GREEN, ms=4, label='chain with dial k = −i  vs  e^(−iHt)')
ax.loglog(Ns, eh, 's-', color=RED, ms=4, label='chain with dial k = −1  vs  e^(−Ht)')
ax.loglog(Ns, et, '^-', color=PURPLE, ms=4, label='Trotter slicing of two non-commuting steps')
ax.loglog(Ns, ec, '--', color=PURPLE, lw=1.5, label='prediction: twist |[A,B]| / (2N)')
ax.set_xlabel('number of time slices N'); ax.set_ylabel('error'); ax.legend(fontsize=9)
ax.set_title('13. All three shrink like 1/N (checks P4, P5, P7). The slicing error is led by the twist itself:\n'
             'order costs exactly the commutator, as on the order-dial page (Trotter 1959, Chernoff 1968).', fontsize=11)
fig.tight_layout(); fig.savefig(OUT + 'fig13-slices-converge.png', bbox_inches='tight'); plt.close(fig)

# ---------- fig14: what the framework supplies vs what Feynman needs ----------
fig, ax = plt.subplots(figsize=(12.5, 6.6)); ax.set_axis_off(); ax.set_xlim(0, 12.5); ax.set_ylim(0, 7.2)
items = [
 ('Compose short steps, whole chain, order kept', 'YES: the weighted circle fold IS the product integral (P1, P2)', GREEN),
 ('Expand into a sum over paths', 'YES: insert positions between steps (P3); known', GREEN),
 ('Order cost per slice', 'YES: the twist [A,B] leads the slicing error (P7)', GREEN),
 ('Phase from oriented area (the p·dq term)', 'YES, once size is kept: Lévy area → e^(i·area/ħ) (P8)', GREEN),
 ('The i', 'PARTLY: size-keeping forces the dial imaginary (P6); that rule is an input', ORANGE),
 ('The size of ħ', 'NO: nothing here fixes |k| = 1/ħ', RED),
 ('The dynamics (H, the potential, the action)', 'NO: supplied from outside; the framework only chains it', RED),
 ('Measure on paths, rigour in real time', 'NO: a known open-ended problem (only imaginary time is fully rigorous)', RED),
]
for i, (l, r, c) in enumerate(items):
    y = 6.5 - 0.82*i
    ax.add_patch(FancyBboxPatch((0.2, y - .32), 4.6, .64, boxstyle='round,pad=0.04', fc='#f4f4f4', ec='#999'))
    ax.add_patch(FancyBboxPatch((5.1, y - .32), 7.2, .64, boxstyle='round,pad=0.04', fc=c, ec='#999', alpha=.28))
    ax.text(2.5, y, l, ha='center', va='center', fontsize=10); ax.text(8.7, y, r, ha='center', va='center', fontsize=10)
ax.text(2.5, 7.05, 'What a path integral needs', ha='center', fontsize=11.5, weight='bold')
ax.text(8.7, 7.05, 'Does the framework supply it?', ha='center', fontsize=11.5, weight='bold')
fig.savefig(OUT + 'fig14-path-integral-scorecard.png', bbox_inches='tight'); plt.close(fig)
print('path figs ok')
