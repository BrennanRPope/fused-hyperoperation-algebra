"""Figure for 'The inverse of the Osgood test' (Nov 25, 2026). E(x)=e^x-1, L(y)=log(1+y)."""
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
BLUE, RED, GREEN, INK, SOFT = '#1d4e89', '#c1121f', '#2f6f3e', '#10131a', '#5a6070'
plt.rcParams.update({'font.family': 'serif', 'font.size': 10.5, 'axes.edgecolor': SOFT, 'axes.labelcolor': INK,
                     'xtick.color': SOFT, 'ytick.color': SOFT, 'axes.spines.top': False, 'axes.spines.right': False})
ln10 = np.log(10)
fig = plt.figure(figsize=(12, 13.2))
gs = fig.add_gridspec(3, 2, hspace=0.42, wspace=0.22)

# (a) the clock: time to reach height y, plotted against u = log10 log10 (1+y)
ax = fig.add_subplot(gs[0, 0])
u = np.geomspace(1e-2, 300, 2000)
Ly = 10.0 ** u * ln10                    # L(y) = log(1+y), with log10(1+y) = 10^u
curves = [
    (np.log1p(Ly), BLUE, '-', r"$f_2=(1+y)(1+L y)$: clock $L^2 y$"),
    (np.log1p(np.log1p(Ly)), BLUE, '--', r"$f_3=f_2\,(1+L^2y)$: clock $L^3 y$"),
    ((1 - (1 + Ly) ** -0.3) / 0.3, RED, '-', r"$f_2\,(1+Ly)^{0.3}$: stops at $1/0.3$"),
    (2 * (1 - np.exp(-0.5 * np.minimum(Ly, 700))), RED, '--', r"$(1+y)^{1.5}$: stops at 2"),
    (1 - np.exp(-np.minimum(Ly, 700)), RED, ':', r"$(1+y)^{2}$: stops at 1"),
]
for T, c, ls, lab in curves:
    ax.plot(u, T, color=c, ls=ls, lw=2, label=lab)
for v in [1, 2, 1 / 0.3]:
    ax.axhline(v, color=RED, lw=0.6, alpha=0.4)
ax.set_xscale('log')
ax.set_xlabel(r"height reached: $\log_{10}\log_{10}(1+y)$  (300 means $y\approx10^{10^{300}}$)")
ax.set_ylabel(r"time to get there, $\int_0^y ds/f(s)$")
ax.set_ylim(0, 7.4)
ax.text(0.02, 7.0, r"$f_1=1+y$: clock $\log(1+y)$, off the top at once", color=BLUE, fontsize=9.3)
ax.text(30, 5.55, "still climbing:\nnever stops", color=BLUE, fontsize=9.3)
ax.legend(loc='upper left', fontsize=8.6, frameon=False, bbox_to_anchor=(0.0, 0.92))
ax.set_title("(a) The Osgood integral is a clock", loc='left', fontweight='bold')

# (b) the solutions
ax = fig.add_subplot(gs[0, 1])
t = np.linspace(0, 3.6, 3000)
def clip(v): v = np.array(v, float); v[v > 300] = np.nan; return v
ax.plot(t, clip(t / ln10), color=BLUE, lw=2, label=r"$E(t)$ solves $y'=1+y$")
ax.plot(t, clip(np.expm1(t) / ln10), color=BLUE, lw=2, ls='--', label=r"$E^2(t)$ solves $y'=f_2(y)$")
with np.errstate(over='ignore'):
    ax.plot(t, clip(np.expm1(np.expm1(t)) / ln10), color=BLUE, lw=2, ls=':', label=r"$E^3(t)$ solves $y'=f_3(y)$")
with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
    t1 = 1 - np.geomspace(1, 1e-60, 3000)
    ax.plot(t1, clip(-np.log10(1 - t1)), color=RED, lw=2, ls=':', label=r"$y'=(1+y)^2$: $y=t/(1-t)$")
    t2 = 2 - 2 * np.geomspace(1, 1e-30, 3000)
    ax.plot(t2, clip(-2 * np.log10(1 - t2 / 2)), color=RED, lw=2, ls='--', label=r"$y'=(1+y)^{1.5}$")
    tt = (1 / 0.3) * (1 - np.geomspace(1, 1e-3, 3000))
    ax.plot(tt, clip(((1 - 0.3 * tt) ** (-1 / 0.3) - 1) / ln10), color=RED, lw=2, label=r"$y'=f_2(y)(1+Ly)^{0.3}$")
for v in [1, 2, 1 / 0.3]:
    ax.axvline(v, color=RED, lw=0.8, alpha=0.5)
ax.set_ylim(0, 40); ax.set_xlim(0, 3.6)
ax.set_xlabel("time $t$"); ax.set_ylabel(r"$\log_{10}(1+y)$")
ax.legend(loc='upper left', fontsize=8.6, frameon=False)
ax.set_title("(b) Finite clock = blow-up at that time", loc='left', fontweight='bold')

# (c) uniqueness fan
ax = fig.add_subplot(gs[1, 0])
t = np.linspace(0, 3, 600)
ax.plot(t, 0 * t, color=RED, lw=2.4)
for c in [0, 0.5, 1.0, 1.5, 2.0]:
    x = np.where(t >= c, (t - c) ** 2, 0)
    ax.plot(t, x, color=RED, lw=1.6, alpha=0.85)
ax.plot(t, 0.02 * np.exp(t) * 0, color=BLUE, lw=0)
ax.text(0.05, 3.5, r"$x'=2\sqrt{x}$, $x(0)=0$: every $x_c=(t-c)^2$" "\n" r"(zero until $c$) is a solution." "\n"
        r"Clock out of 0: $\int_0^r ds/2\sqrt{s}=\sqrt{r}$, finite.", fontsize=9.5, color=RED, va='top')
ax.text(1.55, 0.35, r"Lipschitz $x'=x$ or Yudovich $x'=x\log(1/x)$:" "\n" r"clock out of 0 is infinite, only $x\equiv0$",
        fontsize=9.2, color=BLUE)
ax.set_ylim(-0.2, 4); ax.set_xlabel("time $t$"); ax.set_ylabel("$x$")
ax.set_title("(c) Finite clock out of 0 = non-uniqueness", loc='left', fontweight='bold')

# (d) direction decides
ax = fig.add_subplot(gs[1, 1])
g = np.linspace(-1.6, 1.6, 15); X, Y = np.meshgrid(g, g); R2 = X ** 2 + Y ** 2
m = (R2 > 0.05) & (R2 < 2.4)
U, V = -Y * R2, X * R2; S = np.hypot(U, V) + 1e-9
left = (X < -0.05) & m; right = (X > 0.05) & m
ax.quiver(X[left], Y[left], (U / S)[left], (V / S)[left], color=BLUE, scale=22, width=0.004)
U2, V2 = X * R2, Y * R2; S2 = np.hypot(U2, V2) + 1e-9
ax.quiver(X[right], Y[right], (U2 / S2)[right], (V2 / S2)[right], color=RED, scale=22, width=0.004)
th = np.linspace(np.pi / 2, 3 * np.pi / 2, 200)
ax.plot(np.cos(th), np.sin(th), color=BLUE, lw=2.4)
ax.plot([0.0, 1.55], [0, 0], color=RED, lw=2.4)
ax.annotate('', xy=(1.62, 0), xytext=(1.3, 0), arrowprops=dict(arrowstyle='->', color=RED, lw=2.4))
ax.axvline(0, color=SOFT, lw=0.6)
ax.text(-1.72, 1.85, r"$|x|^2 Jx$: speed $r^3$," "\n" r"all of it sideways;" "\n" r"$x\cdot F=0$, size never changes", color=BLUE, fontsize=8.8)
ax.text(0.12, -1.85, r"$|x|^2 x$: speed $r^3$," "\n" r"all of it outward;" "\n" r"$x\cdot F=r^4$, blow-up", color=RED, fontsize=8.8, va="top")
ax.set_aspect('equal'); ax.set_xlim(-1.75, 1.75); ax.set_ylim(-2.6, 2.6); ax.set_xticks([]); ax.set_yticks([])
ax.spines['left'].set_visible(False); ax.spines['bottom'].set_visible(False)
ax.set_title("(d) Systems: same size, opposite fates", loc='left', fontweight='bold')

# (e) the rung shift preserves the area (the test moves down a rung)
X4 = 4.0
ax = fig.add_subplot(gs[2, 0])
s = np.linspace(0, np.expm1(X4), 800)
yv = 1 / ((1 + s) * (1 + np.log1p(s)))
ax.fill_between(s, yv, color=GREEN, alpha=0.25); ax.plot(s, yv, color=GREEN, lw=2)
ax.set_xlabel(r"$s$  (up to $E(4)=e^4-1\approx 53.6$)"); ax.set_ylabel(r"$1/f(s)$")
ax.text(14, 0.55, r"$f(s)=(1+s)\,g(s)$, $g=1+L(s)$" "\n" r"shaded area $=\log 5\approx1.609$", color=GREEN, fontsize=10)
ax.set_title("(e) The test on rung $s$ …", loc='left', fontweight='bold')
ax = fig.add_subplot(gs[2, 1])
x = np.linspace(0, X4, 400); yv = 1 / (1 + x)
ax.fill_between(x, yv, color=GREEN, alpha=0.25); ax.plot(x, yv, color=GREEN, lw=2)
ax.set_xlabel(r"$x$ with $s=E(x)$, $1+s=e^x$"); ax.set_ylabel(r"$1/g(E(x))$")
ax.text(1.3, 0.62, r"$g(E(x))=1+x$" "\n" r"same area $=\log 5$: $ds=(1+s)\,dx$" "\n" r"the lift factor $1+s$ is used up", color=GREEN, fontsize=10)
ax.set_title("(e) … is the same test one rung down", loc='left', fontweight='bold')

fig.savefig('/mnt/project-files/math-theory/osgood-inverse/figures/osgood_inverse.png', dpi=130, bbox_inches='tight', facecolor='white')
print('ok')
