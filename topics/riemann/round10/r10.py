"""Round 10: a prime's lifecycle, untried prime results on RH, and hyper integrals.
Writes figures/fig_*.png and r10.json."""
import json, math, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mpmath as mp

OUT = {}
G = np.array(json.load(open('../round9/zeros.json'))['gamma'])   # 5500 zeros, mpmath.zetazero
EG = math.exp(-0.5772156649015329)

N = 10**8
s = np.ones(N + 1, bool); s[:2] = False
for i in range(2, int(N**0.5) + 1):
    if s[i]:
        s[i * i::i] = False
P = np.nonzero(s)[0]
print('primes up to 1e8:', len(P), flush=True)

# ---------- L1 birth ----------
fig, ax = plt.subplots(figsize=(12, 3.6))
X = 130
for j, p in enumerate([2, 3, 5, 7]):
    m = np.arange(2 * p, X + 1, p)
    ax.plot(m, np.full(len(m), j), '|', ms=14, color=f'C{j}')
    ax.text(-6, j, f'hops of {p}', va='center', ha='right', fontsize=8)
hit = np.zeros(X + 1, bool)
for p in [2, 3, 5, 7]:
    hit[2 * p::p] = True
free = [n for n in range(8, X + 1) if not hit[n]]
ax.plot(free, np.full(len(free), 4.2), 'k*', ms=8)
ax.text(-6, 4.2, 'untouched', va='center', ha='right', fontsize=8)
ax.annotate('11 is born: the first\nnumber no hop lands on', (11, 4.2), (20, 5.1), arrowprops=dict(arrowstyle='->'), fontsize=8)
ax.annotate('121 = 11²: the first number\nonly 11 can catch', (121, 4.2), (95, 5.1), arrowprops=dict(arrowstyle='->'), fontsize=8)
ax.set_xlim(-30, X + 2); ax.set_ylim(-0.7, 5.9); ax.set_yticks([]); ax.set_xlabel('n')
ax.set_title('Stage 1, birth: a prime is the first number the smaller primes skip (every untouched number below 121 is prime)')
fig.tight_layout(); fig.savefig('figures/fig_life1_birth.png', dpi=130); plt.close(fig)

# ---------- L2 first catch at p^2 ----------
lpf_share = {}
fig, ax = plt.subplots(figsize=(9, 4.5))
for p, c in [(11, 'C0'), (101, 'C1'), (1009, 'C2')]:
    # numbers n <= x with least prime factor p (n composite): n = p*m, m >= p, m has no prime < p
    small = [q for q in P[P < p]]
    xs = np.unique(np.logspace(math.log10(p * p * 0.5), 8, 300).astype(np.int64))
    m = np.arange(p, N // p + 1)
    ok = np.ones(len(m), bool)
    for q in small:
        ok &= (m % q != 0)
    vals = p * m[ok]
    cnt = np.searchsorted(vals, xs, side='right') / xs
    limit = (1 / p) * np.prod([1 - 1 / q for q in small])
    ax.loglog(xs, np.maximum(cnt, 1e-9), color=c, label=f'p = {p}')
    ax.axhline(limit, color=c, ls=':', lw=1)
    ax.axvline(p * p, color=c, ls='--', lw=0.6)
    lpf_share[p] = {'limit': limit, 'at_1e8': float(cnt[-1])}
ax.set_ylim(1e-6, 0.1)
ax.set_xlabel('N'); ax.set_ylabel('share of n ≤ N whose smallest prime factor is p (composites)')
ax.set_title('Stage 2, first catch: p starts sieving only at p² (dashed), and its share climbs toward\n(1/p)·∏_{q<p}(1 − 1/q) (dotted), slowly for large p (it has to wait for p·(next prime))', fontsize=9)
ax.legend(); fig.tight_layout(); fig.savefig('figures/fig_life2_first_catch.png', dpi=130); plt.close(fig)
OUT['lpf_share'] = lpf_share

# ---------- L3 wheel ----------
ps = P[P < 10**5]
mert = np.cumprod(1 - 1 / ps)
fig, ax = plt.subplots(figsize=(9, 4.5))
ax.semilogx(ps, mert * np.log(ps), label='open fraction of the wheel × log p  →  e^(−γ) = 0.5615 (Mertens)')
ax.semilogx(ps, mert * np.log(ps.astype(float)**2), label='open fraction × log(p²): density of the wheel at the range it certifies, over the true 1/log x')
ax.axhline(EG, color='C0', ls=':'); ax.axhline(2 * EG, color='C1', ls=':')
ax.text(12, 2 * EG + 0.02, '2e^(−γ) = 1.1229: the wheel overcounts the primes it certifies by 12%', fontsize=8, color='C1')
ax.set_xlabel('p (newest prime in the wheel)'); ax.set_ylim(0.4, 1.4)
ax.set_title('Stage 3, the wheel: p joins and closes 1/p of the open spokes. The wheel at p certifies every number below p²,\nyet its density there is 2e^(−γ)/log x, not 1/log x (this is the Legendre overcount, P8)', fontsize=9)
ax.legend(fontsize=7.5); fig.tight_layout(); fig.savefig('figures/fig_life3_wheel.png', dpi=130); plt.close(fig)
OUT['mertens_at_1e5'] = float(mert[-1] * math.log(ps[-1])); OUT['wheel_overcount_at_1e5'] = float(mert[-1] * 2 * math.log(ps[-1]))

# ---------- L4 gaps ----------
g = np.diff(P)
rec_i = []
best = 0
for i, v in enumerate(g):
    if v > best:
        best = v; rec_i.append(i)
rp = P[rec_i].astype(float); rg = g[rec_i].astype(float)
fig, ax = plt.subplots(figsize=(9, 4.8))
ax.loglog(rp, rg, 'ko', ms=4, label='record gaps up to 10⁸')
xx = np.logspace(1, 8, 200)
ax.loglog(xx, np.log(xx)**2, label='(log p)²  Cramér 1936')
ax.loglog(xx, 2 * EG * np.log(xx)**2, label='2e^(−γ)(log p)²  Granville 1995 (conjectured lower bound for the limsup)')
ax.loglog(xx, np.sqrt(xx) * np.log(xx), label='√p·log p  what RH guarantees (Cramér 1920)')
ax.loglog(xx, np.log(xx), ':', label='log p  average gap')
ax.set_xlabel('p'); ax.set_ylabel('gap after p'); ax.legend(fontsize=7.5)
ratio = rg / np.log(rp)**2
ax.set_title(f'Stage 4, gaps: a prime\'s gap to the next is log p on average; records stay below (log p)² (largest ratio here {ratio[rp>1e3].max():.3f})\nRH only promises √p log p, far above every record', fontsize=9)
fig.tight_layout(); fig.savefig('figures/fig_life4_gaps.png', dpi=130); plt.close(fig)
OUT['record_gap_max_ratio_above_1e3'] = float(ratio[rp > 1e3].max())

# ---------- L5 heard by the zeros (Landau 1912) ----------
fig, ax = plt.subplots(figsize=(9, 4.5))
spk = {}
for p, c in [(11, 'C0'), (101, 'C1'), (1009, 'C2'), (10007, 'C3')]:
    cs = np.cumsum(np.cos(G * math.log(p)))
    norm = -cs * 2 * math.pi * math.sqrt(p) / (G * math.log(p))
    ax.plot(G, norm, color=c, lw=0.8, label=f'p = {p}')
    spk[p] = float(np.mean(norm[-500:]))
ax.axhline(1, color='k', lw=0.8)
ax.set_xscale('log'); ax.set_ylim(-1, 3)
ax.set_xlabel('height T (zeros used up to T)'); ax.set_ylabel('−(2π√p / (T log p))·Σ_{γ≤T} cos(γ log p)')
ax.set_title('Stage 5, heard by the zeros: Σ cos(γ log x) over zeros up to T tends to −(T/2π)·Λ(x)/√x (Landau 1912).\nLarger primes need more zeros before their spike settles at 1', fontsize=9)
ax.legend(); fig.tight_layout(); fig.savefig('figures/fig_life5_heard.png', dpi=130); plt.close(fig)
OUT['landau_spike_last500'] = spk

# ---------- L6 the prime's share of log|zeta| on the line ----------
mp.mp.dps = 15
ts = np.linspace(5000, 5060, 1500)
lz = np.array([float(mp.log(abs(mp.zeta(mp.mpc(0.5, t))))) for t in ts])
ok = np.isfinite(lz)
cuts = [2, 3, 5, 7, 11, 31, 101, 1009, 10007]
expl = []
for Pc in cuts:
    pp = P[P <= Pc]
    approx = np.array([np.sum(-np.log(np.abs(1 - pp**(-0.5 - 1j * t)))) for t in ts])
    c = np.corrcoef(lz[ok], approx[ok])[0, 1]
    expl.append(c * c)
fig, axs = plt.subplots(1, 2, figsize=(14, 4.5))
axs[0].plot(ts, lz, lw=0.7, label='log|ζ(½+it)|')
pp = P[P <= 1009]
axs[0].plot(ts, [np.sum(-np.log(np.abs(1 - pp**(-0.5 - 1j * t)))) for t in ts], lw=0.7, label='Euler factors p ≤ 1009')
axs[0].set_ylim(-6, 4); axs[0].legend(fontsize=8); axs[0].set_xlabel('t')
axs[0].set_title('log|ζ| on the line near t = 5000 against a truncated Euler product', fontsize=9)
axs[1].semilogx(cuts, expl, 'o-')
axs[1].set_xlabel('largest prime kept'); axs[1].set_ylabel('share of the variance of log|ζ(½+it)| explained (r²)')
axs[1].set_title('Stage 6, one voice among many: each prime adds about 1/(2p) to the variance of log|ζ|. At t ≈ 5000 the\nsmall primes still carry much of it (2 alone: 17%), but the total grows like ½ log log t, so each prime\'s share fades', fontsize=9)
fig.tight_layout(); fig.savefig('figures/fig_life6_share.png', dpi=130); plt.close(fig)
OUT['variance_share'] = dict(zip(map(str, cuts), map(float, expl)))
print('lifecycle done', flush=True)

# ---------- U1 wheel gaps vs zero spacings vs prime gaps ----------
M = 2 * 3 * 5 * 7 * 11 * 13
w = np.ones(M, bool)
for q in [2, 3, 5, 7, 11, 13]:
    w[::q] = False
sp = np.nonzero(w)[0]; wg = np.diff(np.concatenate([sp, [sp[0] + M]])); wg = wg / wg.mean()
sel = (P > 10**7) & (P < 2 * 10**7)
pg = np.diff(P[sel]) / np.log(P[sel][:-1])
# unfold zeros: N(T) ~ T/2π log(T/2πe) + 7/8
Nt = G / (2 * math.pi) * np.log(G / (2 * math.pi * math.e)) + 7 / 8
zs = np.diff(Nt)
fig, ax = plt.subplots(figsize=(9, 4.8))
bins = np.linspace(0, 4, 41)
ax.hist(wg, bins, density=True, histtype='step', lw=1.5, label=f'wheel mod 30030 spoke gaps ({len(wg)})')
ax.hist(pg, bins, density=True, histtype='step', lw=1.5, label='prime gaps near 10⁷, over log p')
ax.hist(zs, bins, density=True, histtype='step', lw=1.5, label='zeta zero spacings (5500 zeros, unfolded)')
u = np.linspace(0, 4, 200)
ax.plot(u, np.exp(-u), 'k:', label='Poisson e^(−u)')
ax.plot(u, 32 / np.pi**2 * u**2 * np.exp(-4 * u**2 / np.pi), 'k--', label='GUE (Wigner surmise)')
ax.set_xlabel('gap / mean gap'); ax.legend(fontsize=8)
ax.set_title('Untried P4/P5: the wheel\'s gaps are neither Poisson nor GUE; prime gaps are near Poisson;\nzero spacings repel (GUE). The wheel\'s gap statistics say nothing about the zeros\' spacing', fontsize=9)
fig.tight_layout(); fig.savefig('figures/fig_u_wheel_vs_zeros.png', dpi=130); plt.close(fig)
OUT['wheel_gap_small_frac'] = float(np.mean(wg < 0.5)); OUT['zero_gap_small_frac'] = float(np.mean(zs < 0.5)); OUT['prime_gap_small_frac'] = float(np.mean(pg < 0.5))

# ---------- U2 Chebyshev bias ----------
bias = {}
fig, ax = plt.subplots(figsize=(9, 4.5))
for mod, a, b, c in [(4, 3, 1, 'C0'), (3, 2, 1, 'C1')]:
    d = np.cumsum(np.where(P % mod == a, 1, np.where(P % mod == b, -1, 0)))
    frac = float(np.mean(d > 0))
    xs = P[::2000]; dd = d[::2000]
    ax.semilogx(xs, dd * np.log(xs) / np.sqrt(xs), color=c, lw=0.7, label=f'mod {mod}: (π(x;{mod},{a}) − π(x;{mod},{b}))·log x/√x; ahead {100*frac:.2f}% of primes')
    bias[mod] = frac
ax.axhline(0, color='k', lw=0.6)
ax.set_xlabel('x'); ax.legend(fontsize=8)
ax.set_title('Untried O14, Chebyshev\'s bias: the non-squares lead almost always. Rubinstein–Sarnak (1994, assuming GRH and\nlinear independence) predict 99.59% (mod 4) and 99.90% (mod 3) in logarithmic density; their zeros are those of L(s, χ), not ζ', fontsize=9)
fig.tight_layout(); fig.savefig('figures/fig_u_chebyshev_bias.png', dpi=130); plt.close(fig)
OUT['chebyshev_ahead_fraction'] = bias

# ---------- U4 pair correlation (Hardy-Littlewood tuples meet the zeros) ----------
Nu = Nt
diffs = []
for k in range(1, 40):
    dd = Nu[k:] - Nu[:-k]
    diffs.append(dd[dd < 3])
diffs = np.concatenate(diffs)
fig, ax = plt.subplots(figsize=(9, 4.5))
bins = np.linspace(0, 3, 61)
h, e = np.histogram(diffs, bins)
ax.plot((e[:-1] + e[1:]) / 2, h / (len(Nu) * (e[1] - e[0])), 'o', ms=3, label='zeta zeros 1–5500')
u = np.linspace(0.001, 3, 300)
ax.plot(u, 1 - (np.sin(np.pi * u) / (np.pi * u))**2, 'k-', label='Montgomery 1973: 1 − (sin πu / πu)²')
ax.set_xlabel('u (distance in mean spacings)'); ax.set_ylabel('pair correlation')
ax.set_title('Untried P9/P10: the Hardy–Littlewood tuple constants do reach the zeros, through pair correlation.\nMontgomery\'s form, and its lower-order terms from the singular series (Bogomolny–Keating 1996), are in print', fontsize=9)
ax.legend(); fig.tight_layout(); fig.savefig('figures/fig_u_pair_correlation.png', dpi=130); plt.close(fig)

# ---------- H hyper integrals ----------
f = lambda x: np.sin(3 * x) + 1.2
I = (1 - math.cos(3)) / 3 + 1.2
res = {}
fig, axs = plt.subplots(1, 3, figsize=(16, 4.6))
ns = np.unique(np.logspace(0.5, 5, 30).astype(int))
for k, c in [(1.0, 'C0'), (-0.5, 'C1'), (2.0, 'C2')]:
    errs = []
    target = (math.exp(k * I) - 1) / k
    for n in ns:
        x = (np.arange(n) + 0.5) / n
        acc = 0.0
        for h in f(x) / n:
            acc = acc + h + k * acc * h
        errs.append(abs(acc - target))
    axs[0].loglog(ns, errs, 'o-', ms=3, color=c, label=f'k = {k}')
    res[f'rung2_k{k}'] = {'target': target, 'err_at_max_n': errs[-1]}
axs[0].loglog(ns, 1 / ns, 'k:', label='1/n')
axs[0].set_xlabel('n (slices)'); axs[0].set_ylabel('|fused sum − (e^(k∫f) − 1)/k|')
axs[0].set_title('Rung 2 (∘_k): the fused Riemann sum converges to the\nVolterra product integral (e^(k∫f) − 1)/k, error ~ 1/n', fontsize=9)
axs[0].legend(fontsize=8)
# rung 3: in log coordinates L = log(1+x), rung 3 multiplies L's
for scale, c, lab in [('slice', 'C3', 'heights × Δ: product of small logs → 0'), ('power', 'C4', 'logs raised to Δ: → exp(∫ log log(1+f))')]:
    vals = []
    for n in ns:
        x = (np.arange(n) + 0.5) / n
        Ls = np.log1p(f(x) / n) if scale == 'slice' else np.log(np.log1p(f(x))) / n
        vals.append(np.sum(np.log(Ls)) if scale == 'slice' else np.sum(Ls))
    vals = np.array(vals)
    if scale == 'slice':
        axs[1].semilogx(ns, vals, 'o-', ms=3, color=c, label=lab)
    else:
        tgt = mp.quad(lambda t: mp.log(mp.log(1 + mp.sin(3 * t) + 1.2)), [0, 1])
        axs[1].semilogx(ns, vals, 'o-', ms=3, color=c, label=lab)
        axs[1].axhline(float(tgt), color=c, ls=':')
        res['rung3_power_target_log'] = float(tgt); res['rung3_power_at_max_n'] = float(vals[-1])
axs[1].set_ylim(-60, 2)
axs[1].set_xlabel('n'); axs[1].set_ylabel('log of the folded log-height')
axs[1].set_title('Rung 3: each rung up is one more log. Fold with Δ-scaled heights\nand it collapses; fold with Δ as an exponent and it is a geometric integral', fontsize=9)
axs[1].legend(fontsize=7.5)
# on primes
ss = [2.0, 3.0]
rows = []
for s_ in ss:
    pp = P[P < 10**6].astype(float)
    a = pp**(-s_)
    r2 = np.prod(1 + a) - 1
    exact = float(mp.zeta(s_) / mp.zeta(2 * s_)) - 1
    r3 = np.sum(np.log(np.log1p(a)))      # log of prod log(1 + p^-s)
    rows.append((s_, r2, exact, r3))
res['primes'] = [{'s': r[0], 'rung2_fold': r[1], 'zeta(s)/zeta(2s)-1': r[2], 'rung3_log_of_fold': r[3]} for r in rows]
PP = P[P < 10**6].astype(float)
cum3 = np.cumsum(np.log(np.log1p(PP**-2.0)))
cum2 = np.cumprod(1 + PP**-2.0) - 1
axs[2].semilogx(PP, cum2, label='rung 2 fold of p^(−2) → ζ(2)/ζ(4) − 1 = 0.5198')
axs[2].axhline(float(mp.zeta(2) / mp.zeta(4)) - 1, color='C0', ls=':')
ax2 = axs[2].twinx(); ax2.semilogx(PP, cum3, color='C3', label='rung 3: log of ∏ log(1 + p^(−2)) → −∞'); ax2.set_ylabel('rung 3 (log scale of the fold)', color='C3')
axs[2].set_xlabel('largest prime folded'); axs[2].legend(fontsize=7.5, loc='center left'); ax2.legend(fontsize=7.5, loc='lower left')
axs[2].set_title('On the primes: rung 2 gives the Euler product ζ(s)/ζ(2s) at k = 1;\nrung 3 multiplies logs, has no Dirichlet series, and runs to 0', fontsize=9)
fig.tight_layout(); fig.savefig('figures/fig_hyper_integrals.png', dpi=130); plt.close(fig)
OUT['hyper'] = res
json.dump(OUT, open('r10.json', 'w'), indent=1, default=float)
print(json.dumps(OUT, indent=1, default=float)[:3000])
