"""Pictures of the parent in its full space, and the bath.

Inputs tracked: Re z, Im z (field), t (heat-flow time), spins N with coupling J, bath field h.
Height: log|H| (the pressure); colour where used: arg H (the phase).
Kernels: Riemann Phi; p-ball p=3, a=2 (holds, a <= p); p-ball p=3, a=4 (breaks, a > p).
Writes figures/fig_*.png and landscape.json.
"""
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import hsv_to_rgb
import kern

OUT = 'figures'
import os
os.makedirs(OUT, exist_ok=True)
RES = {}

R = kern.riemann(500)
PH = kern.pball(3, 2, 600)
PB = kern.pball(3, 4, 800)
KS = [(R, 'Riemann Φ', 0.0, 72.0), (PH, 'p-ball p=3, a=2 (holds)', 0.0, 40.0), (PB, 'p-ball p=3, a=4 (breaks)', 70.0, 125.0)]


def envelope(k, x, t=0.0):
    """Smooth size of H along the real line, to flatten the landscape: mean |H| over a vertical band."""
    ys = np.linspace(-3, 3, 13)
    A = np.abs(k.H(x[None, :] + 1j * ys[:, None], t))
    return np.exp(np.mean(np.log(A + 1e-300), axis=0))


def domain_rgb(F, E):
    hue = (np.angle(F) / (2 * np.pi)) % 1.0
    L = np.log(np.abs(F) / E + 1e-300)
    val = np.clip(0.5 + 0.12 * L, 0.0, 1.0)
    return hsv_to_rgb(np.dstack([hue, 0.75 * np.ones_like(hue), val]))


# ---------- 1. landscape: |H| as height over the plane, zeros touch the floor ----------
fig = plt.figure(figsize=(16, 5.2))
for i, (k, name, x0, x1) in enumerate(KS):
    ax = fig.add_subplot(1, 3, i + 1, projection='3d')
    xs = np.linspace(x0 + 0.1, x1, 260)
    ys = np.linspace(-2.5, 2.5, 121)
    X, Y = np.meshgrid(xs, ys)
    F = k.H(X + 1j * Y)
    E = envelope(k, xs)[None, :]
    Z = np.clip(np.log10(np.abs(F) / E + 1e-300), -3, 1.2)
    ax.plot_surface(X, Y, Z, cmap='viridis', linewidth=0, antialiased=True, rstride=2, cstride=2, alpha=0.95)
    zz = kern.complex_zeros(k, x1, 2.5, nx=int((x1) * 10)) if i == 2 else kern.real_zeros(k, x1) + 0j
    zz = zz[(zz.real > x0) & (zz.real < x1)]
    allz = np.concatenate([zz, np.conj(zz[np.abs(zz.imag) > 1e-6])])
    ax.scatter(allz.real, allz.imag, -3 * np.ones(len(allz)), color='red', s=18, depthshade=False)
    ax.set_xlabel('Re z'); ax.set_ylabel('Im z'); ax.set_zlabel('log10 |H| (flattened)')
    ax.set_title(name + '\nred dots: zeros on the floor')
    ax.view_init(elev=32, azim=-60)
    RES[f'landscape_zeros_{i}'] = [[float(z.real), float(z.imag)] for z in zz]
fig.suptitle('The parent as a landscape at heat time t = 0: height is log|H|, zeros are where the surface touches the floor')
fig.tight_layout()
fig.savefig(f'{OUT}/fig_landscape_3d.png', dpi=130)
plt.close(fig)

# ---------- 2. heat-time dimension: zero curves in (x, t) ----------
fig, axs = plt.subplots(1, 3, figsize=(16, 5.2))
# Riemann: real zeros below x = 72, rewind t to -6
ts = np.linspace(-6, 0.3, 64)
pts = []
for t in ts:
    for z in kern.real_zeros(R, 72, t=t):
        pts.append((z, t))
pts = np.array(pts)
axs[0].plot(pts[:, 0], pts[:, 1], 'k.', ms=2.5)
axs[0].axhline(0, color='tab:blue', lw=1); axs[0].axhspan(0, 0.2, color='tab:green', alpha=0.15)
axs[0].text(45, 0.25, 'Λ in [0, 0.2]', fontsize=8, color='darkgreen')
axs[0].set_title('Riemann Φ: real zeros vs heat time t\n(4th and 5th meet near t ≈ −4 on rewind)')
axs[0].set_xlabel('real zero z (= 2γ)'); axs[0].set_ylabel('heat time t (negative = rewind)')
RES['riemann_heat_counts'] = {f'{t:.2f}': int(np.sum(np.isclose(pts[:, 1], t))) for t in ts[::8]}
# p-ball holds: rewind to below Lambda ~ -2.5
ts = np.linspace(-3.2, 1.0, 50)
for t in ts:
    zr = kern.real_zeros(PH, 40, t=t)
    axs[1].plot(zr, t * np.ones_like(zr), 'k.', ms=2.5)
lam = json.load(open('lam_pball.json')).get('3|2', {}).get('Lambda')
if lam is not None:
    axs[1].axhline(lam, color='tab:red', lw=1)
    axs[1].text(20, lam - 0.25, f'measured Λ ≈ {lam:.2f}', fontsize=8, color='tab:red')
axs[1].axhline(0, color='tab:blue', lw=1)
axs[1].set_title('p-ball p=3, a=2 (holds): real zeros vs t\n(real down to Λ ≈ −2.50, then a pair meets)')
axs[1].set_xlabel('real zero z'); axs[1].set_ylabel('heat time t')
# p-ball breaks: imaginary part of first off-line zero vs t
ts = np.linspace(-1.0, 3.0, 17)
first_im = []
for t in ts:
    zz = kern.complex_zeros(PB, 125, 3, t=t, nx=1300)
    off = zz[np.abs(zz.imag) > 1e-6]
    sc = axs[2].scatter(off.real, t * np.ones(len(off)), c=off.imag, cmap='plasma', vmin=0, vmax=1.3, s=12)
    first_im.append([float(t), float(off.real.min()) if len(off) else None, float(off.imag[np.argmin(off.real)]) if len(off) else None])
plt.colorbar(sc, ax=axs[2], label='Im z of off-line zero')
axs[2].set_title('p-ball p=3, a=4 (breaks): off-line zeros vs t\n(forward flow pushes them past x = 125; Λ = +∞)')
axs[2].set_xlabel('Re z'); axs[2].set_ylabel('heat time t')
RES['broken_first_offline_vs_t'] = first_im
fig.suptitle('Heat-flow time is a real extra dimension: forward flow keeps real zeros real (de Bruijn), rewind loses them at Λ')
fig.tight_layout()
fig.savefig(f'{OUT}/fig_heat_dimension.png', dpi=130)
plt.close(fig)


# ---------- 3. spin dimension: two coupled spins, field h on both ----------
def two_spin_diag(k, hs, J):
    s = np.concatenate([-k.u, k.u]); w = np.concatenate([k.w, k.w])
    W = np.outer(w, w) * np.exp(J * np.outer(s, s))
    E = np.exp(1j * np.outer(hs, s))
    return np.sum((E @ W) * E, axis=1).real


def count_real(vals):
    return int(np.sum(np.sign(vals[:-1]) * np.sign(vals[1:]) < 0))


fig, axs = plt.subplots(1, 2, figsize=(13, 5))
spin = {}
for ax, (k, name, Js, hmax) in zip(axs, [(PH, 'p-ball p=3, a=2', np.linspace(-12, 12, 49), 30.0),
                                         (R, 'Riemann Φ', np.linspace(-300, 300, 49), 40.0)]):
    hs = np.linspace(0.02, hmax, 3000)
    rows = []
    for J in Js:
        v = two_spin_diag(k, hs, J)
        i = np.where(np.sign(v[:-1]) * np.sign(v[1:]) < 0)[0]
        ax.plot(hs[i], J * np.ones(len(i)), 'k.', ms=2.5)
        rows.append([float(J), len(i)])
    ax.axhline(0, color='tab:blue', lw=1)
    ax.axhspan(0, Js.max(), color='tab:green', alpha=0.08)
    ax.axhspan(Js.min(), 0, color='tab:red', alpha=0.06)
    ax.set_ylim(Js.min(), Js.max())
    ax.text(0.5, Js.max() * 0.82, 'J > 0: each double zero splits into two real zeros\n(Lee–Yang–Newman: a theorem for the p-ball; for Φ it would follow from RH)', fontsize=8, color='darkgreen', backgroundcolor='white')
    ax.text(0.5, Js.min() * 0.5, 'J < 0: no real zeros at all; each double zero\nsplits into a complex pair (Z ≈ G² − J·G′²)', fontsize=8, color='darkred', backgroundcolor='white')
    ax.text(hmax * 0.55, Js.max() * 0.03, 'J = 0: Z = G(h)², double zeros', fontsize=8, color='tab:blue')
    ax.set_xlabel('field h on both spins (real)'); ax.set_ylabel('coupling J (negative = antiferromagnetic)')
    ax.set_title(f'{name}: two coupled spins, real zeros of Z₂(h, h; J)')
    spin[name] = rows
RES['two_spin_real_zero_counts'] = spin
fig.suptitle('Spin count as a dimension: a second spin with J > 0 keeps the line, J < 0 loses it (stability passes up the chain, not down)')
fig.tight_layout()
fig.savefig(f'{OUT}/fig_spin_dimension.png', dpi=130)
plt.close(fig)


# ---------- 4. bath averaged away: two-level ----------
def bath2(k, z, f, a, t=0.0):
    return (1 - f) * k.H(z + 1j * a, t) + f * k.H(z - 1j * a, t)


fig, axs = plt.subplots(2, 4, figsize=(17, 7.5))
bathres = {}
for row, (k, name, x1, a) in enumerate([(PH, 'p-ball p=3, a=2', 40.0, 1.0), (R, 'Riemann Φ', 72.0, 4.0)]):
    for col, f in enumerate([0.5, 0.6, 0.75, 0.9]):
        ax = axs[row, col]
        xs = np.linspace(0.1, x1, 500); ys = np.linspace(-3, 3, 161)
        X, Y = np.meshgrid(xs, ys)
        F = bath2(k, X + 1j * Y, f, a)
        E = np.exp(np.mean(np.log(np.abs(F) + 1e-300), axis=0))[None, :]
        ax.imshow(domain_rgb(F, E), extent=[xs[0], xs[-1], ys[0], ys[-1]], origin='lower', aspect='auto')
        A = np.abs(F) / E
        m = (A < np.roll(A, 1, 0)) & (A < np.roll(A, -1, 0)) & (A < np.roll(A, 1, 1)) & (A < np.roll(A, -1, 1)) & (A < 0.3)
        m[0, :] = m[-1, :] = m[:, 0] = m[:, -1] = False
        zz = (X + 1j * Y)[m]
        ax.plot(zz.real, zz.imag, 'wo', ms=5, mec='k')
        ax.axhline(0, color='w', lw=0.6, ls='--')
        kprime = -1 / f
        ax.set_title(f'{name}\noccupation f = {f}, k′ = −1/f = {kprime:.2f}', fontsize=9)
        ax.set_xlabel('Re z'); ax.set_ylabel('Im z')
        bathres[f'{name} f={f}'] = {'zeros_found': len(zz), 'mean_Im': float(np.mean(zz.imag)) if len(zz) else None,
                                    'max_abs_Im': float(np.max(np.abs(zz.imag))) if len(zz) else None}
RES['bath2_averaged'] = bathres
fig.suptitle('Two-level bath averaged away: (1−f)·H(z + ia) + f·H(z − ia). Only half filling (f = 1/2, k′ = −2) keeps the zeros on the line;\n'
             'any other filling sends every zero into one half-plane (white dots). Colour is the phase, brightness the size.')
fig.tight_layout()
fig.savefig(f'{OUT}/fig_bath_two_level_averaged.png', dpi=120)
plt.close(fig)


# ---------- 5. bath kept: a second variable, the bath field h ----------
def bath_kept(k, x, h, a):
    g = k.H(x - 1j * a)              # H(x+ia) = conj(H(x-ia)) on the real line
    return 2 * np.real(np.exp(0.5j * h)[:, None] * g[None, :])


fig, axs = plt.subplots(1, 2, figsize=(13, 5))
for ax, (k, name, x1, a) in zip(axs, [(PH, 'p-ball p=3, a=2', 40.0, 1.0), (R, 'Riemann Φ', 72.0, 4.0)]):
    xs = np.linspace(0.05, x1, 2500)
    hs = np.linspace(-2 * np.pi, 2 * np.pi, 121)
    V = bath_kept(k, xs, hs, a)
    for j, h in enumerate(hs):
        v = V[j]
        i = np.where(np.sign(v[:-1]) * np.sign(v[1:]) < 0)[0]
        ax.plot(xs[i], h * np.ones(len(i)), 'k.', ms=1.8)
    ax.set_xlabel('Re z (real)'); ax.set_ylabel('bath field h (real)')
    ax.set_title(f'{name}: bath kept, F(z, h) = e^(ih/2) H(z − ia) + e^(−ih/2) H(z + ia)\nzeros for every real bath field stay real (each row is a full set)', fontsize=9)
    ax.axhline(0, color='tab:blue', lw=1)
    ax.text(0.3, 0.3, 'h = 0 is half filling, f = 1/2', color='tab:blue', fontsize=8)
fig.suptitle('Bath kept: F is stable in two variables. Averaging at filling f ≠ 1/2 means setting h = i·log((1−f)/f), a non-real field, off this plane')
fig.tight_layout()
fig.savefig(f'{OUT}/fig_bath_kept.png', dpi=130)
plt.close(fig)


# ---------- 6. three-level bath: weights (w, 1-2w, w) ----------
fig, axs = plt.subplots(2, 4, figsize=(17, 7.5))
three = {}
for row, (k, name, x1, a) in enumerate([(PH, 'p-ball p=3, a=2', 40.0, 2.0), (R, 'Riemann Φ', 72.0, 8.0)]):
    for col, w in enumerate([0.05, 0.15, 0.25, 0.35]):
        ax = axs[row, col]
        wt = lambda u, w=w, a=a: 1 - 2 * w + 2 * w * np.cosh(a * u)
        xs = np.linspace(0.1, x1, 500); ys = np.linspace(-3, 3, 161)
        X, Y = np.meshgrid(xs, ys)
        F = k.H(X + 1j * Y, weight=wt)
        E = np.exp(np.mean(np.log(np.abs(F) + 1e-300), axis=0))[None, :]
        ax.imshow(domain_rgb(F, E), extent=[xs[0], xs[-1], ys[0], ys[-1]], origin='lower', aspect='auto')
        zr = kern.real_zeros(k, x1, weight=wt)
        zc = kern.complex_zeros(k, x1, 3, weight=wt, nx=int(x1 * 10))
        off = zc[np.abs(zc.imag) > 1e-6]
        ax.plot(zr, 0 * zr, 'wo', ms=4, mec='k')
        ax.plot(off.real, off.imag, 'ro', ms=5, mec='k'); ax.plot(off.real, -off.imag, 'ro', ms=5, mec='k')
        ax.set_title(f'{name}\nw = {w}  ({"universal factor" if w >= 0.25 else "not a universal factor"})', fontsize=9)
        ax.set_xlabel('Re z'); ax.set_ylabel('Im z')
        three[f'{name} w={w}'] = {'real': len(zr), 'offline_upper': [[float(z.real), float(z.imag)] for z in off]}
RES['bath3_averaged'] = three
fig.suptitle('Three-level symmetric bath averaged away: kernel × (1 − 2w + 2w·cosh(au)). The factor is a Pólya universal factor only for w ≥ 1/4.\n'
             'White: real zeros. Red: zeros off the line in the window shown.')
fig.tight_layout()
fig.savefig(f'{OUT}/fig_bath_three_level.png', dpi=120)
plt.close(fig)

json.dump(RES, open('landscape.json', 'w'), indent=1)
print(json.dumps({k: v for k, v in RES.items() if not k.startswith('landscape_zeros')}, indent=1)[:6000])
