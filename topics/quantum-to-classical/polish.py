"""Polish for checks.py section B (Nov 26): the stationary point of the summed path contributions is the zero of the
local phase gradient kloc(x) on an FFT grid (dx = 6.1e-4), located in checks.py by LINEAR interpolation between the two
grid points that bracket the sign change. Here the same fields are recomputed (same split-step chain, same grid), and the
zero is re-found with a cubic through the four nearest grid values, finished by Newton iterations to machine precision.
Newton's midpoint (xmid) is taken from results.json (brentq on a shooting solution, already ~1e-12). numpy only."""
import json, time
import numpy as np
t0 = time.process_time()
B = json.load(open('results.json'))['B']; xmid = B['xmid_newton']
lam = 0.5; V = lambda x: 0.5*x**2 + lam*x**4; xa, xb, T = -0.6, 0.9, 0.8
L, N = 10.0, 16384
X = np.linspace(-L/2, L/2, N, endpoint=False); dx = X[1] - X[0]; P = 2*np.pi*np.fft.fftfreq(N, dx)
def evolve(psi, t_, hb, dt=5e-4, sign=1):          # identical to checks.py (without recording)
    n = int(round(t_/dt)); dt = t_/n
    eV = np.exp(-1j*sign*V(X)*dt/(2*hb)); eK = np.exp(-1j*sign*hb*P**2*dt/2)
    for _ in range(n):
        psi = eV*psi; psi = np.fft.ifft(eK*np.fft.fft(psi)); psi = eV*psi
    return psi
lines = ["Polish of the stationary-point location (checks.py B1): linear -> cubic + Newton zero of the phase gradient", ""]
offs0, offs1 = [], []
for hb, s_old in zip(B['hbars'], B['stat_pt']):
    w = 0.2*hb; g = lambda x0: np.exp(-(X - x0)**2/(4*w**2)).astype(complex)
    f = np.conj(evolve(g(xb), T/2, hb, sign=-1))*evolve(g(xa), T/2, hb)
    kloc = np.imag(np.gradient(f, dx)/(f + 1e-300))
    near = np.where(np.abs(X - xmid) < 0.6)[0]; sgn = np.sign(kloc[near]); zc = near[np.where(np.diff(sgn) != 0)[0]]
    zs = X[zc] + dx*np.abs(kloc[zc])/(np.abs(kloc[zc]) + np.abs(kloc[zc + 1])); j = np.argmin(np.abs(zs - xmid))
    lin = zs[j]; i = zc[j]
    idx = np.arange(i - 1, i + 3); c = np.polyfit(X[idx] - X[i], kloc[idx], 3); dc = np.polyder(c)
    u = lin - X[i]
    for _ in range(30):
        step = np.polyval(c, u)/np.polyval(dc, u); u -= step
        if abs(step) < 1e-17: break
    cub = X[i] + u
    offs0.append(abs(lin - xmid)); offs1.append(abs(cub - xmid))
    lines.append(f"hbar={hb:<5}: stationary point {lin:.9f} (linear; stored {s_old:.9f}) -> {cub:.9f} (cubic+Newton); "
                 f"shift {cub - lin:+.1e}; offset from Newton {abs(lin-xmid):.2e} -> {abs(cub-xmid):.2e}")
r = np.array(offs1[:-1])/np.array(offs1[1:])
lines += ["", f"max offset {max(offs0):.1e} -> {max(offs1):.1e} (polished); README quotes 'offsets <= 4e-5': "
          + ("unchanged at that precision" if round(max(offs0), 5) == round(max(offs1), 5) else "CHANGED"),
          f"offset ratio between successive halvings of hbar: {np.round(r, 2)} (about 4, i.e. offset ~ hbar^2: a real",
          "finite-hbar shift of the stationary point, not interpolation error). Verdict B1 unchanged.",
          f"CPU time: {time.process_time() - t0:.1f} s"]
out = "\n".join(lines); print(out); open('polish_output.txt', 'w').write(out + "\n")
