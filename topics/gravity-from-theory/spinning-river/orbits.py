"""DRAFT: two spinning bodies in orbit (round 12). Shared by checks.py (section O) and figures.py.

Units G = c = 1, total mass M = 1 unless stated. Run `python3 orbits.py` to make the long spin ensemble
(orbits.json; about 20 minutes on 4 cores): 80 inspirals of mass ratio 0.7 from r = 1000 M to 10 M."""
import json
import multiprocessing
import os

import numpy as np
from scipy.integrate import solve_ivp

HERE = os.path.dirname(os.path.abspath(__file__))


# ---------------------------------------------------------------- planar orbit: slow-orbit (1PN) motion + waves
def orbit_rhs(t, y, nu, pn=True, rr=True):
    """Relative motion of two bodies (harmonic gauge): a = -(M/r^2)[(1 + A) n + B v].
    1PN part (Einstein-Infeld-Hoffmann) and 2.5PN radiation reaction (the waves' back-push)."""
    x, v = y[:2], y[2:]
    r = np.hypot(*x)
    n = x / r
    rd = n @ v
    v2 = v @ v
    A = B = 0.0
    if pn:
        A += (1 + 3 * nu) * v2 - 1.5 * nu * rd**2 - (4 + 2 * nu) / r
        B += -(4 - 2 * nu) * rd
    if rr:
        A += -(8 / 5) * nu / r * rd * (17 / (3 * r) + 3 * v2)
        B += (8 / 5) * nu / r * (3 / r + v2)
    a = -(1 / r**2) * ((1 + A) * n + B * v)
    return [v[0], v[1], a[0], a[1]]


def orbit(nu, a0, e0, T, pn=True, rr=True, rstop=3.0, max_step=2.0):
    rp = a0 * (1 - e0)
    vp = np.sqrt((1 + e0) / (a0 * (1 - e0)))
    stop = lambda t, y, *_: np.hypot(y[0], y[1]) - rstop  # noqa: E731
    stop.terminal = True
    return solve_ivp(orbit_rhs, (0, T), [rp, 0, 0, vp], args=(nu, pn, rr), rtol=1e-11, atol=1e-13,
                     dense_output=True, events=stop, max_step=max_step)


def circular(nu, r0, T):
    v0 = r0 * np.sqrt((1 - (3 - nu) / r0) / r0**3)
    return solve_ivp(orbit_rhs, (0, T), [r0, 0, 0, v0], args=(nu, True, True), rtol=1e-11, atol=1e-13,
                     dense_output=True, max_step=2.0)


def energy_1pn(y, nu):
    x, v = y[:2], y[2:]
    r = np.hypot(*x)
    v2 = v @ v
    rd = (x @ v) / r
    return (0.5 * v2 - 1 / r + 3 / 8 * v2**2 - 9 / 8 * nu * v2**2
            + (0.5 * rd**2 * nu + 1.5 * v2 + 0.5 * nu * v2) / r + 0.5 / r**2)


def periastron_advance(nu, a0, e0, norb=6):
    """Mean turn of the closest-approach point per orbit, against 6 pi M / p."""
    s = orbit(nu, a0, e0, norb * 2 * np.pi * a0**1.5 * 1.3, rr=False)
    t = np.linspace(0, s.t[-1], 600000)
    Y = s.sol(t)
    r = np.hypot(Y[0], Y[1])
    idx = [i for i in range(1, len(r) - 1) if r[i] < r[i - 1] and r[i] <= r[i + 1]]
    ang = np.unwrap(np.arctan2(Y[1][idx], Y[0][idx]))
    p = 2 * r.min() * r.max() / (r.min() + r.max())
    return float(np.diff(np.r_[0.0, ang]).mean()), float(6 * np.pi / p)


def flux_ratio(nu, r0, norb=6):
    """Energy the orbit loses per time, over the quadrupole formula's (32/5) nu^2 M^5 / r^5."""
    P = 2 * np.pi * r0**1.5
    s = circular(nu, r0, norb * P)
    ts = np.linspace(0, norb * P, 61)
    Es = np.array([energy_1pn(s.sol(t), nu) for t in ts])
    rs = np.array([np.hypot(*s.sol(t)[:2]) for t in ts])
    return float(np.polyfit(ts, Es, 1)[0] * nu / (-(32 / 5) * nu**2 / np.mean(rs)**5))


# ---------------------------------------------------------------- spins: orbit-averaged precession + shrinking orbit
def masses(q, chi1, chi2):
    m1, m2 = 1 / (1 + q), q / (1 + q)
    return m1, m2, m1 * m2, chi1 * m1**2, chi2 * m2**2


def spin_rhs(r, y, m1, m2, nu):
    """d/dr of (S1, S2, L-hat). Each spin turns about Omega_k (check O2, orbit-averaged); the orbit's direction turns
    so that the total J = L + S1 + S2 stays fixed; the orbit shrinks by the quadrupole rate dr/dt = -(64/5) nu / r^3."""
    S1, S2, Lh = y[:3], y[3:6], y[6:9] / np.linalg.norm(y[6:9])
    L = nu * np.sqrt(r) * Lh
    O1 = ((2 + 1.5 * m2 / m1) * L + 0.5 * S2 - 1.5 * (Lh @ S2) * Lh) / r**3
    O2 = ((2 + 1.5 * m1 / m2) * L + 0.5 * S1 - 1.5 * (Lh @ S1) * Lh) / r**3
    dS1, dS2 = np.cross(O1, S1), np.cross(O2, S2)
    dLh = -(dS1 + dS2) / np.linalg.norm(L)
    return np.concatenate([dS1, dS2, dLh]) / (-(64 / 5) * nu / r**3)


def inspiral_spins(q, chi1, chi2, s1hat, s2hat, r0=300.0, r1=10.0, npts=2000):
    m1, m2, nu, s1, s2 = masses(q, chi1, chi2)
    y0 = np.concatenate([s1 * np.asarray(s1hat), s2 * np.asarray(s2hat), [0, 0, 1.0]])
    rs = np.geomspace(r0, r1, npts)
    return solve_ivp(spin_rhs, (r0, r1), y0, args=(m1, m2, nu), t_eval=rs, rtol=1e-9, atol=1e-12, method="DOP853")


def spin_angles(y):
    """(angle between the spins, spin 1 to orbit, spin 2 to orbit, angle between the spins' shadows on the orbit
    plane), in degrees."""
    S1, S2, Lh = y[:3], y[3:6], y[6:9] / np.linalg.norm(y[6:9])
    c = lambda a, b: np.clip(a @ b / np.linalg.norm(a) / np.linalg.norm(b), -1, 1)  # noqa: E731
    p1, p2 = S1 - (S1 @ Lh) * Lh, S2 - (S2 @ Lh) * Lh
    return [float(np.degrees(np.arccos(x))) for x in (c(S1, S2), c(S1, Lh), c(S2, Lh), c(p1, p2))]


def unit(th_deg, ph_deg):
    th, ph = np.radians(th_deg), np.radians(ph_deg)
    return np.array([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)])


def _one(args):
    q, seed, r0 = args
    rng = np.random.default_rng(seed)
    iso = lambda: (lambda v: v / np.linalg.norm(v))(rng.normal(size=3))  # noqa: E731
    s = inspiral_spins(q, 0.8, 0.8, iso(), iso(), r0=r0, r1=10.0, npts=3)
    return spin_angles(s.y[:, 0]), spin_angles(s.y[:, -1])


def ensemble(q, n, r0, seed0=0, procs=None):
    args = [(q, seed0 + k, r0) for k in range(n)]
    if procs == 1:
        res = [_one(a) for a in args]
    else:
        with multiprocessing.Pool(procs or multiprocessing.cpu_count()) as P:
            res = P.map(_one, args)
    return np.array([r[0] for r in res]), np.array([r[1] for r in res])


if __name__ == "__main__":
    a0, a1 = ensemble(0.7, 80, 1000.0, seed0=1000)
    json.dump({"q": 0.7, "chi": 0.8, "r0": 1000.0, "r1": 10.0, "start": a0.tolist(), "end": a1.tolist()},
              open(os.path.join(HERE, "orbits.json"), "w"))
    print("orbits.json written")
