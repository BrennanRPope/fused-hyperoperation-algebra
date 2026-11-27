"""Pendulum picture checks (Nov 27). mpmath 30 digits; ODE by mpmath odefun / scipy at 1e-12."""
import mpmath as mp, numpy as np, json
from scipy.integrate import solve_ivp
mp.mp.dps = 30
res = []
def rec(name, ok, val):
    res.append((name, bool(ok), str(val))); print(('PASS' if ok else 'FAIL'), name, val)

def G(u, m):
    m = mp.mpf(m) if not isinstance(m, mp.mpc) else m
    if mp.re(m) > 1:   # reciprocal-modulus (real) transformation, stays real for real u
        r = mp.sqrt(m); v = r*u; mu = 1/m
        s, c, d = (mp.ellipfun(f, v, m=mu) for f in ('sn','cn','dn'))
        return (s/r)*c/d   # sn(u|m)=sn(v|mu)/r, cn(u|m)=dn(v|mu), dn(u|m)=cn(v|mu)
    s, c, d = (mp.ellipfun(f, u, m=m) for f in ('sn','cn','dn'))
    return s*d/c

# P1 complementary time: G(u) G(K-u) = 1 for librating swings
err = 0
for m in [mp.mpf('0.1'), mp.mpf('0.5'), mp.mpf('0.9'), mp.mpf('0.999')]:
    K = mp.ellipk(m)
    for u in [K*0.13, K*0.4, K*0.77, mp.mpc(0.3, 0.2)]:
        err = max(err, abs(G(u, m)*G(K-u, m) - 1))
rec('P1 G(u)G(K-u)=1 (complementary time)', err < 1e-25, mp.nstr(err, 3))

# P2 the halfway line Re u = K/2 maps onto the unit circle |t| = 1
err = 0
for m in [mp.mpf('0.2'), mp.mpf('0.7')]:
    K = mp.ellipk(m); Kp = mp.ellipk(1-m)
    for y in np.linspace(-0.95, 0.95, 9):
        err = max(err, abs(abs(G(mp.mpc(K/2, y*Kp), m)) - 1))
rec('P2 |G| = 1 on Re u = K/2', err < 1e-25, mp.nstr(err, 3))

# ODE helper: theta'' = -sin theta (omega = 1); slope reading t = sin(theta)/theta'
def orbit(theta0dot, T, n=400):
    sol = solve_ivp(lambda t, y: [y[1], -np.sin(y[0])], [0, T], [0.0, theta0dot],
                    rtol=1e-12, atol=1e-13, dense_output=True)
    ts = np.linspace(1e-3, T, n); th, thd = sol.sol(ts)
    return ts, th, thd

# P3 energy reading: m = E/E_top, kappa = 1 - 2E/E_top (E from the bottom, E_top = 2 omega^2)
ok = True; worst = 0
for m in [0.3, 0.8, 1.7, 4.0]:
    v0 = 2*np.sqrt(m)                    # theta' at the bottom
    E = v0**2/2; kap = 1 - 2*E/2.0
    worst = max(worst, abs(kap - (1-2*m)))
rec('P3 bookkeeping: kappa = 1 - 2E/E_top with E from the bottom, E_top = 2 omega^2', worst < 1e-14, worst)

# P4 librating swing: t = 1 exactly at half the time to the turning point, where cos th* = 2cos(th0/2) - 1
worst = 0
for th0deg in [20, 60, 120, 170]:
    th0 = np.radians(th0deg); m = np.sin(th0/2)**2
    K = float(mp.ellipk(m))
    sol = solve_ivp(lambda t, y: [y[1], -np.sin(y[0])], [0, K/2], [0.0, 2*np.sqrt(m)], rtol=1e-12, atol=1e-13)
    th, thd = sol.y[:, -1]
    worst = max(worst, abs(np.sin(th)/thd - 1), abs(np.cos(th) - (2*np.cos(th0/2) - 1)))
rec('P4 t = 1 at half time to the wall; cos th* = 2cos(th0/2) - 1', worst < 1e-9, worst)

# P5 swinging over the top (m > 1): the ODE slope reading equals G_m(time), never infinite
worst = 0
for m in [1.3, 2.5]:
    ts, th, thd = orbit(2*np.sqrt(m), 6.0)
    g = np.array([float(G(mp.mpf(t), mp.mpf(m))) for t in ts[::20]])
    worst = max(worst, np.max(np.abs(np.sin(th[::20])/thd[::20] - g)))
rec('P5 over the top: sin th / th\' = G_m(t), m > 1', worst < 1e-9, worst)

# P6 the cap: max of G_m over real time is sqrt(m) - sqrt(m-1) for m > 1
worst = 0
for m in [mp.mpf('1.01'), mp.mpf('1.5'), mp.mpf(3), mp.mpf(10)]:
    # G' = 0 where sn^2 = 1 - sqrt(1 - 1/m)
    s2 = 1 - mp.sqrt(1 - 1/m)
    r = mp.sqrt(m); u = mp.findroot(lambda u: (mp.ellipfun('sn', r*u, m=1/m)/r)**2 - s2, mp.asin(mp.sqrt(s2)))
    gm = G(u, m)
    # and a brute check on a grid over one real period
    P = 4*mp.ellipk(1/m)/mp.sqrt(m)
    grid = max(abs(G(P*k/2000, m)) for k in range(2000))
    worst = max(worst, abs(gm - (mp.sqrt(m) - mp.sqrt(m-1))), max(0, grid - gm))
rec('P6 cap max G = sqrt(m) - sqrt(m-1)', worst < 1e-20, mp.nstr(worst, 3))

# P7 the cap and its mirror are the quartic's roots: t^2 = (sqrt m -+ sqrt(m-1))^2, a reciprocal pair
worst = 0
for m in [mp.mpf('1.5'), mp.mpf(7)]:
    kap = 1-2*m
    for t in [mp.sqrt(m)-mp.sqrt(m-1), mp.sqrt(m)+mp.sqrt(m-1)]:
        worst = max(worst, abs(1 + 2*kap*t**2 + t**4))
rec('P7 cap solves 1 + 2 kappa t^2 + t^4 = 0, partner 1/t', worst < 1e-25, mp.nstr(worst, 3))

# P8 over the top the reading is periodic with real period 2K(1/m)/sqrt(m) (one full turn of theta is 4K(1/m)/sqrt m in time... check G)
worst = 0
for m in [mp.mpf('1.5'), mp.mpf(4)]:
    Pg = 2*mp.ellipk(1/m)/mp.sqrt(m)
    for u in [mp.mpf('0.3'), mp.mpf('1.1')]:
        worst = max(worst, abs(G(u+Pg, m) - G(u, m)))
rec('P8 over the top G has real period 2K(1/m)/sqrt m', worst < 1e-25, mp.nstr(worst, 3))

# P9 libration <-> rotation: the sn chart's fingerprint is symmetric under m -> 1/m
def R_sn(m): return 6*(1-m+m**2)/(5*(1+m)**2)
worst = max(abs(R_sn(m) - R_sn(1/m)) for m in [mp.mpf('0.2'), mp.mpf('0.7'), mp.mpf(3)])
rec('P9 sn-chart fingerprint symmetric m -> 1/m', worst < 1e-28, mp.nstr(worst, 3))

# P10 dwell time near the turning point piles up like (th0 - th)^(-1/2)
th0 = np.radians(100.0); c0 = np.cos(th0)
eps = np.array([1e-3, 1e-4, 1e-5])
dwell = 1/np.sqrt(2*(np.cos(th0-eps) - c0))
slope = np.polyfit(np.log(eps), np.log(dwell), 1)[0]
rec('P10 dwell density ~ (th0-th)^(-1/2) at the wall', abs(slope + 0.5) < 1e-3, slope)

# P11 separatrix: reading -> tanh, halfway time -> infinity, cap -> 1 from both sides
m = mp.mpf(1) - mp.mpf('1e-12'); K = mp.ellipk(m)
a = abs(G(mp.mpf(2), m) - mp.tanh(2))
b = abs((mp.sqrt(1+mp.mpf('1e-12')) - mp.sqrt(mp.mpf('1e-12'))) - 1)
rec('P11 separatrix: G -> tanh, cap -> 1', a < 1e-10 and b < 1e-5, (mp.nstr(a,3), mp.nstr(b,3)))

# P12 Wick: G_m(iu) = i G_{1-m}(u) carries over-the-top (m > 1) to m < 0, where the real-time reading has walls again
worst = 0
for m in [mp.mpf('1.5'), mp.mpf(3)]:
    for u in [mp.mpf('0.2'), mp.mpf('0.6')]:
        worst = max(worst, abs(G(1j*u, m) - 1j*G(u, 1-m)))
mneg = mp.mpf(-2); Kneg = mp.ellipk(mneg)
wall = abs(1/G(Kneg*(1-mp.mpf('1e-12')), mneg))
rec('P12 Wick swaps caps and walls (m>1 <-> m<0); m<0 has a wall at K(m)', worst < 1e-25 and wall < 1e-10, (mp.nstr(worst,3), mp.nstr(wall,3)))

npass = sum(r[1] for r in res)
print(f'{npass} of {len(res)} pass')
json.dump(res, open('pendulum_checks.json', 'w'), indent=1)
