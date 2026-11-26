"""Polish for the elliptic-modulus page: independent high-precision checks. Run: python3 polish.py > polish_output.txt"""
import mpmath as mp, numpy as np
from heights import height
mp.mp.dps = 50
def G(x, m): return mp.ellipfun('sn', x, m=m)*mp.ellipfun('dn', x, m=m)/mp.ellipfun('cn', x, m=m)
out = []
# P1: the fingerprint read off an orbit, independently of the series. For an attracting odd germ the orbit obeys
# n = -1/(2 a x_n^2) + R log x_n + C + o(1) (Fatou coordinate of a contact-3 germ). Fit R from the orbit of G at m = 0.8 (mpmath, 50 digits).
m = mp.mpf('0.8'); a = (1 - 2*m)/3; Rexact = mp.mpf(6)/5 - mp.mpf(9)/(10*(1 - 2*m)**2)
x = mp.mpf('0.5'); pts = []
for n in range(1, 200001):
    x = G(x, m) if n < 2000 else x + a*x**3 + 2*(1 - m + m*m)/15*x**5 + (-4*m**3/315 + 2*m**2/105 - 4*m/35 + mp.mpf(17)/315)*x**7
    if n in (50000, 100000, 200000): pts.append((n, x))
# solve for (R, C) from two late points, using the next asymptotic order to cancel o(1) by Richardson
def RC(p, q):
    (n1, x1), (n2, x2) = p, q
    A = mp.matrix([[mp.log(x1), 1], [mp.log(x2), 1]]); b = mp.matrix([n1 + 1/(2*a*x1**2), n2 + 1/(2*a*x2**2)])
    return mp.lu_solve(A, b)[0]
R1, R2 = RC(pts[0], pts[1]), RC(pts[1], pts[2])
out.append("P1 fingerprint from the orbit at m = 0.8: %s, %s (later pair) vs formula %s" % (mp.nstr(R1, 8), mp.nstr(R2, 8), mp.nstr(Rexact, 8)))
out.append("   (orbit uses G itself for 2000 steps, then its Taylor polynomial to order 7, accurate below |x| ~ 0.03)")
# P2: the wall distance is the law's logarithm out to t = infinity: K = int_0^inf dt / sqrt(1 + 2 kappa t^2 + t^4)
for mm in [mp.mpf('0.3'), mp.mpf('0.5'), mp.mpf('0.9')]:
    k = 1 - 2*mm; I = mp.quad(lambda t: 1/mp.sqrt(1 + 2*k*t**2 + t**4), [0, 1, mp.inf])
    out.append("P2 m = %s: int_0^inf dt/sqrt(1 + 2kt^2 + t^4) = %s, K(m) = %s, diff %s" % (mm, mp.nstr(I, 20), mp.nstr(mp.ellipk(mm), 20), mp.nstr(abs(I - mp.ellipk(mm)), 3)))
# P3: half height at 30 digits by direct regular iteration in mpmath vs the double-precision routine
mm = mp.mpf('0.8'); a = (1 - 2*mm)/3
import numseries as NS
NS.mp.mp.dps = 50
V = NS.julia(NS.slope(mm, 31), 31); Ps = NS.flow(V, mp.mpf('0.5'), 31)
def P(y): return sum(c*y**j for j, c in enumerate(Ps))
def Ginv(t): return mp.findroot(lambda y: G(y, mm) - t, t)
x0 = mp.mpf('0.6'); y = x0; n = 0
while abs(y) > mp.mpf('0.02'): y = G(y, mm); n += 1
z = P(y)
for _ in range(n): z = Ginv(z)
dbl = height(np.array([0.6]), 0.8, 0.5)[0]
out.append("P3 G_0.8^(1/2)(0.6) at 50 digits: %s (%d steps); double routine %r; diff %s" % (mp.nstr(z, 30), n, dbl, mp.nstr(abs(z - dbl), 3)))
zz = z; y = zz; n = 0
while abs(y) > mp.mpf('0.02'): y = G(y, mm); n += 1
w = P(y)
for _ in range(n): w = Ginv(w)
out.append("   half of the half: %s vs G(0.6) = %s, diff %s" % (mp.nstr(w, 30), mp.nstr(G(x0, mm), 30), mp.nstr(abs(w - G(x0, mm)), 3)))
print("\n".join(out))
