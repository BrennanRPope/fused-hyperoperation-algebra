import mpmath as mp, json
from epstein import Z
mp.mp.dps = 50
def m_of_x(x):
    q = mp.exp(-mp.pi*x); return (mp.jtheta(2, 0, q)/mp.jtheta(3, 0, q))**4
def amp_deg(m): return float(2*mp.asin(mp.sqrt(m))*180/mp.pi)
out = {}
# (1) theta(x) = sqrt(2K/pi) with x = K'/K
for m in (0.1, 0.5, 0.9):
    K = mp.ellipk(m); Kp = mp.ellipk(1-m); x = Kp/K
    th = mp.jtheta(3, 0, mp.exp(-mp.pi*x))
    print("m", m, "theta(x)", mp.nstr(th, 15), "sqrt(2K/pi)", mp.nstr(mp.sqrt(2*K/mp.pi), 15), "x", mp.nstr(x, 8))
# (2) zeta in pendulum form: pi^-s/2 Gamma(s/2) zeta(s) = int_0^1 (sqrt(2K/pi)-1)/2 (K'/K)^(s/2-1) * pi/(4 m(1-m) K^2) dm
# Legendre relation: dx/dm = -pi/(4 m (1-m) K^2), so the Mellin integral over x becomes an integral over swings m
for m in (0.2, 0.5, 0.8):
    xm = lambda mm: mp.ellipk(1-mm)/mp.ellipk(mm)
    print("dx/dm", mp.nstr(mp.diff(xm, m), 12), mp.nstr(-mp.pi/(4*m*(1-m)*mp.ellipk(m)**2), 12))
# (3) real-zero threshold: Z(1/2, y) changes sign
f = lambda y: Z(mp.mpf(0.5)+mp.mpf(10)**-25, y, K=30).real
for y in (5, 10, 30, 49, 50): print("Z(1/2,", y, ") =", mp.nstr(f(y), 8))
ystar = mp.findroot(f, (49, 50), solver="anderson", tol=1e-30)
print("threshold y* =", mp.nstr(ystar, 10), " x* =", mp.nstr(mp.sqrt(ystar), 8), " m* =", mp.nstr(m_of_x(mp.sqrt(ystar)), 8), "amp", amp_deg(m_of_x(mp.sqrt(ystar))))
out["ystar"] = float(ystar); out["amp_star"] = amp_deg(m_of_x(mp.sqrt(ystar)))
for y in (1, 2, 3, 4, 7): print("class-number-one swing y =", y, "amp", round(amp_deg(m_of_x(mp.sqrt(y))), 3), "deg")
json.dump(out, open("pend.json", "w"))
