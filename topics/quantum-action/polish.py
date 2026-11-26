"""Polish step for the one searched number in checks.py (B4, Fermat/Snell).

checks.py finds the least-optical-path crossing point with scipy's bounded Brent search
(minimize_scalar, xatol=1e-12). A minimiser of a smooth function only pins x to about
sqrt(machine eps) in practice, so the Snell check there agreed to ~3e-8. The polish starts
from that answer and runs Newton on the stationarity condition dOPL/dx = n1 sin a1 - n2 sin a2 = 0.
scipy is not installed on this machine, so the 'before' value comes from an exact line-by-line
port of scipy's bounded Brent routine (it reproduces checks_output.txt to all 10 printed digits).
Run: python3 polish.py > polish_output.txt
"""
import math, time
import numpy as np

def fminbound(func, x1, x2, xatol=1e-5, maxfun=500):
    """Port of scipy.optimize._minimize_scalar_bounded (bounded Brent)."""
    sqrt_eps = np.sqrt(2.2e-16); golden_mean = 0.5 * (3.0 - np.sqrt(5.0))
    a, b = x1, x2
    fulc = a + golden_mean * (b - a); nfc, xf = fulc, fulc
    rat = e = 0.0; x = xf; fx = func(x); num = 1
    ffulc = fnfc = fx; xm = 0.5 * (a + b)
    tol1 = sqrt_eps * np.abs(xf) + xatol / 3.0; tol2 = 2.0 * tol1
    while np.abs(xf - xm) > (tol2 - 0.5 * (b - a)):
        golden = 1
        if np.abs(e) > tol1:
            golden = 0
            r = (xf - nfc) * (fx - ffulc); q = (xf - fulc) * (fx - fnfc)
            p = (xf - fulc) * q - (xf - nfc) * r; q = 2.0 * (q - r)
            if q > 0.0: p = -p
            q = np.abs(q); r = e; e = rat
            if (np.abs(p) < np.abs(0.5 * q * r)) and (p > q * (a - xf)) and (p < q * (b - xf)):
                rat = (p + 0.0) / q; x = xf + rat
                if ((x - a) < tol2) or ((b - x) < tol2):
                    si = np.sign(xm - xf) + ((xm - xf) == 0); rat = tol1 * si
            else:
                golden = 1
        if golden:
            e = (a - xf) if xf >= xm else (b - xf)
            rat = golden_mean * e
        si = np.sign(rat) + (rat == 0)
        x = xf + si * np.maximum(np.abs(rat), tol1)
        fu = func(x); num += 1
        if fu <= fx:
            if x >= xf: a = xf
            else: b = xf
            fulc, ffulc = nfc, fnfc; nfc, fnfc = xf, fx; xf, fx = x, fu
        else:
            if x < xf: a = x
            else: b = x
            if (fu <= fnfc) or (nfc == xf):
                fulc, ffulc = nfc, fnfc; nfc, fnfc = x, fu
            elif (fu <= ffulc) or (fulc == xf) or (fulc == nfc):
                fulc, ffulc = x, fu
        xm = 0.5 * (a + b); tol1 = sqrt_eps * np.abs(xf) + xatol / 3.0; tol2 = 2.0 * tol1
        if num >= maxfun: break
    return xf

n1, n2 = 1.0, 1.5; Apt, Bpt = (0, 1), (1, -1)
opl = lambda xi: n1*math.hypot(xi-Apt[0], Apt[1]) + n2*math.hypot(Bpt[0]-xi, Bpt[1])
sines = lambda xi: ((xi-Apt[0])/math.hypot(xi-Apt[0], Apt[1]), (Bpt[0]-xi)/math.hypot(Bpt[0]-xi, Bpt[1]))

xi0 = fminbound(opl, 0, 1, xatol=1e-12)
s1, s2 = sines(xi0); gap0 = n1*s1 - n2*s2
print("B4 Fermat/Snell: crossing point of the least optical path")
print(f"  before (bounded Brent, xatol 1e-12): x = {xi0:.16f}   n1 sin a1 = {n1*s1:.12f}, n2 sin a2 = {n2*s2:.12f}, gap {gap0:.2e}")

t = time.perf_counter()
xi = xi0
for it in range(6):          # Newton on g(x) = dOPL/dx = n1 sin a1 - n2 sin a2
    d1 = math.hypot(xi-Apt[0], Apt[1]); d2 = math.hypot(Bpt[0]-xi, Bpt[1])
    g = n1*(xi-Apt[0])/d1 - n2*(Bpt[0]-xi)/d2
    gp = n1*Apt[1]**2/d1**3 + n2*Bpt[1]**2/d2**3
    step = g/gp; xi -= step
    if abs(step) < 1e-17: break
dt = time.perf_counter() - t
s1, s2 = sines(xi); gap = n1*s1 - n2*s2
print(f"  after  (Newton on dOPL/dx, {it+1} steps): x = {xi:.16f}   n1 sin a1 = {n1*s1:.12f}, n2 sin a2 = {n2*s2:.12f}, gap {gap:.2e}")
dig = math.log10(abs(xi0-xi)/max(abs(xi)*2.2e-16, 1e-300)) if xi0 != xi else 0
print(f"  change in x: {xi-xi0:+.2e}; Snell gap {abs(gap0):.1e} -> {abs(gap):.1e} (about {dig:.0f} more correct digits in x)")
print(f"  wall time of the polish: {dt*1e3:.3f} ms")
print("  verdict unchanged: Snell's law holds at the least optical path (the check's 1e-6 tolerance was already met).")
