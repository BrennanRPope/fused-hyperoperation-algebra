"""Polish step for the one searched number in checks.py: the force of peak muscle power.

checks.py finds it with scipy's bounded Brent search (minimize_scalar, default xatol 1e-5) and
compares it with the closed form F/F0 = (a/F0)(sqrt(1 + F0/a) - 1) to 1e-5. The numbers on the
page (0.309, 0.309, 0.095) come from the closed form, not the search. The polish starts from the
search answer and runs Newton on dP/dF = 0. scipy is not installed here, so the search answer is
reproduced with an exact port of scipy's bounded Brent routine.
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

F0 = 1.0; a = 0.25*F0; vmax = 1.0; b = a*vmax/F0
P = lambda f: f*b*(F0-f)/(f+a)
fopt = a*(np.sqrt(1+F0/a)-1)
opt0 = fminbound(lambda f: -P(f), 0, F0)            # minimize_scalar(..., method="bounded") default xatol 1e-5
t = time.perf_counter()
f = opt0
for it in range(8):          # Newton on P'(f); P' = b[(F0 - 2f)(f + a) - f(F0 - f)]/(f + a)^2, numerator N(f)
    N = (F0 - 2*f)*(f + a) - f*(F0 - f)      # = -f^2 - 2af + aF0
    Np = -2*f - 2*a
    step = N/Np; f -= step
    if abs(step) < 1e-16: break
dt = time.perf_counter() - t
print("Peak power of Hill's muscle law (a/F0 = 0.25)")
print(f"  closed form            F/F0 = {fopt:.16f}")
print(f"  before (bounded Brent) F/F0 = {opt0:.16f}   off by {opt0-fopt:+.1e}")
print(f"  after  (Newton, {it+1} steps) F/F0 = {f:.16f}   off by {f-fopt:+.1e}")
gain = math.log10(abs(opt0-fopt)/max(abs(f-fopt), 2.2e-16*fopt))
print(f"  digits gained: about {gain:.0f}; wall time {dt*1e3:.3f} ms")
print(f"  page numbers (from the closed form): F/F0 = {fopt:.3f}, v/vmax = {b*(F0-fopt)/(fopt+a)/vmax:.3f}, power {P(fopt)/(F0*vmax):.3f} -> no change")
