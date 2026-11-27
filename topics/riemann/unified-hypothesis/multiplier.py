"""Is Phi a Laguerre-multiplier child of Polya's Phi* (1926)?

G_K(z) = int_0^inf K(u) cos(z u) du = sum_k (-1)^k m_{2k}(K) z^{2k}/(2k)!, m_{2k} = int_0^inf u^{2k} K.
Write r_k = m_{2k}(Phi) / m_{2k}(Phi*). If {r_k} were a multiplier sequence (Polya-Schur 1914:
sum r_k x^k/k! in Laguerre-Polya class with zeros of one sign), Laguerre's theorem would carry
Polya's real zeros for Phi* over to Phi, i.e. RH.
Necessary conditions tested: Turan r_k^2 >= r_{k-1} r_{k+1}; all Jensen polynomials
J_n(x) = sum_j C(n,j) r_j x^j real-rooted (for every n this is equivalent to the
multiplier property for a positive sequence, Polya-Schur).
Also run on the p-ball / Polya pair, where the ratio is Gamma(a+1)/Gamma((2k+1)/p + a + 1),
to check the test passes where it should.
"""
import json
import mpmath as mp

mp.mp.dps = 40
KMAX = 40


def Phi(u):
    return mp.fsum((2 * mp.pi**2 * n**4 * mp.e**(9 * u) - 3 * mp.pi * n**2 * mp.e**(5 * u)) * mp.e**(-mp.pi * n**2 * mp.e**(4 * u)) for n in range(1, 9))


def PhiStar(u):   # Polya 1926, in the Rodgers-Tao variable
    return 2 * mp.pi**2 * (mp.e**(9 * u) + mp.e**(-9 * u)) * mp.e**(-mp.pi * (mp.e**(4 * u) + mp.e**(-4 * u)))


def moments(K):
    return [mp.quad(lambda u: u**(2 * k) * K(u), [0, 0.25, 0.5, 0.75, 1, 1.5]) for k in range(KMAX + 1)]


def jensen_check(r, nmax=30):
    bad = []
    for n in range(2, nmax + 1):
        c = [mp.binomial(n, j) * r[j] for j in range(n + 1)]
        roots = mp.polyroots(c[::-1], maxsteps=400, extraprec=400)
        im = max(abs(mp.im(x)) for x in roots)
        if im > mp.mpf(10)**(-20) * max(1, max(abs(x) for x in roots)):
            bad.append((n, float(im)))
    return bad


def report(name, r):
    tur = [k for k in range(1, len(r) - 1) if r[k]**2 - r[k - 1] * r[k + 1] < 0]
    jb = jensen_check(r)
    print(f'{name}: r_0..r_6 = {[mp.nstr(x, 8) for x in r[:7]]}  r_40 = {mp.nstr(r[-1], 10)}', flush=True)
    print(f'   Turan fails at k = {tur[:12]}{"..." if len(tur) > 12 else ""}  ({len(tur)} of {len(r)-2})', flush=True)
    print(f'   first non-real Jensen polynomials (degree, max |Im|): {jb[:6]}', flush=True)
    return {'r': [mp.nstr(x, 15) for x in r], 'turan_fail': tur, 'jensen_fail': jb}


if __name__ == '__main__':
    out = {}
    mP = moments(Phi)
    mS = moments(PhiStar)
    # sanity (float64): Polya's Phi* transform vs Phi, real zeros below z = 72
    import numpy as np, kern
    from scipy.special import roots_legendre
    x, w = roots_legendre(500); u = (x + 1) * 0.95 / 2; w = w * 0.95 / 2
    ps = 2*np.pi**2*(np.exp(9*u)+np.exp(-9*u))*np.exp(-np.pi*(np.exp(4*u)+np.exp(-4*u)))
    KS = kern.Kernel(u, w * ps, 'Phi*')
    out['PhiStar_real_zeros_below_72'] = [float(v) for v in kern.real_zeros(KS, 72)]
    out['Phi_real_zeros_below_72'] = [float(v) for v in kern.real_zeros(kern.riemann(), 72)]
    print('Phi* zeros', np.round(out['PhiStar_real_zeros_below_72'], 2), ' Phi zeros', np.round(out['Phi_real_zeros_below_72'], 2), flush=True)
    r = [mP[k] / mS[k] for k in range(KMAX + 1)]
    r = [x / r[0] for x in r]
    out['Phi_over_PhiStar'] = report('Phi / Phi*', r)
    # control: p-ball over Polya bump, p = 4, a = 2 (should pass: 1/Gamma(k/2 + 3 + 1/4) is a multiplier sequence)
    p, a = 4, 2
    rc = [mp.gamma(a + 1) / mp.gamma(mp.mpf(2 * k + 1) / p + a + 1) for k in range(KMAX + 1)]
    rc = [x / rc[0] for x in rc]
    out['control_pball_over_polya_p4_a2'] = report('control p-ball/Polya p=4 a=2', rc)
    # control that must fail: 1/(k+1) (Sum x^k/((k+1)k!) = (e^x - 1)/x has zeros 2 pi i m)
    rf = [mp.mpf(1) / (k + 1) for k in range(KMAX + 1)]
    out['control_1_over_k_plus_1'] = report('control 1/(k+1)', rf)
    json.dump(out, open('multiplier.json', 'w'), indent=1)
