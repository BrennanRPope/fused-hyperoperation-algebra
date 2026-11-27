"""Turán and Jensen-polynomial tests on the three kernel families.
F(x) = sum_k gamma_k x^k / k!  with gamma_k = k! m_2k / (2k)!, m_2k = int t^2k K(t) dt (two-sided).
Real zeros of int K cos(zt) dt  <=>  all Jensen polynomials J^{d,n}(X) = sum_j C(d,j) gamma_{n+j} X^j
are hyperbolic (Pólya–Schur; Pólya 1927 for xi)."""
import mpmath as mp, json, time, sys
mp.mp.dps = 80

def gammas_pball(p, a, K):
    p = mp.mpf(p); a = mp.mpf(a)
    return [mp.factorial(k) * (2 / p) * mp.beta((2 * k + 1) / p, a + 1) / mp.factorial(2 * k) for k in range(K)]

def gammas_polya(p, K):
    p = mp.mpf(p)
    return [mp.factorial(k) * (2 / p) * mp.gamma((2 * k + 1) / p) / mp.factorial(2 * k) for k in range(K)]

def Phi(t):
    s = mp.mpf(0)
    for n in range(1, 30):
        term = (2 * mp.pi ** 2 * n ** 4 * mp.e ** (9 * t) - 3 * mp.pi * n ** 2 * mp.e ** (5 * t)) * mp.e ** (-mp.pi * n ** 2 * mp.e ** (4 * t))
        s += term
        if abs(term) < mp.mpf(10) ** (-mp.mp.dps):
            break
    return s

def gammas_phi(K):
    out = []
    for k in range(K):
        m = 2 * mp.quad(lambda t: t ** (2 * k) * Phi(t), [0, 0.25, 0.5, 1, 1.5, 2.5])
        out.append(mp.factorial(k) * m / mp.factorial(2 * k))
    return out

def hyperbolic(coeffs):
    # coeffs low->high; all positive here.  Count real roots by Sturm-free approach: roots at high precision.
    c = list(reversed(coeffs))
    r = mp.polyroots(c, maxsteps=400, extraprec=400)
    return max(abs(mp.im(z)) / (1 + abs(z)) for z in r) < mp.mpf(10) ** (-30)

def turan(g):
    return [g[k] ** 2 - g[k - 1] * g[k + 1] for k in range(1, len(g) - 1)]

def first_failure(g, D, N):
    for d in range(2, D + 1):
        for n in range(0, N + 1):
            J = [mp.binomial(d, j) * g[n + j] for j in range(d + 1)]
            if not hyperbolic(J):
                return dict(d=d, n=n)
    return None

if __name__ == '__main__':
    D, N = 16, 16
    K = D + N + 2
    import os
    out = json.load(open('jensen.json')) if os.path.exists('jensen.json') else {}
    fams = []
    for p in [1.5, 2.5, 3, 3.5, 5]:
        for a in [0.5 * p, p, p + 0.5, 1.5 * p, 2 * p, 3 * p]:
            if f'pball p={p} a={a}' not in out:
                fams.append((f'pball p={p} a={a}', gammas_pball(p, a, K), (a <= p)))
    for p in [2, 4]:
        for a in [2 * p, 3 * p]:
            fams.append((f'pball p={p} a={a}', gammas_pball(p, a, K), True))
    for p in [1.5, 2.5, 3, 3.5, 4, 5, 6]:
        fams.append((f'polya p={p}', gammas_polya(p, K), p in (2, 4, 6)))
    fams = [f for f in fams if f[0] not in out]
    fams.append(('riemann Phi', gammas_phi(K), None))
    t0 = time.time()
    for lab, g, real_expected in fams:
        tu = turan(g)
        tfail = [k for k, v in enumerate(tu, 1) if v < 0]
        ff = first_failure(g, D, N)
        out[lab] = dict(real_zeros_expected=real_expected, turan_fail_k=tfail[:5], first_jensen_failure=ff)
        print(f"{lab:28s} expected real={real_expected}  Turán fails at k={tfail[:5]}  first Jensen failure={ff}  [{time.time()-t0:.0f}s]", flush=True)
        json.dump(out, open('jensen.json', 'w'), default=str)
