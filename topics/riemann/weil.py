"""Pass 1, item 3: local factors as points of the dial plane.
A degree-2 local factor 1 + k u + d u^2 (u = p^-s) is the conic dial point (k, d). Its two roots have the same size
|u| = d^(-1/2) (both on the 'critical circle', i.e. Re s = 1/2 when d = p) exactly when Delta = k^2 - 4d <= 0:
the elliptic side, positive curvature -Delta/4.  Over a finite field the whole zeta function of a curve is a finite
product of such factors (genus g gives g factors), so there the Riemann hypothesis reads: every dial point is elliptic.
We count points on a genus-1 and a genus-2 curve and place every factor.  Writes weil.json."""
import json, numpy as np
from sympy import primerange, legendre_symbol

def nonres(p):
    for r in range(2, p):
        if legendre_symbol(r, p) == -1: return r

def count_Fp(coeffs, p):
    # points on y^2 = f(x) over F_p, plus points at infinity (1 for odd degree)
    x = np.arange(p, dtype=np.int64)
    f = np.zeros(p, dtype=np.int64)
    for c in coeffs: f = (f*x + c) % p
    chi = np.array([0 if v == 0 else legendre_symbol(int(v), p) for v in f])
    return p + int(chi.sum()) + 1

def count_Fp2(coeffs, p):
    # F_{p^2} = F_p[w]/(w^2 - r); chi(z) = legendre(N(z)), N(a + b w) = a^2 - r b^2
    r = nonres(p)
    a, b = np.meshgrid(np.arange(p, dtype=np.int64), np.arange(p, dtype=np.int64))
    a = a.ravel(); b = b.ravel()
    fa = np.zeros_like(a); fb = np.zeros_like(a)
    for c in coeffs:  # Horner: f <- f*x + c with x = a + b w
        na = (fa*a + r*fb % p * b + c) % p
        nb = (fa*b + fb*a) % p
        fa, fb = na, nb
    N = (fa*fa - r*(fb*fb % p)) % p
    leg = np.zeros(p, dtype=np.int64)
    for v in range(1, p): leg[v] = legendre_symbol(v, p)
    chi = leg[N]
    return p*p + int(chi.sum()) + 1

out = {}
# genus 1: E: y^2 = x^3 + x + 1  (discriminant -16*31: bad at 2, 31)
E = [1, 0, 1, 1]
g1 = []
for p in primerange(3, 3000):
    if p == 31: continue
    N1 = count_Fp(E, p); a = p + 1 - N1
    g1.append(dict(p=p, a=a, k=-a, d=p, Delta=a*a - 4*p, theta=float(np.arccos(a/(2*np.sqrt(p))))))
out['genus1'] = g1
# genus 2: C: y^2 = x^5 - x + 1  (check good reduction: f squarefree mod p)
from sympy import Poly, symbols, discriminant
X = symbols('x')
disc = int(discriminant(X**5 - X + 1, X))
C = [1, 0, 0, 0, -1, 1]
g2 = []
for p in primerange(3, 400):
    if disc % p == 0: continue
    N1 = count_Fp(C, p); N2 = count_Fp2(C, p)
    s1 = p + 1 - N1; s2 = p*p + 1 - N2          # sums of the four Frobenius roots and of their squares
    c2 = (s1*s1 - s2)//2
    # P(u) = (1 - A u + p u^2)(1 - B u + p u^2): A + B = s1, A B = c2 - 2p
    D = s1*s1 - 4*(c2 - 2*p)
    A = (s1 + np.sqrt(D))/2 if D >= 0 else complex(s1, np.sqrt(-D))/2
    B = (s1 - np.sqrt(D))/2 if D >= 0 else complex(s1, -np.sqrt(-D))/2
    g2.append(dict(p=p, N1=N1, N2=N2, s1=s1, c2=c2, A=float(np.real(A)), B=float(np.real(B)), AB_real=bool(D >= 0),
                   DeltaA=float(np.real(A)**2 - 4*p), DeltaB=float(np.real(B)**2 - 4*p)))
out['genus2'] = g2
out['genus2_disc'] = disc
json.dump(out, open('weil.json', 'w'), indent=1)
print('genus1 primes', len(g1), 'all elliptic:', all(r['Delta'] <= 0 for r in g1))
print('genus2 primes', len(g2), 'real pairs:', all(r['AB_real'] for r in g2), 'all elliptic:', all(r['DeltaA'] <= 1e-9 and r['DeltaB'] <= 1e-9 for r in g2))
print(g2[:3])
