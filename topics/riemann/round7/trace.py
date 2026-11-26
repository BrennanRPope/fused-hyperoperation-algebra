"""Round 7, item 4: the explicit formula as a trace formula; dials as orbit weights.

Primes are primitive periodic orbits of period log p; repetition m has period m log p and
amplitude log p / p^{m/2} (a hyperbolic orbit with Lyapunov exponent 1).  The oscillating
part of the zero density is  -(1/pi) sum_{p,m} log p p^{-m/2} cos(E m log p)  (the minus sign
is the "absorption" sign; Gutzwiller's sum has +).

Odds dial t weights repetition m by t^m:  sum_p -log(1 - t p^-s) = sum_N M(t,N) log zeta(Ns),
M(t,N) = (1/N) sum_{d|N} mu(N/d) t^d  (necklace / Witt numbers; cyclotomic identity
1/(1-tu) = prod_N (1-u^N)^{-M(t,N)}).  So the dressed spectrum is the union of rho/N with
multiplicity M(t,N): positive integers for every N exactly when t is a whole number, where
M(t,N) is the dimension of the degree-N piece of the free Lie algebra on t generators.
Homothety k: (1-ku)^{-1/k} has multiplicities M(k,N)/k; whole numbers for all N only at k = +-1
(Bose zeta; Fermi zeta(s)/zeta(2s)).
"""
import json, math
import numpy as np
import mpmath as mp
from sympy import primerange, mobius, divisors
from rmt import zeros_on_line

OUT = {}
def M(t, N): return sum(mobius(N//d)*t**d for d in divisors(N))/N

# (a) zeros from primes
E = np.linspace(5, 60, 2200)
X = 3000
dens = (1/(2*np.pi))*np.log(E/(2*np.pi)) * 0
osc = np.zeros_like(E)
for p in primerange(2, X):
    lp = math.log(p); m = 1
    while p**m < X:
        w = math.exp(-(m*lp/math.log(X))**2*4)      # smooth cutoff
        osc += -(1/np.pi)*lp*p**(-m/2)*np.cos(E*m*lp)*w
        m += 1
smooth = (1/(2*np.pi))*np.log(E/(2*np.pi))
OUT['zeros_from_primes'] = {'E': E.tolist(), 'smooth': smooth.tolist(), 'osc': osc.tolist()}
rz = [float(mp.im(mp.zetazero(n))) for n in range(1, 13)]
OUT['rz'] = rz
peaks = [E[i] for i in range(1, len(E)-1) if (smooth+osc)[i] > (smooth+osc)[i-1] and (smooth+osc)[i] > (smooth+osc)[i+1] and (smooth+osc)[i] > 0.35]
print("(a) density peaks from primes < 3000:", [round(x, 2) for x in peaks[:10]])
print("    Riemann zeros:                     ", [round(x, 2) for x in rz[:10]])
OUT['peaks'] = peaks

# (b) primes from zeros: sum_gamma cos(gamma x) with a smooth window, x in [0.3, 3.2]
g, _ = zeros_on_line('zeta', 1, 2000)
x = np.linspace(0.3, 3.2, 2000)
G = 2000.0
win = np.exp(-(g/G)**2*3)
S = (np.cos(np.outer(x, g)) * win).sum(axis=1)
OUT['primes_from_zeros'] = {'x': x.tolist(), 'S': S.tolist(), 'nzeros': int(len(g))}
dips = [x[i] for i in range(1, len(x)-1) if S[i] < S[i-1] and S[i] < S[i+1] and S[i] < -20]
print(f"(b) {len(g)} zeros: dips of sum cos(gamma x) at x =", [round(v, 3) for v in dips])
print("    log of prime powers:", [round(math.log(n), 3) for n in (2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 17, 19, 23)])
OUT['dips'] = dips

# (c) cyclotomic identity and multiplicities
mp.mp.dps = 25
rows = []
for t in [0.5, 2, 3, -1]:
    s = mp.mpf(3)
    euler = -sum(mp.log(1 - t*mp.mpf(p)**(-s)) for p in primerange(2, 200000))
    neck = sum(M(t, N)*mp.log(mp.zeta(N*s)) for N in range(1, 60))
    rows.append([t, float(euler), float(neck)])
    print(f"(c) t={t}: sum_p -log(1-t p^-3) = {float(euler):.12f}; sum_N M(t,N) log zeta(3N) = {float(neck):.12f}")
OUT['cyclotomic'] = rows
mult = {str(t): [float(M(t, N)) for N in range(1, 9)] for t in [0.5, 1, 2, 3, -1]}
multk = {str(k): [float(M(k, N)/k) for N in range(1, 9)] for k in [1, -1, 2, 0.5]}
for t, v in mult.items(): print(f"    odds t={t}: multiplicities of rho/N, N=1..8:", [round(a, 4) for a in v])
for k, v in multk.items(): print(f"    homothety k={k}: multiplicities M(k,N)/k:", [round(a, 4) for a in v])
OUT['mult_t'] = mult; OUT['mult_k'] = multk
# integrality: M(t,N) whole for all N <= 30 iff t whole (scan t on a grid)
grid = [i/4 for i in range(-8, 17)]
integral = [t for t in grid if all(abs(M(t, N)-round(M(t, N))) < 1e-9 for N in range(1, 31))]
print("    t in [-2,4] (step 1/4) with whole multiplicities for N<=30:", integral)
OUT['integral_t'] = integral
json.dump(OUT, open('trace.json', 'w'))
