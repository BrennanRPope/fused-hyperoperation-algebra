"""Checks for the even-p proof of the Lee-Yang rule (section N3, even roundness from Polya's end).

Claim: G_a(s) = int_0^1 (1-x^p)^a cos(sx) dx has Taylor coefficients equal to those of
Polya's F_p(s) = int_0^inf exp(-y^p) cos(sy) dy multiplied by q(k) = Gamma(a+1)/Gamma(2k/p + a + 1 + 1/p).
q is 1/Gamma of a linear function of k, an entire function in the Laguerre-Polya class with only
negative zeros, so {q(k)} is a multiplier sequence (Laguerre). Polya 1923: F_p has only real zeros
iff p is an even whole number. So for even p, G_a has only real zeros for every a > -1.

H1  coefficient identity, series vs direct quadrature (several p, a, s, including non-even p)
H2  q(z) zeros: first few lie at -(beta+n)/rho < 0 (sanity), Wright function W = sum z^k q(k)/k!
    has only negative real zeros in |z| < 200 (argument principle vs sign changes), p = 4, 8, 12
H3  even p past 8 (p = 10, 12) and a on both sides of p: every zero of G in the box |Re s| < S,
    |Im s| < 6 is real (argument principle vs real sign changes)
H4  contrast: non-even p above a = p has zeros off the line in the same kind of box (p = 2.5, a = 5 and
    p = 1.5, a = 2.5, where the predicted crossover is inside the box; at p = 9 it is near s = 580)
"""
import json, time, sys
import mpmath as mp
mp.mp.dps = 90

def coeffs(p, a, K):
    p = mp.mpf(p); a = mp.mpf(a)
    return [mp.gamma(a+1)*mp.gamma((2*k+1)/p)/(p*mp.factorial(2*k)*mp.gamma((2*k+1)/p+a+1)) for k in range(K)]

def G_series(c, s):
    w = -s*s; acc = mp.mpf(0)
    for ck in reversed(c):
        acc = acc*w + ck
    return acc

def G_quad(p, a, s):
    p = mp.mpf(p); a = mp.mpf(a)
    return mp.quad(lambda x: (1-x**p)**a*mp.cos(s*x), [0, 0.5, 0.9, 1])

def count_box(c, S, H, N=2400):
    # argument principle on rectangle [-S,S]x[-H,H]
    pts = []
    for j in range(N): pts.append(mp.mpc(-S + 2*S*j/N, -H))
    for j in range(N): pts.append(mp.mpc(S, -H + 2*H*j/N))
    for j in range(N): pts.append(mp.mpc(S - 2*S*j/N, H))
    for j in range(N): pts.append(mp.mpc(-S, H - 2*H*j/N))
    vals = [G_series(c, z) for z in pts]
    tot = mp.mpf(0)
    for i in range(len(vals)):
        tot += mp.arg(vals[(i+1) % len(vals)]/vals[i])
    return float(tot/(2*mp.pi))

def count_real(c, S, N=4000):
    vals = [G_series(c, mp.mpf(S)*(2*j/N - 1) + mp.mpf('1e-9')) for j in range(N+1)]
    return sum(1 for i in range(N) if vals[i]*vals[i+1] < 0)

def W(z, rho, beta):
    s = mp.mpf(0); k = 0; mx = 0
    while True:
        t = z**k/(mp.factorial(k)*mp.gamma(rho*k+beta)); s += t; mx = max(mx, abs(t))
        if k > 30 and abs(t) < mx*mp.mpf(10)**(-mp.mp.dps+10): break
        k += 1
    return s

out = {}; t0 = time.time()
# H1
gaps = []
for (p, a, s) in [(4, 0.3, 5), (4, 7.5, 13), (10, 2.2, 20), (10, 25, 30), (9, 12, 17), (3.3, 1.7, 9), (12, 40, 25)]:
    c = coeffs(p, a, 160)
    g1 = G_series(c, mp.mpf(s)); g2 = G_quad(p, a, mp.mpf(s))
    gaps.append(float(abs(g1-g2)/max(abs(g2), mp.mpf('1e-30'))))
out['H1_relative_gaps'] = gaps
out['H1_pass'] = max(gaps) < 1e-8
# H2  (needs about 150 digits: W is exponentially small on the negative axis)
h2 = []
with mp.workdps(150):
    for p, R in [(4, 200), (8, 200), (12, 100)]:
        rho = mp.mpf(2)/p; beta = mp.mpf('1.3') + mp.mpf(1)/p
        R = mp.mpf(R); N = 3000
        prev = W(R, rho, beta); tot = mp.mpf(0)
        for j in range(1, N+1):
            cur = W(R*mp.expj(2*mp.pi*j/N), rho, beta); tot += mp.arg(cur/prev); prev = cur
        vals = [W(-R*mp.mpf(j)/1500, rho, beta) for j in range(1, 1501)]
        neg = sum(1 for i in range(1499) if vals[i]*vals[i+1] < 0)
        h2.append({'p': p, 'radius': int(R), 'total_in_disc': round(float(tot/(2*mp.pi)), 3), 'negative_real': neg})
        print(h2[-1], flush=True)
out['H2'] = h2
out['H2_pass'] = all(abs(d['total_in_disc'] - d['negative_real']) < 0.01 for d in h2)
# H3, H4
rows = []
for (p, a) in [(10, 3), (10, 10), (10, 25), (12, 6), (12, 30), (2.5, 5), (1.5, 2.5)]:
    c = coeffs(p, a, 220)
    S = mp.mpf(36); H = mp.mpf(6)
    tb = count_box(c, S, H); tr = count_real(c, S)
    rows.append({'p': p, 'a': a, 'box_total': round(tb, 3), 'real': tr})
    print(rows[-1], flush=True)
out['H3_H4'] = rows
out['H3_pass'] = all(abs(r['box_total']-r['real']) < 0.01 for r in rows if float(r['p']) % 2 == 0)
out['H4_pass'] = any(r['box_total']-r['real'] > 0.5 for r in rows if float(r['p']) % 2 != 0 and r['a'] > r['p'])
out['seconds'] = round(time.time()-t0, 1)
json.dump(out, open('even_p_proof_check.json', 'w'), indent=1)
print(json.dumps(out, indent=1))
