"""Round 8, item 4: re-centring the blow-up of zeta at s = 1.

(a) Coordinate moves (relabel only): affine w = s - 1 and the Moebius inversion w = 1/(s-1), which sends the
    pole to infinity and the critical line Re s = 1/2 to the circle |w + 1| = 1.  Residue in w: a simple pole of
    residue 1 at s = 1 stays a simple pole (affine) or becomes the point at infinity (inversion); zeros map to zeros.
(b) Fused centre.  a o_k b = a + b + k ab; phi_k(a) = 1 + k a turns o_k into x; the absorbing point is -1/k.
    The fused Euler factor (1 - k u)^(-1/k), u = p^-s, has its local hole where 1 - k p^-s = 0, i.e.
    Re s = log|k| / log p: a different place for every prime.  zeta_k(s) = prod_p (1 - k p^-s)^(-1/k).
    Its global blow-up stays at s = 1 for every k (the u^1 term does not depend on k); the residue changes:
    C_k = prod_p (1 - k/p)^(-1/k) (1 - 1/p)   (|k| < 2).  At k = 2 the local hole of p = 2 lands on s = 1 and the
    blow-up exponent goes from 1 to 3/2.
(c) Standard pole removals: (s-1) zeta(s) (same zeros), xi(s) (same zeros), zeta(s) - 1/(s-1) (Stieltjes; its
    zeros are where zeta = 1/(s-1), a moving a-point, off the line by a little), Hurwitz zeta(s, a) (pole fixed at
    s = 1, residue 1, for every a).
(d) The order dial does move the blow-up: q = 1 (Boltzmann) gives 1/(1 - P(beta)); s = 1 becomes a zero of it
    (P -> infinity logarithmically) and a new simple pole sits at beta = 1.3994333 (Hagedorn).
"""
import json, numpy as np, mpmath as mp
from sympy import primerange
mp.mp.dps = 25
out = {}
PR = list(primerange(2, 200000))

# (a) Moebius image of zeros
zs = [mp.zetazero(n) for n in range(1, 101)]
w = [1/(z-1) for z in zs]
dev = max(abs(abs(x+1)-1) for x in w)
print(f"(a) first 100 zeros under w = 1/(s-1): max | |w+1| - 1 | = {float(dev):.1e}  (on the circle |w+1| = 1)")
out['mobius_w'] = [[float(mp.re(x)), float(mp.im(x))] for x in w]
# residue under inversion: zeta(1 + 1/w) ~ w as w -> 0 along the real axis... check zeta(1+1/w)/ (1/w)^{-1}... i.e. -> w
for W in [0.1, 0.01, 0.001]:
    print(f"    zeta(1 + 1/w) at w = {W}: {mp.nstr(mp.zeta(1+1/mp.mpf(W)), 12)}  (tends to 1 as w -> 0+: the pole now sits at w = infinity)")

# (b) fused residues
def logCk(k):
    return sum(-mp.log(1 - mp.mpf(k)/p)/k + mp.log(1 - mp.mpf(1)/p) for p in PR)
ks = [-1.0, -0.5, 0.5, 1.0, 1.5, 1.9]
Ck = {}
for k in ks:
    tail_bound = sum(abs(k-1)/2/p**2 for p in PR[-1000:])  # crude
    if k == 1.0: Ck[k] = 1.0
    else: Ck[k] = float(mp.e**logCk(k))
    print(f"(b) k = {k:+.1f}: residue C_k = prod_p (1-k/p)^(-1/k)(1-1/p) = {Ck[k]:.6f}")
out['Ck'] = Ck
# check against (s-1) zeta_k(s) at s = 1.0001 for k = 0.5 (Euler product to p < 2e5 is too short; use zeta * correction)
def zeta_k_over_zeta(s, k):
    return mp.e**sum(-mp.log(1 - k*mp.power(p, -s))/k + mp.log(1 - mp.power(p, -s)) for p in PR)
for k in (0.5, -1.0):
    s = mp.mpf('1.000001')
    val = (s-1)*mp.zeta(s)*zeta_k_over_zeta(s, k)
    print(f"    (s-1) zeta_k(s) at s = 1+1e-6, k = {k}: {mp.nstr(val, 8)}")
# k = 2: exponent of the blow-up
k = 2.0
xs = [mp.mpf(10)**(-j) for j in (2, 3, 4, 5)]
vals = []
for x in xs:
    s = 1 + x
    f = mp.zeta(s)*mp.e**sum(-mp.log(1 - k*mp.power(p, -s))/k + mp.log(1 - mp.power(p, -s)) for p in PR)
    vals.append(f)
slopes = [float(mp.log(vals[i+1]/vals[i])/mp.log(xs[i+1]/xs[i])) for i in range(len(xs)-1)]
print(f"(b) k = 2: log-log slope of zeta_2 near s = 1: {[round(x, 4) for x in slopes]}  (pole order 3/2)")
out['k2_slopes'] = slopes
# local holes: Re s = log|k|/log p for k = 3
out['local_holes_k3'] = [[p, float(mp.log(3)/mp.log(p))] for p in PR[:10]]

# (c) Stieltjes removal: zeros of zeta(s) - 1/(s-1) near the first zeros
f = lambda s: mp.zeta(s) - 1/(s-1)
st = []
for n in range(1, 11):
    r = mp.findroot(f, zs[n-1])
    st.append([float(mp.re(r)), float(mp.im(r))])
    print(f"(c) zero of zeta - 1/(s-1) near rho_{n}: {mp.nstr(r, 12)}   (Re - 1/2 = {float(mp.re(r)-0.5):+.2e})")
out['stieltjes_zeros'] = st
for a in (0.3, 0.5, 0.8):
    s = mp.mpf('1.000001')
    print(f"    Hurwitz zeta(s, {a}) * (s-1) at s = 1+1e-6: {mp.nstr((s-1)*mp.zeta(s, a), 10)}")

# (d) order dial q = 1: 1/(1 - P(beta)) on the real axis
def P(b): return mp.primezeta(b)
bs = np.linspace(1.0005, 2.5, 600)
Zb = [float(1/(1-P(mp.mpf(b)))) for b in bs]
wall = mp.findroot(lambda b: P(b) - 1, 1.4)
print(f"(d) q = 1: new pole at beta = {mp.nstr(wall, 10)}; 1/(1-P) at beta = 1.0005: {Zb[0]:.4f} (tends to 0 at beta = 1)")
res = 1/(-mp.diff(P, wall))
print(f"    its residue: {mp.nstr(res, 8)}")
zb = [float(mp.zeta(b)) for b in bs]
out.update({'beta': bs.tolist(), 'Z_boltz': Zb, 'zeta_real': zb, 'wall': float(wall), 'wall_residue': float(res)})
json.dump(out, open('recentre.json', 'w'))
print('saved')
