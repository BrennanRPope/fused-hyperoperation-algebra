"""Part C: does a real hat patch see any prime of the Eisenstein integers?
Every hat corner sits on the triangular lattice; we write a lattice point as the
Eisenstein integer z = x + y*w6 (w6 = e^{i pi/3}, norm N = x^2 + x y + y^2).
We test (1) whether hat anchors and hat corners are spread evenly over the residue
classes mod each small Eisenstein prime, and (2) the square-free sieve on hat anchors
against 1/zeta_K(2) and against the k = -1 (union) composite of local removal rates."""
import json, numpy as np, sympy as sp
import hatgen as h

LEVEL = 4
L = h.levels(LEVEL)
patch = h.constructPatch(*L[LEVEL])
# keep only the patch children that make up the four next-level metatiles (H, T, P, F):
# the other seven children are scaffolding, and one of them overlaps a neighbour
KEEP = set(h.NEW_H_CH + h.NEW_P_CH + h.NEW_F_CH + h.NEW_T_CH)
hs = [x for i, (T, g) in enumerate(patch.children) if i in KEEP for x in g.hats(T)]
Ts = np.array([T for T, k in hs]); kinds = [k for T, k in hs]
anc = np.array([h.transPt(T, h.HAT[0]) for T in Ts])
cen = np.array([h.hat_polygon(T).mean(axis=0) for T in Ts])
# origin: the anchor nearest the centre of the patch
c0 = cen.mean(axis=0)
o = anc[np.argmin(((anc - c0)**2).sum(1))]
def lat(P):
    xy = np.array([h.to_hex(p - o) for p in P])
    r = np.round(xy); assert np.abs(xy - r).max() < 1e-6
    return r.astype(np.int64)
A0 = lat(anc)                                  # hat anchors, half-size lattice coordinates
# Every anchor turns out to be a multiple of 2(1+w6): the anchors sit on the centres of the
# kite grid's hexagons. Divide by 2(1+w6) to get hexagon-centre coordinates:
# z/(2(1+w6)) = z(2-w6)/6, and z(2-w6) = (2x+y) + (y-x) w6.
X2, Y2 = 2*A0[:, 0] + A0[:, 1], A0[:, 1] - A0[:, 0]
assert (X2 % 6 == 0).all() and (Y2 % 6 == 0).all(), 'anchors not all on the hexagon-centre lattice'
A = np.stack([X2 // 6, Y2 // 6], 1)             # hat anchors as Eisenstein integers (hex centres)
V = np.unique(lat(np.array([h.transPt(T, p) for T in Ts for p in h.HAT])), axis=0)  # all corners
ang = (np.round(np.degrees(np.arctan2(Ts[:, 3], Ts[:, 0]))).astype(int)) % 360
refl = (Ts[:, 0]*Ts[:, 4] - Ts[:, 1]*Ts[:, 3]) < 0
orient = [f"{a}{'R' if r else ''}" for a, r in zip(ang, refl)]

def eis_primes(maxnorm):
    """One Eisenstein prime per prime ideal of norm <= maxnorm, as (p, r, kind):
    split/ramified: residue map z -> x + y r mod p with r^2 - r + 1 = 0 mod p (w6 -> r);
    inert q: residue (x mod q, y mod q), norm q^2."""
    out = []
    for p in sp.primerange(2, maxnorm + 1):
        if p == 3: out.append((3, 2, 'ramified'))   # w6 = 2 mod (1+w6): 4-2+1=3
        elif p % 3 == 1:
            rs = [r for r in range(p) if (r*r - r + 1) % p == 0]
            out.append((p, rs[0], 'split')); out.append((p, rs[1], 'split'))
        elif p*p <= maxnorm: out.append((p, None, 'inert'))
    return out

def classes(P, pr):
    p, r, kind = pr
    if kind == 'inert': return (P[:, 0] % p) * p + (P[:, 1] % p), p*p
    return (P[:, 0] + P[:, 1]*r) % p, p

def flatness(P, pr):
    c, n = classes(P, pr)
    cnt = np.bincount(c, minlength=n).astype(float)
    exp = len(P) / n
    chi2 = float(((cnt - exp)**2 / exp).sum())
    return dict(norm=n, kind=pr[2], p=int(pr[0]), r=(None if pr[1] is None else int(pr[1])),
                min_share=float(cnt.min()/len(P)*n), max_share=float(cnt.max()/len(P)*n),
                chi2=chi2, dof=n - 1, empty_classes=int((cnt == 0).sum()))

PR = eis_primes(50)
res = {'n_hats': len(hs), 'n_corners': len(V), 'anchors_on_hex_centres_of_index12': True,
       'anchors': [], 'corners': [], 'by_orientation': {}}
for pr in PR:
    res['anchors'].append(flatness(A, pr)); res['corners'].append(flatness(V, pr))
for o_ in sorted(set(orient)):
    m = np.array([q == o_ for q in orient])
    res['by_orientation'][o_] = dict(count=int(m.sum()),
        tests=[flatness(A[m], pr) for pr in PR])

# square-free sieve on hat anchors (exclude the origin)
Az = A[(A != 0).any(1)]
def divisible_by_square(P, pr):
    p, r, kind = pr
    if kind == 'inert':  # q^2 | z  <=>  q^2 | x and q^2 | y
        return (P[:, 0] % (p*p) == 0) & (P[:, 1] % (p*p) == 0)
    if kind == 'ramified':  # (1+w6)^2 = 3 w6, so pi^2 | z <=> 3 | z <=> 3|x, 3|y
        return (P[:, 0] % 3 == 0) & (P[:, 1] % 3 == 0)
    # split: pi^2 | z  <=>  z = 0 in Z[w6]/pi^2 = Z/p^2 with w6 -> r2 (Hensel lift of r)
    q = p*p
    r2 = (r - (r*r - r + 1) * pow(2*r - 1, -1, q)) % q   # Hensel lift of w6 -> r
    assert (r2*r2 - r2 + 1) % q == 0
    return (P[:, 0] + P[:, 1]*r2) % (p*p) == 0

maxN = int((Az[:, 0]**2 + Az[:, 0]*Az[:, 1] + Az[:, 1]**2).max())
PRsq = eis_primes(int(np.sqrt(maxN)) + 2)   # pi^2 | z needs N(pi)^2 <= N(z)
removed = np.zeros(len(Az), bool); local = []
for pr in PRsq:
    d = divisible_by_square(Az, pr)
    removed |= d
    local.append(float(d.mean()))
frac_sqfree = 1 - removed.mean()
# union (k = -1) composite of the measured local rates: a o_{-1} b = a + b - ab
u = 0.0
for a in local: u = u + a - u*a
# 1/zeta_K(2) for K = Q(sqrt(-3)): zeta(2) L(2, chi_-3)
nn = np.arange(1, 20_000_001, dtype=np.float64); chi = np.where(nn % 3 == 1, 1.0, np.where(nn % 3 == 2, -1.0, 0.0))
L2 = float((chi / nn**2).sum())
zK2 = (np.pi**2 / 6) * L2
# the same sieve on every lattice point of a hexagon of equal size, as a control
R = int(np.abs(Az).max())
gx, gy = np.mgrid[-R:R+1, -R:R+1]
G = np.stack([gx.ravel(), gy.ravel()], 1); G = G[(G != 0).any(1)]
G = G[(G[:, 0]**2 + G[:, 0]*G[:, 1] + G[:, 1]**2) <= (Az[:, 0]**2 + Az[:, 0]*Az[:, 1] + Az[:, 1]**2).max()]
remG = np.zeros(len(G), bool)
for pr in PRsq: remG |= divisible_by_square(G, pr)
# prime anchors: anchor is an Eisenstein prime
NMAX = 2 * int(max((A[:,0]**2+A[:,0]*A[:,1]+A[:,1]**2).max(), (G[:,0]**2+G[:,0]*G[:,1]+G[:,1]**2).max())) + 100
isp = np.ones(NMAX + 1, bool); isp[:2] = False
for i in range(2, int(NMAX**0.5) + 1):
    if isp[i]: isp[i*i::i] = False
def eis_prime_mask(P):
    x, y = P[:, 0], P[:, 1]; n = x*x + x*y + y*y
    m1 = isp[n]                                   # norm a rational prime (split or 3)
    q = np.round(np.sqrt(n)).astype(np.int64)
    m2 = (q*q == n) & isp[q] & (q % 3 == 2) & (x % np.maximum(q, 1) == 0) & (y % np.maximum(q, 1) == 0)
    return m1 | m2
primeA = eis_prime_mask(A)
primeG = eis_prime_mask(G)
res['sieve'] = dict(n_anchors=len(Az), frac_squarefree_hat_anchors=float(frac_sqfree),
                    union_composite_of_local_rates=float(u), one_minus_union=float(1 - u),
                    inv_zetaK2=float(1/zK2), zetaK2=float(zK2), L2_chi_minus3=float(L2),
                    frac_squarefree_all_lattice_points=float(1 - remG.mean()), n_lattice_points=len(G),
                    local_rates_first=[(int(pr[0]), pr[2], local[i]) for i, pr in enumerate(PRsq[:12])],
                    frac_prime_hat_anchors=float(primeA.mean()), frac_prime_lattice=float(primeG.mean()))
# share of hexagon centres (inside a disc well within the patch) that carry a hat anchor
nA = A[:, 0]**2 + A[:, 0]*A[:, 1] + A[:, 1]**2
rr = np.sqrt(nA).max() * 0.25   # a disc well inside the patch
inside = nA <= rr**2
Gd = G[(G[:, 0]**2 + G[:, 0]*G[:, 1] + G[:, 1]**2) <= rr**2]
res['anchor_share_of_hex_centres'] = float(inside.sum() / len(Gd[(Gd != 0).any(1)]))   # expect 3/4 (8 kites per hat, 6 per hexagon)
# matched control: hexagon centres in the same region as the hats (within one step of an anchor)
Aset = set(map(tuple, A.tolist()))
nbr = [(0, 0), (1, 0), (-1, 0), (0, 1), (0, -1), (1, -1), (-1, 1)]
region = {(x + dx, y + dy) for x, y in Aset for dx, dy in nbr}
Rg = np.array(sorted(region)); Rg = Rg[(Rg != 0).any(1)]
pr_reg = eis_prime_mask(Rg)
sq_reg = np.zeros(len(Rg), bool)
for pr in PRsq: sq_reg |= divisible_by_square(Rg, pr)
holes = np.array([tuple(p) not in Aset for p in Rg.tolist()])
res['region_control'] = dict(n_hex_centres=len(Rg), share_with_anchor=float(1 - holes.mean()),
    prime_share_all=float(pr_reg.mean()), prime_share_anchors=float(pr_reg[~holes].mean()),
    prime_share_empty_centres=float(pr_reg[holes].mean()),
    squarefree_share_all=float(1 - sq_reg.mean()), squarefree_share_anchors=float(1 - sq_reg[~holes].mean()),
    squarefree_share_empty=float(1 - sq_reg[holes].mean()))
json.dump(res, open('lattice_primes.json', 'w'), indent=1)
np.savez_compressed('patch_level4.npz', T=Ts, kinds=np.array(kinds), anchors=A, corners=V,
                    orient=np.array(orient), primeA=primeA, origin=o)
print('hats', len(hs), 'corners', len(V))
for a, c in zip(res['anchors'], res['corners']):
    print(f"N={a['norm']:4d} {a['kind']:8s} anchors share {a['min_share']:.3f}-{a['max_share']:.3f} chi2={a['chi2']:.1f}/{a['dof']} empty={a['empty_classes']} | corners {c['min_share']:.3f}-{c['max_share']:.3f} chi2={c['chi2']:.1f} empty={c['empty_classes']}")
print(json.dumps(res['sieve'], indent=1))
print('region control', res['region_control'])
print('orientations', {k: v['count'] for k, v in res['by_orientation'].items()})
