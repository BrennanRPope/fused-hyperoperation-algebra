"""DRAFT. How much of a general flow is already 'twisted' (quarter-turn invariant), and which twisted parts
Navier-Stokes keeps closed.  g = quarter turn about y + shift pi/2 along y;  g^4 = 1 on the box.
Sectors: u = u0 + u1 + u2 + u3 with g u_j = i^j u_j  (u0 lives on the quarter-turn space; u0+u2 on the half-turn space).
Run: python3 twist_sectors.py -> twist_sectors.json"""
import json, numpy as np
N = 32
k1 = np.fft.fftfreq(N, 1/N); KX, KY, KZ = np.meshgrid(k1, k1, k1, indexing='ij'); K = np.array([KX, KY, KZ])
K2 = (K**2).sum(0); K2s = np.where(K2 == 0, 1, K2)
fft = lambda f: np.fft.fftn(f, axes=(1, 2, 3)); ifft = lambda F: np.fft.ifftn(F, axes=(1, 2, 3))
def leray(U): return U - K*(K*U).sum(0)/K2s
def g_image(u):
    idx = np.arange(N); Ii, Jj, Ll = np.meshgrid(idx, idx, idx, indexing='ij')
    s = u[:, Ll, (Jj - N//4) % N, (-Ii) % N]
    return np.array([-s[2], s[1], s[0]])
def sectors(u):
    imgs = [u]
    for _ in range(3): imgs.append(g_image(imgs[-1]))
    return [sum((1j**(-j*m))*imgs[m] for m in range(4))/4 for j in range(4)]
def nonlin(u):   # -P[(u.grad)u] for complex fields (bilinear, no conjugation)
    U = fft(u); gu = np.array([[np.fft.ifftn(1j*K[j]*U[i]) for j in range(3)] for i in range(3)])
    return ifft(leray(fft(-np.einsum('jxyz,ijxyz->ixyz', u, gu))))
rng = np.random.default_rng(3)
R = (rng.normal(size=(3, N, N, N)) + 1j*rng.normal(size=(3, N, N, N)))*np.exp(-K2/18)*(K2 > 0)
u = np.real(ifft(leray(R)))
S = sectors(u); E = lambda v: float(np.mean(np.abs(v)**2))
out = {'random field: energy share by sector (1, i, -1, -i)': [E(s)/E(u) for s in S],
       'sectors add back to the field (error)': float(np.abs(sum(S) - u).max())}
# closure: N(u0) stays in sector 0; N(u2) lands in sector 0 (i^2 * i^2 = 1); N(u1) lands in sector 2
for j in range(4):
    n = nonlin(S[j]); Sn = sectors(n); tot = E(n)
    out[f'nonlinear term of sector {j}: share landing in sectors 0..3'] = [E(s)/tot for s in Sn]
json.dump(out, open('twist_sectors.json', 'w'), indent=1)
for k, v in out.items(): print(k, np.round(v, 4) if isinstance(v, list) else v)
