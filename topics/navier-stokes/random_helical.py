"""Random velocity fields with a chosen handedness mix. Measures enstrophy
production P = <w.(w.grad)u> split into handedness corners, normalised by
enstrophy^(3/2), averaged over seeds. h = share of energy in right-handed waves.
Run: python3 random_helical.py"""
import numpy as np, itertools
N = 32
k1 = np.fft.fftfreq(N, 1/N)
KX, KY, KZ = np.meshgrid(k1, k1, k1, indexing='ij')
K = np.array([KX, KY, KZ]); K2 = (K**2).sum(0); K2[0,0,0] = 1; KA = np.sqrt(K2)
ifs = lambda f: np.real(np.fft.ifftn(f))
def curlh(fh): return 1j*np.cross(K, fh, axis=0)
def basis():
    # helical basis h+(k): i k x h+ = |k| h+
    e = np.cross(K, np.array([np.ones_like(KX)*0.3, np.ones_like(KX)*0.5, np.ones_like(KX)*0.7]), axis=0)
    e /= np.sqrt((e**2).sum(0))+1e-30
    f = np.cross(K, e, axis=0)/KA
    hp = (e - 1j*f)/np.sqrt(2)
    return hp
hp = basis()
def field(h, rng):
    spec = KA**4*np.exp(-2*(KA/3)**2)  # energy peaked near k~3
    amp = np.sqrt(spec)/KA
    ap = amp*np.sqrt(h)*np.exp(2j*np.pi*rng.random(KA.shape))
    am = amp*np.sqrt(1-h)*np.exp(2j*np.pi*rng.random(KA.shape))
    hm = np.conj(hp)  # left-handed partner for this k
    uh = ap*hp + am*hm
    uh[:,0,0,0] = 0
    u = np.real(np.fft.ifftn(uh, axes=(1,2,3)))  # take real part -> real field
    return np.fft.fftn(u, axes=(1,2,3))
def split(uh):
    c = curlh(uh)/KA
    return {'+': 0.5*(uh+c), '-': 0.5*(uh-c)}
def prod(uh):
    us = split(uh)
    w = {s: np.real(np.fft.ifftn(curlh(us[s]), axes=(1,2,3))) for s in us}
    G = {s: np.array([[ifs(1j*K[j]*us[s][i]) for j in range(3)] for i in range(3)]) for s in us}
    out = {}
    for a,b,c in itertools.product('+-', repeat=3):
        out[a+b+c] = np.mean((w[a]*np.einsum('jxyz,ijxyz->ixyz', w[b], G[c])).sum(0))
    wt = w['+']+w['-']; ens = np.mean((wt**2).sum(0))
    Hp = np.mean((np.real(np.fft.ifftn(us['+'],axes=(1,2,3)))*w['+']).sum(0))
    Hm = np.mean((np.real(np.fft.ifftn(us['-'],axes=(1,2,3)))*w['-']).sum(0))
    return out, ens, Hp, Hm
for h in [0.5, 0.75, 0.9, 1.0]:
    tot=[]; single=[]; mixed=[]; hel=[]
    for seed in range(12):
        out, ens, Hp, Hm = prod(field(h, np.random.default_rng(seed)))
        s = out['+++']+out['---']; T = sum(out.values())
        tot.append(T/ens**1.5); single.append(s/ens**1.5); mixed.append((T-s)/ens**1.5); hel.append(Hp/(Hp-Hm))
    print(f"h={h:4.2f} right-handed share of helicity={np.mean(hel):.2f}  P/Omega^1.5 total={np.mean(tot):+.4f}±{np.std(tot)/np.sqrt(12):.4f}  single-handed={np.mean(single):+.4f}  mixed={np.mean(mixed):+.4f}")
