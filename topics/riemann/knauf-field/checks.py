import numpy as np, json, mpmath as mp
from knauf import chain, levels
R = []
def ck(name, ok, info=""):
    R.append((name, bool(ok), info)); print(("PASS" if ok else "FAIL"), name, info)
ck("K1 chain(3) = [1 2 3 3 4 5 5 4] (Stern–Brocot denominators)", list(chain(3)) == [1,2,3,3,4,5,5,4])
q, c = levels(20); N = 22
ph = np.arange(N, dtype=float)
for p in range(2, N):
    if ph[p] == p: ph[p::p] *= (1 - 1/p)
m = q < N
ck("K2 level multiplicity = Euler totient φ(q) for q ≤ 21 (n = 20)", np.allclose(c[m], ph[q[m].astype(int)]))
q22, c22 = levels(22); z3 = (c22*q22**-3.0).sum(); t3 = float(mp.zeta(2)/mp.zeta(3))
ck("K3 Z_22(3) ≈ ζ(2)/ζ(3)", abs(z3 - t3) < 0.005, f"{z3:.5f} vs {t3:.5f}")
nn = np.arange(1, 200001.0); p = 2.5
ck("K4 dictionary: Σ n^-p = (p-norm of 1/n)^p", abs((nn**-p).sum() - (np.sum((1/nn)**p)**(1/p))**p) < 1e-12)
h = 0.37 + 0.81j; E = np.log(q22)
lhs = (c22*q22**-(1.5 - h)).sum(); rhs = (c22*q22**-1.5*np.exp(h*E)).sum()
ck("K5 Wick turn of temperature: Z(3/2 − h) = Σ φ^(−3/2) e^(h·E), E = ln φ", abs(lhs - rhs) < 1e-8*abs(lhs))
r = json.load(open("refined.json")); g = r["gam"]
dev = [(row[j][0] - 1.5, row[j][1] - g[j], row[j][2]) for row in r["rows"].values() for j in range(5)]
dev = np.array(dev)
ck("K6 finite chains n = 10..24 have a zero within 0.09 (Re) and 0.06 (Im) of each 3/2 + iγ, γ ≤ 33", np.all(np.abs(dev[:,0]) < 0.09) and np.all(np.abs(dev[:,1]) < 0.06) and np.all(dev[:,2] < 1e-8), f"sd Re = {dev[:,0].std():.4f}")
sd_first = np.std([row[j][0] for n in range(10, 15) for row in [r["rows"][str(n)]] for j in range(5)])
sd_last = np.std([row[j][0] for n in range(20, 25) for row in [r["rows"][str(n)]] for j in range(5)])
ck("K7 the wobble does not shrink: sd(Re) at n = 20..24 is more than half of n = 10..14", sd_last > 0.5*sd_first, f"{sd_first:.4f} -> {sd_last:.4f}")
right = int((dev[:,0] > 0).sum())
ck("K8 zeros fall on both sides of the line (finite chains give no zero-free region)", 20 < right < 55, f"{right} of 75 right of 3/2")
e = json.load(open("errscan.json"))
ck("K9 RMS error at Re β = 1.5 falls at every step n = 10..24", np.all(np.diff(e["1.5"]) < 0), f"{e['1.5'][0]:.4f} -> {e['1.5'][-1]:.4f}")
ck("K10 at Re β = 1.2 the error does not visibly fall (last/first > 0.8)", e["1.2"][-1]/e["1.2"][0] > 0.8, f"ratio {e['1.2'][-1]/e['1.2'][0]:.2f}")
L = json.load(open("l2edge.json"))
primes = [pp for pp in range(2, 100000) if all(pp % d for d in range(2, int(pp**0.5) + 1))]
Ck = float(np.prod([1 - 2/pp**2 + 1/pp**3 for pp in primes]))
ck("K11 ℓ² edge, ζ: Σ_{n≤X} 1/n grows with slope 1 in ln X", abs(L["zeta"]["log_slope"] - 1) < 0.005, f"{L['zeta']['log_slope']:.4f}")
ck("K12 ℓ² edge, Knauf: Σ φ(q)²/q³ grows with slope Π(1 − 2/p² + 1/p³)", abs(L["knauf"]["log_slope"] - Ck) < 0.002, f"{L['knauf']['log_slope']:.4f} vs {Ck:.4f}")
ck("K13 ℓ² edge, Ramanujan τ: Σ τ(n)²/n¹² grows linearly in ln X (slope > 0, steady)", 0.3 < L["tau"]["log_slope"] < 0.5, f"{L['tau']['log_slope']:.4f}")
cc = json.load(open("circle.json"))
ck("K14 magnetic field zeros of the n = 12 chain miss the unit circle at β = 1.5, 3, 5", all(cc["12"][s]["absmin"] > 1.1 for s in ("1.5", "3.0", "5.0")))
f = json.load(open("finite_zeros.json"))
first = [min((z for z in f[str(n)]["zeros"]), key=lambda z: z[1]) for n in (12, 14, 16, 18)]
ck("K15 the lowest finite-size zero creeps toward the transition β = 2 as n grows", all(first[i+1][1] < first[i][1] and first[i+1][0] > first[i][0] for i in range(3)), str([(round(a,2), round(b,2)) for a, b, w in first]))
# control: Davenport-Heilbronn has coefficients of size 1 (same l2 edge 1/2) but a zero off the line
mp.mp.dps = 30
kap = (mp.sqrt(10 - 2*mp.sqrt(5)) - 2)/(mp.sqrt(5) - 1)
chi = [0, 1, 1j, -1j, -1]; chib = [0, 1, -1j, 1j, -1]
def DH(s): return 0.5*(1 - 1j*kap)*mp.dirichlet(s, chi) + 0.5*(1 + 1j*kap)*mp.dirichlet(s, chib)
z = mp.findroot(DH, mp.mpc(0.8085, 85.6993))
ck("K16 control: Davenport–Heilbronn (no Euler product, same ℓ² edge ½) has a zero at Re ≈ 0.8085, off the line", abs(z.real - 0.808517) < 1e-4 and abs(z.imag - 85.699348) < 1e-4, f"{complex(z):.6f}")
print(f"{sum(o for _, o, _ in R)}/{len(R)} pass")
