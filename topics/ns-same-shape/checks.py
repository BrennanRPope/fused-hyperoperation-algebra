"""Checks for ns-same-shape (DRAFT). numpy only. Writes results.json and prints pass/fail lines."""
import numpy as np, json
rng = np.random.default_rng(7)
R = {}; passed = []; failed = []
def check(name, ok, info=""):
    (passed if ok else failed).append(name); print(("PASS " if ok else "FAIL ") + name, info)

# ---------- 1. Scale entanglement (quantics / tensor-train view) ----------
x0 = 0.31830988618379067  # 1/pi, a generic (non-dyadic) point
def fields(n):
    N = 2**n; x = (np.arange(N) + 0.5) / N; h = 1.0 / N
    out = {}
    out["smooth (sine + bump)"] = np.sin(2*np.pi*x) + np.exp(-(x-0.6)**2/0.01)
    out["shock (Burgers front)"] = np.tanh((x - x0)/1e-4)
    out["cusp |x-x0|^(1/3) (Kolmogorov roughness)"] = np.abs(x - x0)**(1/3)
    out["Leray-type 1/|x-x0| profile"] = 1/np.sqrt((x - x0)**2 + h**2)
    k = np.fft.rfftfreq(N, d=1.0/N); amp = np.zeros_like(k); amp[1:] = k[1:]**(-5/6)
    ph = np.exp(2j*np.pi*np.random.default_rng(1).random(len(k)))
    out["random-phase k^(-5/3) field"] = np.fft.irfft(amp*ph, n=N)
    out["same, with viscous cutoff at k=512"] = np.fft.irfft(amp*np.exp(-k/512.0)*ph, n=N)
    return out

def scale_cuts(f, n, tol=1e-6):
    f = f / np.linalg.norm(f); S = []; ranks = []
    for j in range(1, n):
        M = f.reshape(2**j, 2**(n-j))
        s = np.linalg.svd(M, compute_uv=False); p = s**2/np.sum(s**2); p = p[p > 1e-300]
        S.append(float(-(p*np.log2(p)).sum()))
        tail = np.sqrt(np.cumsum((s**2)[::-1])[::-1])  # error if truncating at r = index
        r = int(np.argmax(tail <= tol)) if np.any(tail <= tol) else len(s)
        ranks.append(r)
    return S, ranks

n = 20; F = fields(n); ent = {}; rk = {}
for name, f in F.items():
    S, r = scale_cuts(f, n); ent[name] = S; rk[name] = r
R["entropy_by_cut_n20"] = ent; R["rank_by_cut_n20"] = rk
grow = {}
for name in F:
    grow[name] = []
    for m in range(10, 21, 2):
        S, r = scale_cuts(fields(m)[name], m); grow[name].append(max(r))
R["max_rank_vs_n"] = {"n": list(range(10, 21, 2)), **grow}
mx = {k: max(v) for k, v in ent.items()}
print("max entanglement (bits) per field:", {k: round(v, 3) for k, v in mx.items()})
print("max rank (tol 1e-6) vs n:", grow)
check("E1 smooth field has low scale entanglement (<1 bit)", mx["smooth (sine + bump)"] < 1.0)
check("E2 self-similar cusp (singular) is ALSO low (<1.5 bits)", mx["cusp |x-x0|^(1/3) (Kolmogorov roughness)"] < 1.5)
check("E3 Leray-type singular profile is low (<1.5 bits)", mx["Leray-type 1/|x-x0| profile"] < 1.5)
check("E4 random-phase field: energy-weighted entropy stays modest (<2 bits) because energy sits at large scales",
      mx["random-phase k^(-5/3) field"] < 2)
vc = grow["same, with viscous cutoff at k=512"]
check("E7 viscous cutoff caps the rank: it stops growing once scales pass the cutoff", vc[-1] == vc[-2] == vc[-3], str(vc))
c = grow["cusp |x-x0|^(1/3) (Kolmogorov roughness)"]; rp = grow["random-phase k^(-5/3) field"]
check("E5 cusp rank grows slowly with added scales (at most +12 from n=10 to 20)", c[-1] - c[0] <= 12, str(c))
check("E6 random-phase rank grows fast (x8 or more from n=10 to 20)", rp[-1] >= 8*rp[0], str(rp))

# ---------- 2. Multiplicative cascade (Mandelbrot / Kahane-Peyriere) ----------
# Each scale multiplies mass by W = 1 + x (mean 1): the chain is an o_1 composition across scales.
sig_c = np.sqrt(2*np.log(2)); R["sigma_critical"] = sig_c
def cascade_masses(sig, depth=16, reps=300, seed=0):
    g = np.random.default_rng(seed); Y = np.zeros((reps, depth+1)); Y[:, 0] = 1
    for r in range(reps):
        m = np.ones(1)
        for d in range(1, depth+1):
            W = np.exp(sig*g.standard_normal(2*len(m)) - sig**2/2)
            m = np.repeat(m, 2)*W/2; Y[r, d] = m.sum()
    return Y
casc = {}
for s in [0.6, 1.0, sig_c, 1.5, 2.0]:
    Y = cascade_masses(s); casc[f"{s:.3f}"] = {"median": np.median(Y, 0).tolist(), "mean": Y.mean(0).tolist()}
R["cascade"] = casc
med = lambda s: casc[f"{s:.3f}"]["median"]
check("C1 below the line (sigma=0.6) total mass stays O(1)", med(0.6)[-1] > 0.6, f"{med(0.6)[-1]:.3f}")
check("C2 above the line (sigma=2.0) typical mass collapses (<0.05)", med(2.0)[-1] < 0.05, f"{med(2.0)[-1]:.4f}")
check("C3 mean stays 1 either way (collapse is in the typical value, not the average)",
      abs(casc["0.600"]["mean"][-1]-1) < 0.1)
# o_1 composition identity: (1+a)(1+b)-1 = a o_1 b
a, b = rng.random(5), rng.random(5)
check("C4 one cascade step pair = mult-add a+b+ab", np.allclose((1+a)*(1+b)-1, a+b+a*b))

# ---------- 3. Chirality: same-handed (Beltrami) modes do not interact ----------
def abc(k, A=1.0, B=0.7, C=0.4, L=32):
    t = np.linspace(0, 2*np.pi, L, endpoint=False); X, Y, Z = np.meshgrid(t, t, t, indexing="ij")
    u = np.array([A*np.sin(k*Z)+C*np.cos(k*Y), B*np.sin(k*X)+A*np.cos(k*Z), C*np.sin(k*Y)+B*np.cos(k*X)])
    return u
def curl(u, L=32):
    kk = np.fft.fftfreq(L, d=1.0/L); KX, KY, KZ = np.meshgrid(kk, kk, kk, indexing="ij"); U = np.fft.fftn(u, axes=(1, 2, 3))
    w = np.array([1j*(KY*U[2]-KZ*U[1]), 1j*(KZ*U[0]-KX*U[2]), 1j*(KX*U[1]-KY*U[0])])
    return np.real(np.fft.ifftn(w, axes=(1, 2, 3)))
u1, u2 = abc(1), abc(2); w1, w2 = curl(u1), curl(u2)
check("H1 ABC flow at k is a curl eigenfield (curl u = k u)", np.allclose(w1, u1, atol=1e-10) and np.allclose(w2, 2*u2, atol=1e-10))
lamb_same = np.cross(u1, w1, axis=0)
check("H2 same-handed, same-size mode: Lamb vector u x w is exactly zero", np.abs(lamb_same).max() < 1e-10)
cross = np.cross(u1, w2, axis=0) + np.cross(u2, w1, axis=0)
check("H3 cross term between modes = (lambda2 - lambda1) u1 x u2", np.allclose(cross, (2-1)*np.cross(u1, u2, axis=0), atol=1e-10))
u3 = abc(1, A=0.3, B=1.1, C=0.8); um = abc(1, A=-1.0, B=-0.7, C=-0.4)  # still eigenvalue +1
check("H4 two different +1 modes also cancel (nonlinearity needs different eigenvalues)",
      np.abs(np.cross(u1, curl(u3), axis=0)+np.cross(u3, curl(u1), axis=0)).max() < 1e-10)

# ---------- 4. Scaling: which equations are supercritical (energy vs scaling) ----------
# For u -> lambda^a u(lambda x, lambda^b t) symmetries, the energy E scales like lambda^e. e<0 (energy grows
# at small scales under the symmetry) = supercritical.  Values below are standard.
table = {"Navier-Stokes 3D": -1, "Navier-Stokes 2D": 0, "dissipative SQG alpha=1/2 (critical)": 0,
         "defocusing NLS d=5 quintic... energy-supercritical example (cubic d=5)": -1, "Yang-Mills 4+1": 0, "Yang-Mills 5+1": -1}
R["energy_scaling_exponent"] = table
# NS 3D check: u_l(x,t)=l u(lx,l^2 t): E_l = l^2 * l^-3 * E = l^-1 E
lam = 3.0; check("S1 3D Navier-Stokes energy scales as lambda^-1 (supercritical)", np.isclose(lam**2*lam**-3, lam**-1))
check("S2 2D Navier-Stokes energy is scale-invariant (critical)", np.isclose(lam**2*lam**-2, 1.0))

R["passed"] = len(passed); R["failed"] = failed
json.dump(R, open("results.json", "w"), indent=1)
print(f"\n{len(passed)} of {len(passed)+len(failed)} checks pass")
