# DRAFT. Bridges study, candidates C1-C3 (Nov 28, 2026). Small exact or numeric checks.
import json, numpy as np, sympy as sp
rng = np.random.default_rng(1128)
res, out = {}, []
def log(name, ok, note): res[name] = bool(ok); out.append(f"{name} {'PASS' if ok else 'FAIL'}  {note}")
def ts(p, q):
    p = np.asarray(p, float).ravel(); p = p[p > 1e-300]
    return float(-(p*np.log(p)).sum()) if abs(q-1) < 1e-12 else float((1-(p**q).sum())/(q-1))

# ---------- C1: non-additive horizon entropies on the law's dial ----------
# D1 Power-of-area entropies (Tsallis-Cirto gamma*A^delta; Barrow delta = 1 + Delta/2) are not a fixed-k law:
#    split one horizon into patches A1, A2; k_eff = (S12 - S1 - S2)/(S1 S2) changes with the split and the size.
def keff_pow(A1, A2, d, g=1.0): S = lambda A: g*A**d; return (S(A1+A2)-S(A1)-S(A2))/(S(A1)*S(A2))
vals = [keff_pow(A1, A2, 1.5) for A1 in (1, 10, 100) for A2 in (1, 10, 100)]
ok1 = (max(vals)/min(vals) > 10) and all(abs(keff_pow(a, b, 1.0)) < 1e-12 for a in (1, 10) for b in (3, 7))
log("D1", ok1, f"A^1.5 patches: k_eff spans {min(vals):.3g} to {max(vals):.3g} (no single k); delta=1 gives k_eff=0")
res["D1_keff_range"] = [min(vals), max(vals)]
# D2 Barrow: G_eff = A^(1-delta)/(4 delta gamma) with delta = 1 + Delta/2 varies with area unless Delta = 0.
A, gam, Dl = sp.symbols('A gamma Delta', positive=True)
Geff = 1/(4*sp.diff(gam*A**(1+Dl/2), A))
ok2 = sp.simplify(sp.diff(Geff.subs(Dl, 0), A)) == 0 and sp.simplify(sp.diff(Geff, A)) != 0
log("D2", ok2, f"Barrow G_eff = {sp.simplify(Geff)}: constant only at Delta = 0")
# D3 The Renyi-of-black-hole construction S_R = ln(1 + lam S_BH)/lam is exactly the law's unlift:
#    if S_BH composes by the law at k = lam, S_R adds.
lam, s1, s2 = sp.symbols('lam s1 s2', positive=True)
SR = lambda s: sp.log(1+lam*s)/lam
ok3 = sp.simplify(sp.expand_log(SR(s1+s2+lam*s1*s2) - SR(s1) - SR(s2), force=True)) == 0
log("D3", ok3, "ln(1+lam*S)/lam turns a+b+lam*ab into plain addition (the law's lift and unlift)")
# D4 Linear mixing: if S_BH were a state entropy composing at k = lam, an unrelated coin would multiply
#    what is left to learn by 1 + lam*S(coin) (round 8 of spin-i-tie, re-run here at three lams).
def left(P, q): return ts(P, q) - ts(P.sum(1), q)
def gcoin(P, c): return np.einsum('ab,c->acb', P, c).reshape(-1, P.shape[1])
r4 = []
for lmb in (-0.3, 0.2, 0.5):
    q = 1-lmb
    for _ in range(20):
        P = rng.random((3, 4)); P /= P.sum(); c = np.array([.3, .7])
        r4.append(abs(left(gcoin(P, c), q)/left(P, q) - (1+lmb*ts(c, q))))
log("D4", max(r4) < 1e-10, f"coin multiplies the leftover by 1 + lam*S(coin), so linear mixing needs lam = 0 (max err {max(r4):.1e})")
# D5 Dong's dictionary: a law-dial lam on area means Renyi order n = 1 - lam and brane tension
#    T = (n-1)/(4nG) = -lam/(4G(1-lam)): lam > 0 needs NEGATIVE tension, lam < 0 positive tension.
G = sp.symbols('G', positive=True); n = 1-lam
T = sp.simplify((n-1)/(4*n*G))
ok5 = all(sp.sign(T.subs({lam: v, G: 1})) == s for v, s in [(sp.Rational(1, 5), -1), (sp.Rational(-1, 5), 1)]) and T.subs(lam, 0) == 0
log("D5", ok5, f"tension T = {T}: lam>0 -> T<0 (excess angle), lam<0 -> T>0 (deficit), lam=0 -> no brane")
# D6 The usual motive (S ~ M^2 is non-additive when black holes merge) is not a fixed-k law either:
#    S(M1+M2) = S1 + S2 + 2 sqrt(S1 S2), so k_eff = 2/sqrt(S1 S2) depends on the sizes.
ke6 = [2/np.sqrt(4*np.pi*m1**2*4*np.pi*m2**2) for m1, m2 in [(1, 1), (1, 10), (10, 10)]]
ok6 = all(abs((4*np.pi*(m1+m2)**2 - 4*np.pi*m1**2 - 4*np.pi*m2**2)/(16*np.pi**2*m1**2*m2**2) - 2/np.sqrt(16*np.pi**2*m1**2*m2**2)) < 1e-12
          for m1, m2 in [(1, 1), (1, 10), (10, 10)]) and ke6[0]/ke6[2] > 50
log("D6", ok6, "merger 'non-additivity' gives k_eff = 2/sqrt(S1 S2), different for every pair of sizes")


# D7 The power-of-area entropies sit on the existing roundness dial p (the l^p sum: cross-polytope at p=1, ball at 2):
#    gamma*A^delta composes patches as S12 = (S1^p + S2^p)^(1/p) with p = 1/delta. Area law: p = 1.
#    Tsallis-Cirto delta = 3/2 and Barrow at its largest Delta = 1 both give p = 2/3 (the non-convex side).
e7 = 0.0
for d in (1.0, 1.25, 1.5, 2.0):
    p = 1/d
    for _ in range(50):
        A1, A2, g = rng.uniform(0.1, 50), rng.uniform(0.1, 50), rng.uniform(0.5, 2)
        S1, S2 = g*A1**d, g*A2**d
        e7 = max(e7, abs((S1**p+S2**p)**(1/p) - g*(A1+A2)**d)/(g*(A1+A2)**d))
log("D7", e7 < 1e-12 and abs(1/1.5 - 2/(2+1)) < 1e-15, f"A^delta = the l^p sum with p = 1/delta (max rel err {e7:.1e}); Tsallis-Cirto and Barrow-max both p = 2/3")

# ---------- C2: Euclidean Kerr from the river's shift by i ----------
# E1 Lorentzian: the Kerr-Schild scalar M r^3/(r^4 + a^2 z^2) is Re[M/R] for the source shifted to z = i a.
def r_obl(rho, z, a):
    b = rho**2+z**2-a**2; return np.sqrt((b+np.sqrt(b*b+4*a*a*z*z))/2)
e1 = 0.0
for _ in range(200):
    rho, z, a = rng.uniform(0.2, 3), rng.uniform(-3, 3), rng.uniform(0.1, 1.5)
    r = r_obl(rho, z, a); R = np.sqrt(complex(rho**2+(z-1j*a)**2))
    e1 = max(e1, abs(abs((1/R).real) - r**3/(r**4+a*a*z*z)))
log("E1", e1 < 1e-12, f"Kerr's scalar = Re of a Coulomb field at z = i a (max err {e1:.1e})")
# E2 Euclidean (a = i alpha): the same scalar becomes an ordinary real Coulomb field of a point moved to
#    z = -alpha on the axis: the two i's multiply to -1 and the source is real.
def r_pro(rho, z, al):
    b = rho**2+z**2+al**2; return np.sqrt((b+np.sqrt(b*b-4*al*al*z*z))/2)
e2 = 0.0
for _ in range(200):
    rho, z, al = rng.uniform(0.2, 3), rng.uniform(-3, 3), rng.uniform(0.1, 1.5)
    r = r_pro(rho, z, al); Rr = np.sqrt(rho**2+(z+al)**2); Rl = np.sqrt(rho**2+(z-al)**2)
    target = r**3/(r**4-al*al*z*z)
    e2 = max(e2, abs(0.5*(1/Rr+1/Rl) - target))
ok2e = e2 < 1e-12
log("E2", ok2e, f"Euclidean scalar = average of real Coulomb fields at z = +alpha and -alpha (max err {e2:.1e})")
# E3 Smoothness of the Euclidean horizon fixes the period: kappa from ds^2 ~ N^2 dtau^2 + dr^2/F,
#    computed on the Euclidean section and independent of the polar angle; continued back it is Kerr's.
r, M, al, th, a = sp.symbols('r M alpha theta a', positive=True)
DeltaE = r**2 - 2*M*r - al**2; rho2 = r**2 - al**2*sp.cos(th)**2
Xi = (r**2-al**2)**2 + DeltaE*al**2*sp.sin(th)**2
N2 = DeltaE*rho2/Xi; F = DeltaE/rho2
rp = M + sp.sqrt(M**2+al**2)
kap2 = sp.simplify((sp.diff(N2, r)*sp.diff(F, r)/4).subs(r, rp))
kapE = (rp-(M-sp.sqrt(M**2+al**2)))/(2*(rp**2-al**2))
ok3e = sp.simplify(kap2 - kapE**2) == 0
rpL = M+sp.sqrt(M**2-a**2); rmL = M-sp.sqrt(M**2-a**2)
kapL = (rpL-rmL)/(2*(rpL**2+a**2))
ok3e = ok3e and sp.simplify(kapE.subs(al, sp.I*a) - kapL) == 0
log("E3", ok3e, "Euclidean smoothness gives kappa_E with no theta dependence; alpha -> -i a gives Kerr's kappa")
# E4 Entropy from the Euclidean action I = beta*M/2 (Gibbons-Hawking 1977): S = beta(M - Omega J) - I = A/4.
T = kapL/(2*sp.pi); beta = 1/T; Om = a/(rpL**2+a**2); J = M*a
S = sp.simplify(beta*(M-Om*J) - beta*M/2)
ok4e = sp.simplify(S - sp.pi*(rpL**2+a**2)) == 0
ok4e = ok4e and sp.simplify(2*T*sp.pi*(rpL**2+a**2) + 2*Om*J - M) == 0
log("E4", ok4e, "S = A/4 = pi(r+^2 + a^2) from I = beta M/2, and Smarr M = 2TS + 2 Omega J")
res["E_numbers_M1_a05"] = {"T": float(T.subs({M: 1, a: 0.5})), "S": float(sp.pi*(rpL**2+a**2).subs({M: 1, a: 0.5}))}

# ---------- C3: the law's extra term across the Page time ----------
# F1 Average purity of a random pure state on X (dim dX) and the rest (dim dY): (dX+dY)/(dX dY + 1)  (Lubkin 1978).
def haar(d): v = rng.normal(size=d)+1j*rng.normal(size=d); return v/np.linalg.norm(v)
dX, dY = 4, 8; ps = []
for _ in range(4000):
    psi = haar(dX*dY).reshape(dX, dY); rho = psi@psi.conj().T; ps.append(np.trace(rho@rho).real)
lub = (dX+dY)/(dX*dY+1)
log("F1", abs(np.mean(ps)-lub) < 4*np.std(ps)/np.sqrt(len(ps)), f"Haar mean purity {np.mean(ps):.4f} vs Lubkin {lub:.4f}")
# F2 Two saddles: purity ~ 1/dX (Hawking, disconnected) + 1/dY (replica wormhole, connected), as in PSSY 2019.
dX, dY = 2.0**12, 2.0**9
log("F2", abs(((dX+dY)/(dX*dY+1)) - (1/dX+1/dY)) < 1e-6*(1/dX+1/dY)*10, "Lubkin = Hawking term 1/d_X + wormhole term 1/d_rest, to 1e-6")
# F3 Evaporation: N qubits; after r emitted, early radiation R1 (first half) and late R2 (second half), black hole N - r.
#    q = 2 Tsallis: S = 1 - purity. Cross term c = k_eff - (1 - q). Hawking saddle alone: c = 0 at every step.
N = 24
def purity(nX, mode):
    dX, dY = 2.0**nX, 2.0**(N-nX)
    return 1/dX if mode == "hawking" else (dX+dY)/(dX*dY+1)
def cross(rq, mode):
    n1 = rq//2; n2 = rq-n1
    s1, s2, s12 = (1-purity(n1, mode)), (1-purity(n2, mode)), (1-purity(rq, mode))
    return (s12-s1-s2)/(s1*s2) - (1-2)
steps = list(range(4, N+1, 2))
cH = [cross(x, "hawking") for x in steps]; cW = [cross(x, "full") for x in steps]
early = [c for x, c in zip(steps, cW) if x <= N//2]; late = [c for x, c in zip(steps, cW) if x >= N-2]
ok3c = max(abs(c) for c in cH) < 1e-12 and max(abs(c) for c in early) < 1e-3 and min(abs(c) for c in late) > 1e-3
log("F3", ok3c, f"cross term: Hawking-only exactly 0; with wormholes |c| <= {max(abs(c) for c in early):.1e} before the Page time, {cW[-1]:.3g} at the end")
res["F3_track"] = {"steps": steps, "hawking": cH, "full": cW, "N": N}
n_ok = sum(1 for v in res.values() if v is True); n_all = sum(1 for v in res.values() if isinstance(v, bool))
out.append(f"{n_ok}/{n_all} passed")
print("\n".join(out)); json.dump(res, open("results.json", "w"), indent=1)
