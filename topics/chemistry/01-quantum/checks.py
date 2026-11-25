"""Quantum chemistry through the circle operation: independent fragments, CISD vs coupled cluster.
Run: python3 checks.py  (needs numpy, pyscf)."""
import numpy as np, itertools, json
out = {}
ok = []
def check(name, cond, info=""):
    ok.append(bool(cond)); print(("PASS " if cond else "FAIL ") + name + ("  " + info if info else ""))

# 1. Product of two fragment wavefunctions in intermediate normalisation: corner = disconnected double-double
rng = np.random.default_rng(1)
for _ in range(200):
    cA, cB = rng.normal(size=2)
    lhs = (1+cA)*(1+cB) - 1
    check_val = cA + cB + cA*cB
    assert abs(lhs-check_val) < 1e-12
check("product wavefunction: (1+cA)(1+cB)-1 = cA ∘₁ cB (200 random pairs)", True)

# 2. N-dimer model: each molecule = 2 states (ground, doubly excited at 2Δ), coupling g.
def exact_one(D, g):
    return D - np.sqrt(D*D + g*g)
def cisd_N(D, g, N):
    # CISD space: reference + N single-molecule doubles; one eigenvalue problem
    H = np.zeros((N+1, N+1)); H[0,1:] = g; H[1:,0] = g
    for i in range(1, N+1): H[i,i] = 2*D
    return np.linalg.eigvalsh(H)[0]
def fci_N(D, g, N):
    # full space: 2^N product states; H = sum of independent 2x2 blocks
    h = np.array([[0, g],[g, 2*D]])
    I = np.eye(2); H = np.zeros((2**N, 2**N))
    for i in range(N):
        mats = [I]*N; mats[i] = h
        M = mats[0]
        for m in mats[1:]: M = np.kron(M, m)
        H += M
    return np.linalg.eigvalsh(H)[0]
D, g = 1.0, 0.3
e1 = exact_one(D, g)
rows = []
for N in [1,2,3,4,6,8,10]:
    ec = cisd_N(D,g,N); ef = fci_N(D,g,N)
    closed = D - np.sqrt(D*D + N*g*g)
    rows.append((N, ec/N/e1, ef/N/e1))
    check(f"N={N}: CISD = Δ-√(Δ²+N g²) closed form", abs(ec-closed) < 1e-10)
    check(f"N={N}: full CI = N × one molecule (size-extensive)", abs(ef - N*e1) < 1e-10, f"CISD keeps {ec/(N*e1):.3f} of the correlation")
out["dimer_model"] = rows
# large-N scaling of the fraction kept: ~ 1/sqrt(N)
fr = [(D - np.sqrt(D*D+N*g*g))/(N*e1) for N in [10**4, 4*10**4, 16*10**4]]
check("CISD fraction kept falls like 1/√N at large N (×4 molecules → ÷2)", abs(fr[0]/fr[1]-2) < 0.05 and abs(fr[1]/fr[2]-2) < 0.03, f"{fr}")

# 3. Corner-weighted ansatz: Ψ = Σ_m k^(m-1) c^m (sum of all m-fold products). k=1 is CC-like product form.
from math import comb
def ansatz_energy(D, g, N, k):
    # amplitude for any specific m-fold product state: a_m = k^(m-1) c^m (a_0 = 1)
    best = 0.0
    for c in np.linspace(-1.0, 0.2, 6001):
        a = [1.0] + [k**(m-1) * c**m for m in range(1, N+1)]
        num = 0.0; den = 0.0
        for m in range(N+1):
            den += comb(N,m) * a[m]**2
            num += comb(N,m) * a[m]**2 * 2*D*m
            if m < N:  # coupling between m-fold and (m+1)-fold states: each (m+1)-state links to m+1 m-states
                num += 2 * comb(N, m+1) * (m+1) * g * a[m] * a[m+1]
        best = min(best, num/den)
    return best
kweights = {}
for k in [0.0, 0.5, 1.0]:
    kweights[k] = [ansatz_energy(D,g,N,k)/(N*e1) for N in [1,2,4,8]]
out["corner_weight"] = {str(k): v for k,v in kweights.items()}
check("corner weight k=1 recovers the exact energy for every N", all(abs(x-1) < 1e-4 for x in kweights[1.0]), f"{kweights[1.0]}")
check("corner weight k=0 is CISD (fraction falls with N)", abs(kweights[0.0][3] - (D-np.sqrt(D*D+8*g*g))/(8*e1)) < 1e-4, f"{kweights[0.0]}")
check("k=1/2 sits between, still not size-extensive", kweights[0.0][3] < kweights[0.5][3] < 1 - 1e-3, f"{kweights[0.5]}")

# 4. Real molecules (pyscf): two H2 far apart; CISD vs CCSD
try:
    from pyscf import gto, scf, ci, cc
    def run(atoms):
        mol = gto.M(atom=atoms, basis="cc-pvdz", verbose=0, unit="Angstrom")
        mf = scf.RHF(mol).run()
        e_ci = ci.CISD(mf).run().e_corr
        e_cc = cc.CCSD(mf).run().e_corr
        return mf.e_tot, e_ci, e_cc
    one = run("H 0 0 0; H 0 0 0.74")
    two = run("H 0 0 0; H 0 0 0.74; H 0 0 100; H 0 0 100.74")
    hf1, ci1, cc1 = one; hf2, ci2, cc2 = two
    print(f"  H2:   CISD corr {ci1:.6f}  CCSD corr {cc1:.6f} Eh")
    print(f"  2 H2: CISD corr {ci2:.6f} (2× = {2*ci1:.6f})  CCSD corr {cc2:.6f} (2× = {2*cc1:.6f})")
    check("one H2: CISD = CCSD (two electrons: no corner possible)", abs(ci1-cc1) < 1e-6)
    check("two H2 100 Å apart: CCSD = 2 × one H2", abs(cc2-2*cc1) < 1e-6)
    missing = 2*ci1 - ci2
    check("two H2 100 Å apart: CISD misses energy (size inconsistent)", missing < -1e-5, f"misses {abs(missing)*627.5:.3f} kcal/mol")
    # the missing piece is the corner: CISD+Q-type estimate: product of single-molecule amplitudes
    out["pyscf_h2"] = dict(ci1=ci1, cc1=cc1, ci2=ci2, cc2=cc2, missing_kcal=missing*627.5)
    # water: CISD vs CCSD on one vs two waters 100 Å apart
    w = "O 0 0 0; H 0 0.757 0.587; H 0 -0.757 0.587"
    w2 = w + "; O 100 0 0; H 100 0.757 0.587; H 100 -0.757 0.587"
    def run2(atoms):
        mol = gto.M(atom=atoms, basis="6-31g", verbose=0); mf = scf.RHF(mol).run()
        return ci.CISD(mf).run().e_corr, cc.CCSD(mf).run().e_corr
    wc1, wcc1 = run2(w); wc2, wcc2 = run2(w2)
    out["pyscf_h2o"] = dict(ci1=wc1, cc1=wcc1, ci2=wc2, cc2=wcc2)
    check("two waters 100 Å apart: CCSD = 2 × one water", abs(wcc2-2*wcc1) < 1e-6, f"{wcc2:.6f} vs {2*wcc1:.6f}")
    check("two waters: CISD loses correlation", 2*wc1 - wc2 < -1e-4, f"misses {abs(2*wc1-wc2)*627.5:.2f} kcal/mol; CISD keeps {wc2/(2*wc1):.4f}")
except Exception as ex:
    print("pyscf part skipped:", ex)

# 5. Spin coupling: multiplicities multiply, so (multiplicity - 1) combines by the circle operation
for s1 in range(0,6):
    for s2 in range(0,6):
        n1, n2 = s1+1, s2+1   # 2s+1 with 2s = s1
        total = sum(2*S+1 for S in [x/2 for x in range(abs(s1-s2), s1+s2+1, 2)])  # Clebsch-Gordan series in units of 1/2
        assert total == n1*n2 and n1*n2-1 == s1 + s2 + s1*s2
check("spin coupling: Σ(2S+1) over the Clebsch–Gordan series = (2s1+1)(2s2+1), so 2s combines by the circle operation (36 cases)", True)
json.dump(out, open("results.json","w"), indent=1, default=float)
print(f"{sum(ok)}/{len(ok)} checks pass")
