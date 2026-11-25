"""Bonding and valence through mult-add.
1. Overlap: |φa+φb|² = a² + b² + 2ab·S is the invariant-mass law with cosine S.
2. Hückel chains: the transfer recurrence P(n+2) = x P(n+1) − P(n) is the landscape point (k,q) = (x,1);
   band = turning (|x|<2), band edge = border, gap = two walls. Cassini invariant = q^n.
3. Rings: Hückel 4n+2 / Möbius 4n.  4. Bond alternation (polyacetylene) opens a gap = two-walls country.
5. MO bond orders.  6. Bond-energy ladders and the dial (with nitrogen's triple bond).
Run: python3 checks.py"""
import json, math
import numpy as np
from scipy.linalg import eigh
ok = []; R = {}
def check(name, cond, info=""):
    ok.append(bool(cond)); print(("PASS " if cond else "FAIL ") + name + ("  " + info if info else ""))

# 1. overlap law of cosines, with the hydrogen 1s overlap at the H2 bond length
def S1s(r): return math.exp(-r)*(1 + r + r*r/3)
S = S1s(1.4011)
R["S_H2"] = S
rng = np.random.default_rng(1)
# random vectors with prescribed overlap: check |a u + b v|² = a²+b²+2ab S
u = np.array([1.0, 0.0]); v = np.array([S, math.sqrt(1-S*S)])
good = all(abs(np.dot(a*u+b*v, a*u+b*v) - (a*a + b*b + 2*a*b*S)) < 1e-12 for a, b in rng.normal(size=(50, 2)))
check("|a φa + b φb|² = a² + b² + 2ab·S, the same law as M² = m1² + m2² + 2 m1 m2 c", good, f"H2 1s overlap S = {S:.4f}")
check("S = 1 is the border (orbitals identical, lengths add); |S| < 1 is the turning country", abs(S1s(0) - 1) < 1e-15 and 0 < S < 1)

# asymmetric splitting from the generalized eigenproblem H c = E S c
for alpha, s in [(-13.6, S), (-11.4, 0.25), (-10.0, 0.1)]:
    beta = 1.75*s*alpha          # Wolfsberg–Helmholz estimate, as in extended Hückel
    H = np.array([[alpha, beta], [beta, alpha]]); M = np.array([[1, s], [s, 1]])
    E = eigh(H, M, eigvals_only=True)
    Eb, Ea = (alpha+beta)/(1+s), (alpha-beta)/(1-s)
    assert abs(E[0]-Eb) < 1e-9 and abs(E[1]-Ea) < 1e-9, (E, Eb, Ea)
    rise, fall = Ea - alpha, alpha - Eb
    assert abs(rise/fall - (1+s)/(1-s)) < 1e-9
check("bonding (α+β)/(1+S), antibonding (α−β)/(1−S): antibonding rises more by exactly (1+S)/(1−S) (β = 1.75 S α)", True,
      f"H2-like: ratio {(1+S)/(1-S):.2f}")
R["split_ratio"] = (1+S)/(1-S)

# 2. Hückel chains
def chainP(n, x, k=None, q=1.0):
    k = x if k is None else k
    a, b = 1.0, k
    if n == 0: return a
    for _ in range(n-1): a, b = b, k*b - q*a
    return b
allgood = True
for N in range(1, 13):
    A = np.diag(np.ones(N-1), 1) + np.diag(np.ones(N-1), -1)
    ev = np.sort(np.linalg.eigvalsh(A))
    pred = np.sort([2*math.cos(j*math.pi/(N+1)) for j in range(1, N+1)])
    allgood &= np.allclose(ev, pred)
    allgood &= all(abs(chainP(N, x)) < 1e-9 for x in pred)
check("chain of N carbons: levels x_j = 2cos(jπ/(N+1)) (in units of β), the roots of P(N+2)=xP(N+1)−P(N)", allgood, "N = 1…12")
cass = all(abs(chainP(n+1, x)**2 - chainP(n+2, x)*chainP(n, x) - 1) < 1e-8 for n in range(0, 15) for x in np.linspace(-3, 3, 13))
check("Cassini invariant P(n+1)² − P(n+2)P(n) = 1 for every n: the gluing number q = 1", cass)
casq = all(abs(chainP(n+1, 0, k, q)**2 - chainP(n+2, 0, k, q)*chainP(n, 0, k, q) - q**(n+1)) < 1e-6*max(1, abs(q)**(n+1))
           for n in range(0, 12) for k in (-1.5, 0.3, 2.2) for q in (0.5, 1.0, 1.7))
check("general recurrence P(n+2) = k P(n+1) − q P(n) has invariant q^(n+1): the landscape point (k, q)", casq)
# countries of the infinite chain: eigenvalues of transfer matrix [[x,-1],[1,0]]
def mult(x): return np.linalg.eigvals(np.array([[x, -1.0], [1.0, 0.0]]))
inband = all(np.allclose(abs(mult(x)), 1) and abs(np.imag(mult(x)[0])) > 0 for x in np.linspace(-1.99, 1.99, 50))
gap = all(np.all(np.isreal(mult(x))) and abs(abs(mult(x)[0])*abs(mult(x)[1]) - 1) < 1e-12 and max(abs(mult(x))) > 1
          for x in [2.01, 2.5, 3, -2.5, -4])
edge = np.allclose(mult(2.0), [1, 1], atol=1e-6)
check("|x| < 2: transfer eigenvalues on the unit circle (Bloch waves, turning); |x| > 2: real pair λ, 1/λ (decay, two walls); x = ±2: double root (border)",
      inband and gap and edge)
decay = {x: 1/math.acosh(x/2) for x in (2.2, 3.0, 4.0)}
R["decay"] = decay
check("in the gap the state decays over 1/arccosh(|x|/2) sites, the two-walls rapidity", abs(decay[3.0] - 1/math.log((3+math.sqrt(5))/2)) < 1e-12,
      ", ".join(f"x={x}: {d:.2f} sites" for x, d in decay.items()))

# 3. rings: Hückel (periodic) vs Möbius (antiperiodic)
def ring(N, mob=False):
    return sorted([2*math.cos((2*j + (1 if mob else 0))*math.pi/N) for j in range(N)], reverse=True)
def closed_shell_counts(levels):
    # electron counts that exactly fill degenerate groups
    groups = []; 
    for e in levels:
        if groups and abs(groups[-1][0]-e) < 1e-9: groups[-1][1] += 1
        else: groups.append([e, 1])
    c, tot = [], 0
    for e, g in groups:
        if e < 1e-9: break
        tot += 2*g; c.append(tot)
    return c
hk = all(all(n % 4 == 2 for n in closed_shell_counts(ring(N))) for N in range(3, 20))
mb = all(all(n % 4 == 0 for n in closed_shell_counts(ring(N, True))) for N in range(3, 20))
check("Hückel rings close their bonding shells at 4n+2 electrons; Möbius (one twist) rings at 4n", hk and mb, "N = 3…19")
def pi_energy(levels, ne):
    occ = []; left = ne
    for e in levels:
        t = min(2, left); occ.append(t*e); left -= t
        if left == 0: break
    return sum(occ)
benz = pi_energy(ring(6), 6); cbd = pi_energy(ring(4), 4); bd = pi_energy(sorted([2*math.cos(j*math.pi/5) for j in range(1, 5)], reverse=True), 4)
R.update(benzene=benz, cyclobutadiene=cbd, butadiene=bd)
check("benzene π energy 8β = three ethylenes (6β) + 2β of aromatic stabilization", abs(benz-8) < 1e-12)
check("cyclobutadiene 4β = two ethylenes exactly: no stabilization (4n, Möbius-like in a plain ring)", abs(cbd-4) < 1e-12)
check("butadiene 4.472β: delocalization 0.472β", abs(bd - 2*math.sqrt(5)) < 1e-12, f"{bd:.4f}")

# 4. bond alternation: two-step transfer matrix, still det 1
def alt_gap(t1, t2):
    xs = np.linspace(-3, 3, 60001)
    tr = (xs**2 - t1**2 - t2**2)/(t1*t2)        # trace of two-step transfer matrix
    inb = np.abs(tr) <= 2
    xin = xs[inb]
    return xin[xin > 0].min(), xin.max()
for t1, t2 in [(1.1, 0.9), (1.2, 0.8), (1.0, 1.0)]:
    lo, hi = alt_gap(t1, t2)
    assert abs(lo - abs(t1-t2)) < 2e-4 and abs(hi - (t1+t2)) < 2e-4
check("alternating bonds t1, t2 (polyacetylene): band |t1−t2| ≤ |x| ≤ t1+t2, gap 2|t1−t2| inside the two-walls country", True)
N = 200; t1, t2 = 1.1, 0.9
A = np.zeros((N, N))
for i in range(N-1): A[i, i+1] = A[i+1, i] = t1 if i % 2 == 0 else t2
ev = np.linalg.eigvalsh(A)
R["poly_gap_finite"] = float(np.min(ev[ev > 1e-6]) - np.max(ev[ev < -1e-6])) if True else 0
check("finite 200-site alternating chain shows the same gap 2|t1−t2| = 0.4", abs(R["poly_gap_finite"] - 0.4) < 0.03, f"{R['poly_gap_finite']:.3f}")

# 5. MO bond orders (second row homonuclear, s-p mixed order for B2,C2,N2)
def bond_order(ne):
    order = [("s2s", 1), ("s2s*", -1), ("p2p", 1), ("p2p", 1), ("s2p", 1), ("p2p*", -1), ("p2p*", -1), ("s2p*", -1)]
    b = a = 0; left = ne
    # fill π pair together (Hund) — counts only matter for bond order
    for name, sgn in order:
        t = min(2, left); left -= t
        if sgn > 0: b += t
        else: a += t
    return (b - a)/2
mols = {"Li2": 2, "Be2": 4, "B2": 6, "C2": 8, "N2": 10, "O2": 12, "F2": 14, "Ne2": 16}
bo = {m: bond_order(n) for m, n in mols.items()}
R["bond_orders"] = bo
check("valence MO bond orders: Li2 1, Be2 0, B2 1, C2 2, N2 3, O2 2, F2 1, Ne2 0",
      [bo[m] for m in mols] == [1, 0, 1, 2, 3, 2, 1, 0], str(bo))

# 6. bond-energy ladders and the dial (LibreTexts appendix 3, kJ/mol)
lad = {"C–C": (348, 614, 839), "N–N": (163, 418, 941), "C–N": (293, 615, 891), "C–O": (358, 799, 1072)}
fits = {}
for name, (e1, e2, e3) in lad.items():
    k = (e2 - 2*e1)/(e1*e1)               # dial from the double bond
    pred = ((1 + k*e1)**3 - 1)/k           # associative dial: triple = e1 ✳_k e1 ✳_k e1
    fits[name] = dict(k=k, pred=pred, meas=e3, err=(pred - e3)/e3, inc=(e1, e2-e1, e3-e2))
R["ladders"] = fits
for n, f in fits.items():
    print(f"   {n}: increments {f['inc']}, k = {f['k']:.2e} /(kJ/mol), dial predicts triple {f['pred']:.0f}, measured {f['meas']} ({100*f['err']:+.0f}%)")
check("the dial fixed by single and double bonds misses the triple bond by −13% for N–N, +9% for C–N, +25% for C–O",
      abs(fits['N–N']['err']) > 0.1 and abs(fits['C–N']['err']) > 0.05 and abs(fits['C–O']['err']) > 0.2)
check("only the C–C ladder lands within 3%, so one dial does not describe bond ladders", abs(fits["C–C"]["err"]) < 0.03)
check("nitrogen's ladder accelerates (163, 255, 523 kJ/mol): the third bond is worth 3.2 single bonds; carbon's slows (348, 266, 225)",
      fits["N–N"]["inc"][2] > fits["N–N"]["inc"][1] > fits["N–N"]["inc"][0] and fits["C–C"]["inc"][0] > fits["C–C"]["inc"][1] > fits["C–C"]["inc"][2])
R["N_first_step_share"] = (941-418)/941
check("breaking N≡N to N=N costs 523 of 941 kJ/mol (56%): the first step is the hard one in nitrogen fixation", abs(R["N_first_step_share"] - 0.556) < 0.001)

print(f"\n{sum(ok)}/{len(ok)} checks pass")
json.dump(R, open("results.json", "w"), indent=1, default=float)
