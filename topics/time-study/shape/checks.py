"""DRAFT shape checks. Units c = 1, linear order. Brennan: "Derive the shape". Can the chain give the four-leaf
(transverse and balanced) for a wave, without the field equations? If not, what exactly is the missing input?
Input from the chain: the pattern moves at the speed limit c along z (../tilt S1-S3), so every metric piece is a function
of u = t - z only. Nothing else is assumed: all ten metric pieces are free."""
import json, itertools, numpy as np, sympy as sp
R = {}
def ok(name, cond, info=""):
    R[name] = {"pass": bool(cond), "info": info}; print(("PASS " if cond else "FAIL ") + name, info)

idx = [(a, b) for a in range(4) for b in range(a, 4)]
H = {ab: sp.Symbol(f"H{ab[0]}{ab[1]}") for ab in idx}                 # second u-derivatives of the ten metric pieces
Hm = sp.Matrix(4, 4, lambda a, b: H[(min(a, b), max(a, b))])
k = [1, 0, 0, -1]                                                     # d/dx^mu of f(t - z) = k_mu f'
eta = sp.diag(-1, 1, 1, 1)
def Rm(a, b, c, d):   # linear Riemann for h(u): 1/2 (k_b k_c H_ad + k_a k_d H_bc - k_a k_c H_bd - k_b k_d H_ac)
    return sp.Rational(1, 2)*(k[b]*k[c]*Hm[a, d] + k[a]*k[d]*Hm[b, c] - k[a]*k[c]*Hm[b, d] - k[b]*k[d]*Hm[a, c])
Hs = [H[ab] for ab in idx]
comps = [Rm(a, b, c, d) for a, b, c, d in itertools.product(range(4), repeat=4)]
J = sp.Matrix([[sp.diff(e, h) for h in Hs] for e in comps])
rank_R = J.rank()
ok("D1 the speed limit alone cuts the curvature to 6 independent pieces (out of 20)", rank_R == 6, f"rank {rank_R}")

Ric = sp.Matrix(4, 4, lambda b, d: sum(eta[a, c]*Rm(a, b, c, d) for a in range(4) for c in range(4)))
Jr = sp.Matrix([[sp.diff(Ric[b, d], h) for h in Hs] for b in range(4) for d in range(b, 4)])
JJ = sp.Matrix.vstack(J, Jr)
ricci_rank = Jr.rank(); weyl_left = rank_R - ricci_rank
ok("D2 of those 6, 4 are a source riding along (Ricci); 2 carry nothing with them", ricci_rank == 4 and weyl_left == 2,
   f"Ricci rank {ricci_rank}, source-free pieces {weyl_left}")

E = sp.Matrix(3, 3, lambda i, j: Rm(0, i+1, 0, j+1))
# D3: impose 'nothing rides along' (Ricci = 0) and read off E: transverse and balanced, i.e. the four-leaf
sol = sp.solve(list(Ric), Hs, dict=True)[0]
E0 = sp.simplify(E.subs(sol))
ok("D3 with nothing riding along, E is exactly the four-leaf: E_zj = 0 and E_xx = -E_yy",
   all(sp.simplify(E0[2, j]) == 0 for j in range(3)) and sp.simplify(E0[0, 0] + E0[1, 1]) == 0 and E0 != sp.zeros(3), f"E = {E0.tolist()}")

# D4: 'balanced' for one observer only (E traceless, R_00 = 0) is not enough: longitudinal pieces survive.
sol1 = sp.solve([sp.expand(E.trace())], Hs, dict=True)[0]
E1 = sp.simplify(E.subs(sol1))
r1 = sp.Matrix([[sp.diff(Ric.subs(sol1)[b, d], h) for h in Hs] for b in range(4) for d in range(b, 4)]).rank()
ok("D4 balanced for one observer leaves 3 riding-along pieces (E_xz, E_yz can survive)", r1 == 3 and (E1[0, 2] != 0 or E1[1, 2] != 0),
   f"Ricci rank left {r1}; E_xz = {E1[0,2]}")

# D5: balanced for EVERY observer (R(u,u) = 0 for all moving observers u) forces Ricci = 0, hence the four-leaf.
rng = np.random.default_rng(5); eqs = []
for _ in range(10):
    v = rng.normal(size=3); v = 0.6*v/np.linalg.norm(v)*rng.random(); g = 1/np.sqrt(1 - v@v)
    uvec = [g] + list(g*v)
    eqs.append(sp.expand(sum(sp.nsimplify(round(uvec[a], 6))*sp.nsimplify(round(uvec[b], 6))*Ric[a, b] for a in range(4) for b in range(4))))
Ja = sp.Matrix([[sp.diff(e, h) for h in Hs] for e in eqs])
left = sp.Matrix.vstack(Ja, J).rank() - Ja.rank()   # curvature pieces not fixed to zero by the observer conditions
ok("D5 balanced for every observer leaves exactly the 2 four-leaf pieces", left == 2, f"curvature pieces left {left}")

# D6: helicity. Turning by psi about z: the four-leaf turns by 2 psi (spin 2); the riding-along pieces show up in E as
#     E_xz, E_yz (turn by psi, spin 1) and E_zz, E_xx + E_yy (unchanged, spin 0).
psi = sp.Rational(1, 7); Rz = sp.Matrix([[sp.cos(psi), -sp.sin(psi), 0], [sp.sin(psi), sp.cos(psi), 0], [0, 0, 1]])
a_, b_, c_, d_, e_, f_ = sp.symbols('a b c d e f')
Eg = sp.Matrix([[a_ + b_, c_, d_], [c_, a_ - b_, e_], [d_, e_, f_]])
Er = Rz*Eg*Rz.T
spin2 = sp.simplify((Er[0, 0] - Er[1, 1])/2 - sp.I*Er[0, 1] - sp.exp(-2*sp.I*psi)*(b_ - sp.I*c_))
spin1 = sp.simplify(Er[0, 2] - sp.I*Er[1, 2] - sp.exp(-sp.I*psi)*(d_ - sp.I*e_))
spin0 = sp.simplify((Er[0, 0] + Er[1, 1]) - 2*a_) == 0 and sp.simplify(Er[2, 2] - f_) == 0
ok("D6 the four-leaf is the spin-2 part; what rides along is spin 1 and spin 0",
   sp.simplify(sp.expand_complex(spin2)) == 0 and sp.simplify(sp.expand_complex(spin1)) == 0 and spin0)

# D7 (chain reading of 'balanced'): his neighbour chain from ../tilt S1 averages each place with its neighbours. Its rest
#     states are the balanced patterns (each place = average of its neighbours). In one frame that keeps the four-leaf
#     unchanged but moves the trace pattern x^2 + y^2: it fixes only the trace, matching D4.
g1 = np.arange(-6, 7); Xg, Yg, Zg = np.meshgrid(g1, g1, g1, indexing="ij")
def avg_step(F):   # each place becomes the average of itself and its six neighbours (interior only)
    G_ = F.copy(); G_[1:-1, 1:-1, 1:-1] = (F[1:-1, 1:-1, 1:-1] + F[2:, 1:-1, 1:-1] + F[:-2, 1:-1, 1:-1] + F[1:-1, 2:, 1:-1]
                     + F[1:-1, :-2, 1:-1] + F[1:-1, 1:-1, 2:] + F[1:-1, 1:-1, :-2])/7; return G_
leaf = (Xg**2 - Yg**2).astype(float); spin1 = (Xg*Zg).astype(float); trace = (Xg**2 + Yg**2).astype(float)
chg = {n: float(np.abs(avg_step(F) - F)[1:-1, 1:-1, 1:-1].max()) for n, F in [("four-leaf", leaf), ("spin 1", spin1), ("trace", trace)]}
ok("D7 neighbour averaging leaves the four-leaf (and spin 1) at rest but moves the trace: one frame fixes only the trace",
   chg["four-leaf"] < 1e-12 and chg["spin 1"] < 1e-12 and chg["trace"] > 0.5, str(chg))

json.dump({"checks": R, "funnel": {"all curvature": 20, "moving at c": int(rank_R), "balanced for one observer": int(r1) + 2,
           "balanced for every observer": int(left)}, "avg_change": chg}, open(__file__.replace("checks.py", "results.json"), "w"), indent=1)
print(f"{sum(r['pass'] for r in R.values())} of {len(R)} pass")
