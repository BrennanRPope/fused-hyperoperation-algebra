"""Checks for section 6 of the dial names (Nov 26, 2026): the odds t, similarity ratio r, growth degree beta and holonomy (twist).
Run: python3 checks_dials6.py > checks_dials6_output.txt"""
import numpy as np, sympy as sp, mpmath as mp
res = []
def check(n, ok, detail=""):
    res.append(bool(ok)); print(("PASS " if ok else "FAIL ") + n + (("  | " + detail) if detail else ""))
circ = lambda a, b, k: a + b + k*a*b
t, s, w, a, b, c, e, k = sp.symbols('t s w a b c e k')

# odds t (base dial, one per element); the chance t/(1+t) is a reading
ch = lambda x: x/(1 + x)
check("T1 independent union: chance 1-(1-p)(1-q) has odds t o_1 s, the homothety law at k = 1",
      sp.simplify(sp.simplify((1 - (1 - ch(t))*(1 - ch(s))))/(1 - (1 - (1 - ch(t))*(1 - ch(s)))) - circ(t, s, 1)) == 0)
check("T2 complement is t -> 1/t, the chance goes to 1 - chance; with t = tan(theta) this is the mirror in the diagonal theta -> 90 deg - theta",
      sp.simplify(ch(1/t) - (1 - ch(t))) == 0 and abs(np.tan(np.pi/2 - 0.3) - 1/np.tan(0.3)) < 1e-12)
th = np.arctan(np.sqrt(2) - 1)
check("T3 the union half of the full set is t = sqrt2 - 1 = tan 22.5 deg: halves add to 45 deg", abs(circ(np.sqrt(2) - 1, np.sqrt(2) - 1, 1) - 1) < 1e-12 and abs(np.degrees(th) - 22.5) < 1e-12)
bet, eps, mu = 1.7, 0.9, 0.3
n = 1/(np.exp(bet*(eps - mu)) + 1)
check("T4 Fermi-Dirac: a mode's odds of being filled n/(1-n) = z e^(-beta eps), fugacity times Boltzmann factor (Lewis's fugacity)",
      abs(n/(1 - n) - np.exp(bet*mu)*np.exp(-bet*eps)) < 1e-12)

# similarity ratio r (base dial, one per piece); dimension D is a reading (Moran)
def moran(rs):
    return float(mp.findroot(lambda x: sum(mp.mpf(r)**x for r in rs) - 1, 1.0))
D = [moran([0.25, 0.25]), moran([1/3, 1/3]), moran([0.5]*3)]
check("R1 Moran's equation sum r_i^D = 1 gives D as an output: half box 1/2, Cantor log2/log3, gasket log3/log2",
      abs(D[0] - 0.5) < 1e-12 and abs(D[1] - np.log(2)/np.log(3)) < 1e-12 and abs(D[2] - np.log(3)/np.log(2)) < 1e-12,
      f"{D[0]:.6f}, {D[1]:.6f}, {D[2]:.6f}")
check("R2 ratios compose by multiplying (a piece of a piece has ratio r1 r2) and the weight r^D multiplies with them",
      abs((0.25*0.5)**D[0] - 0.25**D[0]*0.5**D[0]) < 1e-15)
check("R3 a weight w acts as the homothety dial: w (a o_w b) = (w a) o_1 (w b), the library's wAdd_scale",
      sp.expand(w*circ(a, b, w) - circ(w*a, w*b, 1)) == 0)
check("R4 the size measure's doubling weight 2^(-beta) is ratio 1/2 at dimension D = beta", all(abs(0.5**bb - 2.0**(-bb)) < 1e-15 for bb in (0, 1, 2, 3)))

# growth degree beta (a setting of the size measure)
vals = [float(mp.polylog(bb + 1, -0.5)/mp.polylog(bb + 1, -1)) for bb in (0, 1, 2, 3)]
check("B1 plain half keeps Li_(beta+1)(-1/2)/Li_(beta+1)(-1): 58.50, 54.52, 52.42, 51.29 % for beta = 0..3",
      all(abs(100*v - x) < 0.006 for v, x in zip(vals, (58.50, 54.52, 52.42, 51.29))), ", ".join(f"{100*v:.2f}" for v in vals))
check("B2 the repetition share 2^(-beta) is 1 - eta(beta+1)/zeta(beta+1): 7/8 Stefan-Boltzmann at beta = 3",
      all(abs((1 - mp.altzeta(bb + 1)/mp.zeta(bb + 1)) - 2.0**(-bb)) < 1e-12 for bb in (1, 2, 3)) and abs(mp.altzeta(4)/mp.zeta(4) - 7/8) < 1e-12)

# holonomy (twist): a reading of the curvature dial and the directions; it measures the oriented area a ^ e - b ^ c
comp = lambda f, g: (circ(f[0], g[0], 1), f[1] + (1 + f[0])*g[1])     # (a, c): x -> (1 + a) x + c ; f after g
fg, gf = comp((a, c), (b, e)), comp((b, e), (a, c))
check("W1 two stretch-and-shift maps (a, c), (b, e): stretches commute; the shifts differ by a e - b c, the oriented area (a wedge)",
      sp.expand(fg[0] - gf[0]) == 0 and sp.expand(fg[1] - gf[1] - (a*e - b*c)) == 0)
q, u1, u2, v1, v2 = 0.37, 0.21, -0.12, 0.08, 0.33
U, V = complex(u1, u2), complex(v1, v2)
X, Y = 1 - q*(u1*v1 + u2*v2), q*(u1*v2 - u2*v1)
om = 2*np.arctan2(Y, X)
lhs = (U + V)/(1 - q*np.conj(U)*V); rhs = (V + U)/(1 - q*np.conj(V)*U)
check("W2 the twist is the phase between the two orders of the dial boost law (u + v)/(1 - q conj(u) v), q = d at k = 0: "
      "tan(omega/2) = q (u x v)/(1 - q u.v)", abs(lhs - np.exp(1j*om)*rhs) < 1e-12 or abs(rhs - np.exp(1j*om)*lhs) < 1e-12,
      f"omega = {om:.6f}")
print(f"\n{sum(res)} of {len(res)} pass")
