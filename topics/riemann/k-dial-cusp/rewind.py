"""With the k dial on, the chain sum is F_k(b) = sum_m phi(m) [m]_q^(-b), q = 1+k: it converges for every Re b > 0
(geometric tail, no cusp). Exactly, with L = log(1+k) and h(x) = (x/(e^x - 1))^b (Bose/Bernoulli),
  (L/k)^b F_k(b) = sum_m phi(m) m^(-b) h(L m),
a smooth cutoff. Mellin: = F(b) + H(2-b) L^(b-2)/zeta(2) + O(L) on Re b = 3/2, H = Mellin transform of h.
Rewind k -> 0 subtracting the one known term, with no fitting, and compare with F = zeta(b-1)/zeta(b). DRAFT."""
import numpy as np, mpmath as mp, json, sys
def totients(M):
    ph = np.arange(M+1, dtype=np.int64)
    for p in range(2, M+1):
        if ph[p] == p: ph[p::p] -= ph[p::p]//p
    return ph
PH = totients(400000)
def Fk(b, k):
    """Returns (L/k)^b * sum phi(m) [m]_q^(-b) = sum phi(m) m^(-b) h(L m), and L."""
    L = np.log1p(k); M = int(min(400000, 40/(1.5*L) + 50)); m = np.arange(1, M+1); x = L*m
    h = np.exp(b*(np.log(x) - np.log(np.expm1(x))))
    return complex((PH[1:M+1]*np.exp(-b*np.log(m))*h).sum()), L
def H(s, b): return complex(mp.quad(lambda x: x**(s-1)*(x/mp.expm1(x))**b, [0, 1, 10, mp.inf]))
def G(b, k, Hv):
    f, L = Fk(b, k); return f - Hv*L**(b-2)/float(mp.zeta(2)), f
if __name__ == "__main__":
    ks = [0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001, 0.0005, 0.0002]; out = dict(ks=ks, rows={})
    for t in [5.0, 10.0, 14.134725141734695, 18.0, 21.022039638771556, 25.01085758014569, 30.0, 40.0]:
        b = 1.5+1j*t; Hv = H(2-b, b); tr = complex(mp.zeta(b-1)/mp.zeta(b)); row = []
        for k in ks:
            g, f = G(b, k, Hv); row.append([abs(g-tr), abs(f - tr)])
        out["rows"]["%.4f" % t] = row
        print("t %6.2f  rewind error: %s" % (t, " ".join("%.1e" % r[0] for r in row)))
        print("          raw (no subtraction): %s" % " ".join("%.1e" % r[1] for r in row))
    json.dump(out, open("rewind.json", "w"), indent=1)
