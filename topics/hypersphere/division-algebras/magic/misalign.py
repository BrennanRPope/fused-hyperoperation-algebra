"""Probe bodies with mixed charge directions orbiting four-charge black holes (STU, verified exact in stu.py).
Probe = components j (weights w_j, signs s_j) riding one worldline: m(phi) = sum w_j e^{-a_j.phi/2},
q A_t = sum w_j s_j (1 - 1/H_j). Aligned (s=+1) components feel no static force from their own charge."""
import sympy as sp, itertools, json
from stu import AV
from precess import T_series, T_numeric, u, L, E
def W_stu(ks, ws, ss):
    Hs = [1 + k*u for k in ks]; H = sp.prod(Hs)
    phi = sum((AV[i]*sp.log(Hs[i]) for i in range(4)), sp.zeros(3, 1))/2
    m = sum(w*sp.exp(-(AV[j].T*phi)[0]/2) for j, w in enumerate(ws))
    qA = sum(w*s*(1 - 1/Hs[j]) for j, (w, s) in enumerate(zip(ws, ss)))
    A = H**sp.Rational(-1, 2); B = H**sp.Rational(1, 2); C = H**sp.Rational(1, 2)
    return sp.simplify((C**2/B)*((E - qA)**2/A - m**2 - L**2*u**2/C))
if __name__ == "__main__":
    one = sp.Integer(1); rows = []
    bgs = {"1 charge": [one, 0, 0, 0], "2 charges": [one, one, 0, 0], "3 charges": [one, one, one, 0], "4 charges": [one, one, one, one]}
    for bname, ks in bgs.items():
        for pat in itertools.product([0, 1, -1], repeat=4):
            if all(p == 0 for p in pat): continue
            n = sum(1 for p in pat if p); ws = [sp.Rational(abs(p), n) for p in pat]; ss = [p if p else 1 for p in pat]
            P = sp.Matrix([w*s for w, s in zip(ws, ss)]); Q = sp.Matrix(ks)
            # static pull: only anti-aligned components on charged directions attract
            if not any(p == -1 and ks[j] != 0 for j, p in enumerate(pat)): continue
            W = W_stu(ks, ws, ss)
            try:
                T, K = T_series(W)
            except Exception as ex:
                T, K = "err", str(ex)
            rows.append(dict(bg=bname, probe=pat, T=str(T), K=str(K), cos=str(sp.nsimplify((Q.dot(P))/(Q.norm()*P.norm())))))
            print(bname, pat, "T =", T, " K =", K)
    json.dump(rows, open("results_misalign_grid.json", "w"), indent=1)
