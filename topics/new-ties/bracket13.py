"""Round 13 (DRAFT): the certified bracket for the original break with runs of up to K = 10^6 copies of lam.

Same weights as original_bracket.py. Write E_k = W_inf - W_k = L^(k+1)/(1-L) * [(Q-1) + (M-1) nu^s Q] > 0. With s_k, s_F the
roots of W_k = 1 and W_inf = 1, the mean value theorem gives s_F - s_k = E_k(s_k) / |W_inf'(xi)| for some xi in [s_k, s_F].
Both s_k and xi lie in J = [s_F - 1e-35, s_F + 1e-35] (Round 11 put s_10000, hence s_k, within 1.8e-4811 of s_F), so
    s_F - s_k <= sup_J E_k / inf_J |W_inf'|,
evaluated in interval arithmetic at 60 digits (mpmath iv; exponents of any size are exact in mpmath).
Writes bracket13.json. Usage: python3 bracket13.py [K]"""
import json, sys
from mpmath import iv, mp, mpf, log10
K = int(sys.argv[1]) if len(sys.argv) > 1 else 10 ** 6
iv.dps = 60; mp.dps = 60
SF = "0.4704834398630828674262577144580257083077"     # s_F from original_bracket.json (certified to ~4,950 digits)
J = iv.mpf([mpf(SF) - mpf("1e-35"), mpf(SF) + mpf("1e-35")])
lam, mu, nu, kap = iv.mpf(19) / 200, iv.mpf(1) / 20, iv.mpf(1) / 19, iv.mpf(1) / 10
def parts(s):
    L, Mu, N, Kp = iv.exp(s * iv.log(lam)), iv.exp(s * iv.log(mu)), iv.exp(s * iv.log(nu)), iv.exp(s * iv.log(kap))
    dL, dMu, dN, dKp = L * iv.log(lam), Mu * iv.log(mu), N * iv.log(nu), Kp * iv.log(kap)
    Lam = 1 / (1 - L); dLam = dL / (1 - L) ** 2
    M = 1 / (1 - Mu); dM = dMu / (1 - Mu) ** 2
    Q = 1 / (1 - N - Kp); dQ = (dN + dKp) / (1 - N - Kp) ** 2
    W = (Lam - 1) * (Q - 1) + Lam * (M - 1) * N * Q
    dW = dLam * (Q - 1) + (Lam - 1) * dQ + dLam * (M - 1) * N * Q + Lam * dM * N * Q + Lam * (M - 1) * dN * Q + Lam * (M - 1) * N * dQ
    return L, Lam, M, N, Q, W, dW
L, Lam, M, N, Q, W, dW = parts(J)
assert dW.b < 0, "W_inf must be decreasing on J"
assert (W - 1).a <= 0 <= (W - 1).b, "s_F must lie in J"
E = iv.exp((K + 1) * J * iv.log(lam)) * Lam * ((Q - 1) + (M - 1) * N * Q)
gap = mpf((iv.mpf(E.b) / (-iv.mpf(dW.b))).b)      # upper end of the interval quotient
out = dict(k=K, log10_gap_bound=float(log10(gap)), gap_bound=mp.nstr(gap, 6), s_F=SF,
           lower=mp.nstr(mpf(SF) - gap, 40) + " (equal to s_F at 40 digits)", round11_log10_gap=-4810.7516)
json.dump(out, open("bracket13.json", "w"), indent=1); print(out)
