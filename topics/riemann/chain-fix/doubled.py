"""Doubled chain: add a compass spin tau. Energy H(tau, s) = H_even(s) + tau*H_odd(s). Symmetric under flipping every spin,
ferromagnetic, and Z_two(beta) = 2 Z_Knauf(beta) exactly. Field polynomial Z_two(beta, w) = P(w) + w^(n+1) P(1/w)."""
import numpy as np, sys, json, mpmath as mp
sys.path.insert(0, "../knauf-field")
from knauf import chain
mp.mp.dps = 50
def P_coef(n, beta, phi, up):
    w = phi**(-beta); return [mp.fsum([mp.mpf(x) for x in w[up == k]]) for k in range(n+1)]
def minroot(n, beta, phi, up):
    c = P_coef(n, beta, phi, up); r = mp.polyroots(c[::-1], maxsteps=500, extraprec=300)
    return min(abs(z) for z in r)
def two_roots(n, beta, phi, up):
    c = P_coef(n, beta, phi, up); full = [mp.mpf(0)]*(n+2)
    for k in range(n+1): full[k] += c[k]; full[n+1-k] += c[k]
    return mp.polyroots(full[::-1], maxsteps=500, extraprec=300)
def main():
  out = {}
  for n in range(4, 23):
      phi = chain(n).astype(float); up = np.array([bin(k).count("1") for k in range(2**n)])
      if n <= 12:
          r = two_roots(n, 1.5, phi, up); dev = max(abs(abs(z)-1) for z in r)
          print("   n=%d beta=1.5 doubled-chain field zeros: max ||w|-1| = %s, min|root P| = %s" % (n, mp.nstr(dev, 3), mp.nstr(minroot(n, 1.5, phi, up), 6)))
      lo, hi = 0.05, 4.0
      for _ in range(36):
          mid = 0.5*(lo+hi)
          if minroot(n, mid, phi, up) > 1: hi = mid
          else: lo = mid
      out[n] = hi; print(n, "P zero-free in the closed unit disk (so doubled chain is Lee-Yang) for beta >= %.6f" % hi, flush=True)
      json.dump(out, open("doubled_threshold.json", "w"))
if __name__ == '__main__':
  main()
