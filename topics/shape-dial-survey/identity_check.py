"""The step identity linking a and a+1 for the Lee-Yang field function, and what it does.
G_a(s) = int_0^1 (1 - x^p)^a cos(s x) dx,  theta = s d/ds,  kappa = (a+1) p + 1.
Identity:  (a+1) p G_a  =  (theta + kappa) G_(a+1)  =  s^(1-kappa) d/ds [ s^kappa G_(a+1) ].
"""
import json, time
import numpy as np
import mpmath as mp
exec(open("interlace_check.py").read().split("if __name__")[0])     # coeffs, ev, zeros, interlace, S, DS
RES = []; OUT = {}
def check(name, ok, detail=""):
    RES.append(bool(ok)); print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))
t0 = time.perf_counter()
# K1 identity on the functions, by quadrature (independent of the series)
mp.mp.dps = 30; worst = 0
for p, a, s in [(3.0, 1.5, 4.2), (1.5, 1.2, 9.1), (5.0, 3.3, 2.7), (2.0, 2.0, 6.0), (4.0, 0.4, 11.3)]:
    G = lambda t, a=a: mp.quad(lambda x: (1 - x**p)**a*mp.cos(t*x), [0, 1])
    kap = (a + 1)*p + 1
    lhs = (a + 1)*p*G(s, a); rhs = s*mp.diff(lambda t: G(t, a + 1), s) + kap*G(s, a + 1)
    worst = max(worst, abs(lhs - rhs))
check("K1 (a+1) p G_a = (s d/ds + (a+1) p + 1) G_(a+1), by quadrature at 5 (p, a, s)", worst < 1e-12, f"max gap {float(worst):.1e}")
# K2 the same identity coefficient by coefficient (exact Beta-function step)
mp.mp.dps = 50; worst = 0
for p, a in [(3.0, 1.5), (1.5, 1.2), (6.0, 4.0), (2.5, 0.0)]:
    mp.mp.dps = 50; P, A = mp.mpf(p), mp.mpf(a)
    c0 = coeffs(P, A, K=40, dps=50); c1 = coeffs(P, A + 1, K=40, dps=50); kap = (A + 1)*P + 1
    worst = max(worst, max(abs((A + 1)*P*c0[k] - (2*k + kap)*c1[k])/abs(c0[k]) for k in range(40)))
check("K2 coefficient form: (a+1) p c_k(a) = (2k + kappa) c_k(a+1) for k < 40 (relative)", worst < 1e-40, f"max {float(worst):.1e}")
# K3 Rolle: every gap (0, z1), (z1, z2), ... between zeros of G_(a+1) holds exactly one zero of G_a
il = json.load(open("interlace_check.json"))["interlace"]; n_ok = n = 0
for key, v in il.items():
    if v["n_b"] == 0: continue
    p = float(key.split(",")[0][2:]); a = float(key.split("a=")[1])
    za, _ = zeros(p, a, polish=True); zb, _ = zeros(p, a + 1, polish=True)
    edges = [0.0] + zb; good = all(sum(1 for z in za if edges[i] < z < edges[i + 1]) == 1 for i in range(len(zb)))
    n += 1; n_ok += good
check("K3 Rolle: each gap between zeros of G_(a+1) holds exactly one zero of G_a", n_ok == n, f"{n_ok} of {n} pairs")
# K4 why the step can't simply go up: F(s) = int_0^s t^(kappa-1) G_a(t) dt is s^kappa G_(a+1)/((a+1) p).
# Its values at the zeros of G_a must alternate in sign for G_(a+1) to keep a zero in every gap.
lobes = {}
for p, a in [(3.0, 1.0), (1.5, 1.0), (1.5, 0.4)]:
    kap = (a + 1)*p + 1; za, ca = zeros(p, a, polish=True); cb = coeffs(p, a + 1)
    Fz = [float(mp.mpf(z)**kap*ev(cb, z)/((a + 1)*p)) for z in za]        # F at the extrema, via the identity
    alt = [Fz[i]*Fz[i + 1] < 0 for i in range(len(Fz) - 1)]
    first_fail = next((i + 1 for i, ok_ in enumerate(alt) if not ok_), None)
    zb, _ = zeros(p, a + 1)
    lobes[f"p={p:g}, a={a:g} -> {a+1:g}"] = dict(zeros_a=[round(z, 3) for z in za], F_at_zeros=Fz, alternates=all(alt),
                                               first_fail=first_fail, n_zeros_next=len(zb), above_p=bool(a + 1 > p))
    # cross-check F against direct quadrature of t^(kappa-1) G_a at the first extremum
    if p == 3.0:
        mp.mp.dps = 30; z1 = za[0]
        direct = mp.quad(lambda t: t**(kap - 1)*ev(ca, t), [0, z1])
        OUT["K4_direct_vs_identity"] = float(abs(direct - Fz[0])/abs(Fz[0]))
OUT["lobes"] = lobes
print("    ", {k: (v["alternates"], v["first_fail"], v["n_zeros_next"]) for k, v in lobes.items()})
check("K4 the up-step works exactly where the rule says: values alternate for p = 3 (a+1 = 2 <= p) and p = 1.5, a = 0.4 -> 1.4, "
      "and stop alternating for p = 1.5, a = 1 -> 2 (a+1 > p)",
      lobes["p=3, a=1 -> 2"]["alternates"] and lobes["p=1.5, a=0.4 -> 1.4"]["alternates"]
      and not lobes["p=1.5, a=1 -> 2"]["alternates"] and OUT["K4_direct_vs_identity"] < 1e-10,
      f"identity vs direct integral {OUT['K4_direct_vs_identity']:.1e}")
OUT["seconds"] = round(time.perf_counter() - t0, 1)
json.dump(OUT, open("identity_check.json", "w"), indent=1, default=float)
print(f"\n{sum(RES)} of {len(RES)} checks pass  ({OUT['seconds']} s)")
