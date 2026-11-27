"""Chain-fix checks. DRAFT. Each line prints PASS/FAIL."""
import numpy as np, json, sys, mpmath as mp
sys.path.insert(0, "../knauf-field"); from knauf import chain
mp.mp.dps = 50; res = []
def ok(name, cond): res.append(bool(cond)); print(("PASS " if cond else "FAIL ") + name)
def walsh(f):
    a = f.astype(float).copy(); h = 1
    while h < len(a):
        for i in range(0, len(a), 2*h):
            x = a[i:i+h].copy(); y = a[i+h:i+2*h].copy(); a[i:i+h] = x+y; a[i+h:i+2*h] = x-y
        h *= 2
    return a/len(a)
up_of = lambda n: np.array([bin(k).count("1") for k in range(2**n)])
# C1 Knauf's chain is ferromagnetic (known: Knauf 1993; Contucci-Knauf)
ok("C1 every Walsh link of -ln(phi) is >= 0 (n = 6, 8, 10)", all((-walsh(np.log(chain(n).astype(float))))[1:].min() > -1e-12 for n in (6, 8, 10)))
# C2 but it is lopsided under flipping every spin
ok("C2 phi(all down) = 1 and phi(all up) = n+1 (n = 4..16)", all(chain(n)[0] == 1 and chain(n)[-1] == n+1 for n in range(4, 17)))
# doubled chain on spins (tau, s): H2 = He(s) + tau*Ho(s)
n = 10; phi = chain(n).astype(float); H = np.log(phi); He = 0.5*(H+H[::-1]); Ho = 0.5*(H-H[::-1])
H2 = np.concatenate([He+Ho, He-Ho])        # first half tau = +1 (bit 0), second half tau = -1 (bit 1)
for b in (1.5, 2.5, 3.0):
    ok("C3 doubled chain Z = 2 x Knauf Z exactly (n = 10, beta = %g)" % b, abs(np.exp(-b*H2).sum() - 2*np.exp(-b*H).sum()) < 1e-10*np.exp(-b*H).sum())
ok("C4 doubled chain unchanged by flipping every spin, compass included", np.allclose(H2, H2[::-1]))
J2 = walsh(H2); ok("C5 doubled chain is still ferromagnetic (all 2^11 Walsh links of -H2 >= 0)", (-J2)[1:].min() > -1e-12)
ok("C6 doubled chain has no odd-order links left (no built-in field)", np.abs(J2[up_of(n+1) % 2 == 1]).max() < 1e-12)
# C7 field polynomial of doubled chain = P(w) + w^(n+1) P(1/w)
b = 1.5; u = up_of(n); u2 = np.concatenate([u, u+1])   # tau down (bit 1) adds one up-spin in w-counting convention
direct = np.bincount(u2, weights=np.exp(-b*H2), minlength=n+2)
P = np.bincount(u, weights=phi**(-b), minlength=n+1); Q = np.zeros(n+2); Q[:n+1] += P; Q[1:] += P[::-1]
ok("C7 doubled field polynomial = P(w) + w^(n+1) P(1/w)", np.allclose(direct, Q))
# C8 self-inversive (Cohn type) lemma at work
def Pc(n, b):
    p = chain(n).astype(float); u = up_of(n); return [mp.fsum([mp.mpf(x) for x in p[u == k]**(-b)]) for k in range(n+1)]
def Qroots(n, b):
    c = Pc(n, b); f = [mp.mpf(0)]*(n+2)
    for k in range(n+1): f[k] += c[k]; f[n+1-k] += c[k]
    return mp.polyroots(f[::-1], maxsteps=500, extraprec=300)
c = Pc(12, 1.5); mr = min(abs(z) for z in mp.polyroots(c[::-1], maxsteps=500, extraprec=300))
ok("C8 n = 12, beta = 1.5: P has no zeros in |w| <= 1 (min |root| = %s) and every Q zero is on |w| = 1" % mp.nstr(mr, 5), mr > 1 and max(abs(abs(z)-1) for z in Qroots(12, 1.5)) < 1e-30)
thr = json.load(open("doubled_threshold.json")); dd = json.load(open("doubled_direct.json"))
ok("C9 Cohn threshold is close to where the circle really fails (within 0.07, n = 6..14)", all(0 <= thr[str(k)] - dd[str(k)]["last_fail"] < 0.07 for k in range(6, 15)))
ok("C10 n = 12: zeros leave the circle just below the true edge (beta = 1.14)", max(abs(abs(z)-1) for z in Qroots(12, 1.14)) > 1e-6)
ts = json.load(open("threshold.json")); ts.update(json.load(open("threshold_hi.json")))
for name, d in (("flip-averaged", ts), ("doubled", thr)):
    v = [d[k] for k in sorted(d, key=int)]; ok("C11 %s threshold rises with every added spin" % name, all(np.diff(v) > 0))
ok("C12 doubled chain passes the Cohn test at beta = 3/2 up to n = 20 (threshold %.4f), not from n = 21 (%.4f)" % (thr["20"], thr["21"]), thr["20"] < 1.5 < thr["21"] < thr["22"])
# C13 flip-averaged chain: circle exact at n = 12, beta = 1.5 but its Z drains
p = chain(12).astype(float); u = up_of(12); Hs = 0.5*(np.log(p)+np.log(p[::-1]))
r = np.roots(np.bincount(u, weights=np.exp(-1.5*Hs), minlength=13)[::-1])
ok("C13 flip-averaged n = 12, beta = 1.5: all field zeros on the circle (dev %.1e)" % np.abs(np.abs(r)-1).max(), np.abs(np.abs(r)-1).max() < 1e-8)
zs = [((chain(k).astype(float)*chain(k).astype(float)[::-1])**(-1.5)).sum() for k in (10, 16, 22)]
ok("C14 flip-averaged Z at beta = 3 keeps falling: %s" % ["%.4f" % z for z in zs], zs[0] > zs[1] > zs[2])
# C15 per-spin smoothing for +-1 spins is just a field
g = np.array([0.37, 2.9]); h = 0.5*np.log(g[0]/g[1]); cst = np.sqrt(g[0]*g[1])
ok("C15 any weight g(s) on s = +-1 equals c*exp(h*s): a field, nothing smoother", np.allclose(g, cst*np.exp(h*np.array([1, -1]))))
# C16 Knauf limit check
zo = (chain(22).astype(float)**(-3.0)).sum(); tgt = float(mp.zeta(2)/mp.zeta(3))
ok("C16 Knauf Z_22(3) = %.5f approaches zeta(2)/zeta(3) = %.5f from below" % (zo, tgt), 0 < tgt - zo < 0.01)
print("%d/%d passed" % (sum(res), len(res)))
