"""Checks for pendulum-swap. Run: python3 checks.py. DRAFT."""
import numpy as np, json, sys, glob, mpmath as mp
sys.path.insert(0, "../knauf-field"); from knauf import chain
from frames import to_bonds, parts, odd_share, walsh, doubled
mp.mp.dps = 30; ok = []
def check(name, cond): ok.append(bool(cond)); print(("PASS " if cond else "FAIL ") + name)
# --- the swap on spins ---
for n in (8, 10):
    x = np.arange(2**n); y = x ^ (x >> 1); check("Gray (bond) map is one-to-one, n=%d" % n, len(set(y.tolist())) == 2**n)
    H = np.log(chain(n).astype(float)); Hb = to_bonds(H, n)
    check("partition function unchanged by the swap of frames, n=%d" % n, np.isclose(np.exp(-1.5*H).sum(), np.exp(-1.5*Hb).sum()))
K = 0.7; Ks = -0.5*np.log(np.tanh(K)); check("Kramers-Wannier dual is an involution (e^-2K* = tanh K and back)", np.isclose(-0.5*np.log(np.tanh(Ks)), K))
sz = lambda n: np.array([bin(k).count("1") for k in range(2**n)])
H = np.log(chain(12).astype(float)); He, _ = parts(H); J = walsh(He); s = sz(12)
check("even half is >= 99% pair links, n=12", (J[s == 2]**2).sum()/(J[1:]**2).sum() > 0.99)
bf = json.load(open("bestflip.json"))
check("staggered flip (dual frame) is the best of all flip choices, n=8..16", all(bf[n]["best"] == bf[n]["alternating"] for n in bf))
check("dual frame odd share in [0.80, 0.82], n=8..16", all(0.80 <= bf[n]["alternating"] <= 0.82 for n in bf))
check("flip-all odd share falls with length", all(bf[str(a)]["all"] > bf[str(b)]["all"] for a, b in ((8, 10), (10, 12), (12, 14), (14, 16))))
fr = json.load(open("frames.json"))
check("dual frame passes the field condition at beta = 3/2, n=8,10,12", all(fr[n][1]["cohn_minroot"] > 1 for n in fr))
check("doubled dual chain has every field zero on the circle, n=8,10,12", all(fr[n][1]["circle_dev"] < 1e-9 for n in fr))
check("dual frame pushes off the line more slowly than the spin frame (first order), n=8,10,12", all(fr[n][1]["speed_per_height"] < fr[n][0]["speed_per_height"] for n in fr))
lc = json.load(open("linecount_8_10.json"))
check("odd-half rewind, spin frame: no zero leaves the line, n=8,10, every degree", lc["spin frame n=8"] == [] and lc["spin frame n=10"] == [])
check("odd-half rewind, dual frame n=8: no zero leaves", lc["dual frame n=8"] == [])
check("odd-half rewind, dual frame n=10: zeros leave only between 5 and 23 degrees", [r[0] for r in lc["dual frame n=10"]] == list(range(5, 24)))
lcm = json.load(open("linecount.json"))
check("odd-half rewind, spin frame n=12: no zero leaves the line, every degree", lcm["spin frame n=12"] == [])
check("odd-half rewind, dual frame n=12: zeros leave only between 4 and 27 degrees", [r[0] for r in lcm["dual frame n=12"]] == list(range(4, 28)))
rows = [json.loads(l) for f in glob.glob("lc12_*.jsonl") for l in open(f)]
if rows:
    for fr_ in ("spin", "dual"):
        rr = sorted(json.loads(l) for f in glob.glob("lc12_%s*.jsonl" % fr_) for l in open(f))
        print("     n=12 %s frame: %d angles done, off-line at:" % (fr_, len(rr)), [r for r in rr if r[1] != r[2]])
# --- the area finding ---
g = lambda x: (-mp.mpf(1)/12 + x**2/240 - x**4/6048) if abs(x) < mp.mpf('1e-4') else 1/(4*mp.sinh(x/2)**2) - 1/x**2
check("area of the gap = -1/2 = zeta(0)", abs(mp.quad(g, [0, 1, 10, mp.inf]) + mp.mpf(1)/2) < 1e-20 and mp.zeta(0) == -mp.mpf(1)/2)
s = mp.mpc(1.5, 5); lhs = mp.quad(lambda x: x**(s-1)*g(x), [0, 1, 10, 40, mp.inf])
check("Mellin of the gap = Gamma(s) zeta(s-1) on the line Re s = 3/2 (s = 3/2 + 5i)", abs(lhs - mp.gamma(s)*mp.zeta(s-1)) < 1e-15)
x = mp.mpf(1.3); S = mp.nsum(lambda k: 2*(x**2 - 4*mp.pi**2*k**2)/(x**2 + 4*mp.pi**2*k**2)**2, [1, mp.inf])
check("Poisson swap: gap = sum over frequencies 2 pi k, k != 0", abs(S - g(x)) < 1e-20)
check("gap curve is negative everywhere (one sign)", all(g(mp.mpf(10)**(e/10.)) < 0 for e in range(-30, 30)))
X = mp.mpf(10)**6; check("zero-point: quantum minus classical energy -> 1/2 as it cools (x = 1e6)", abs(mp.coth(X/2)/2 - 1/X - mp.mpf(1)/2) < 2e-6)
t1 = mp.im(mp.zetazero(1)); check("area kernel's Fourier transform vanishes at the first zeta zero", abs(mp.gamma(mp.mpc(1.5, t1))*mp.zeta(mp.mpc(0.5, t1))) < 1e-20)
xi = lambda s: s*(s-1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
check("xi(0) = xi(1) = 1/2: the area's 1/2 is a value at the strip's edges", abs(xi(mp.mpf('1e-15')) - 0.5) < 1e-10 and abs(xi(1 + mp.mpf('1e-15')) - 0.5) < 1e-10)
Phi = lambda u: mp.nsum(lambda n: (2*mp.pi**2*n**4*mp.e**(4.5*u) - 3*mp.pi*n**2*mp.e**(2.5*u))*mp.e**(-mp.pi*n**2*mp.e**(2*u)), [1, mp.inf])
check("Riemann's kernel is even (the theta / pendulum swap), u = 0.4", abs(Phi(0.4) - Phi(-0.4)) < 1e-15)
uk = lambda u: mp.e**(1.5*u)*g(mp.e**u)
check("the area kernel is NOT even (no swap built in), u = 1", abs(uk(1) - uk(-1)) > 0.1*abs(uk(1)))
print("%d/%d passed" % (sum(ok), len(ok)))
