"""Pass 1, item 5: the kernel ladder exp(-|u|^p) on the roundness dial, and where Riemann's kernel sits.
G_p(x) = int_0^inf exp(-u^p) cos(xu) du.  Our N3 reading: a kink of |u|^p at the centre (p not even) pushes with size
Gamma(p+1) sin(pi p/2) x^(-p-1); with no rim the smooth part decays faster than any power, so the push wins and only
finitely many zeros are real.  At even p there is no kink and no push.  Writes ladder.json."""
import json, mpmath as mp
mp.mp.dps = 40
def G(p, x):
    p = mp.mpf(p); x = mp.mpf(x)
    U = mp.mpf(110)**(1/p)
    n = max(4, int(x*U/mp.pi) + 2)
    pts = mp.linspace(0, U, n+1)
    return mp.quad(lambda u: mp.exp(-u**p)*mp.cos(x*u), pts)
out = {}
for p in (2.5, 3, 3.5, 4, 5, 6):
    xs = [mp.mpf(i)/5 for i in range(0, 201)]   # x in [0, 40]
    vals = [G(p, x) for x in xs]
    ch = [float(xs[i]) for i in range(len(xs)-1) if vals[i]*vals[i+1] < 0]
    push = [float(mp.gamma(p+1)*mp.sin(mp.pi*p/2)*x**(-p-1)) if x > 0 else None for x in xs]
    out[str(p)] = dict(sign_changes=ch, x=[float(x) for x in xs],
                       val=[float(v) for v in vals], push=push,
                       ratio_at_40=float(vals[-1]/(mp.gamma(p+1)*mp.sin(mp.pi*p/2)*xs[-1]**(-p-1))) if p not in (4, 6) else None)
    print(p, len(ch), [round(c, 1) for c in ch][:14], out[str(p)]['ratio_at_40'], flush=True)
json.dump(out, open('ladder.json', 'w'))
