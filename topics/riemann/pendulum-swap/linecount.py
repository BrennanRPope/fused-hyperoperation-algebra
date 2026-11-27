"""Fine angle scan: are all zeros in the box on the line? Compare (a) zeros on the line = sign changes of the real
function F(iy), with (b) all zeros in the box = winding of F around the box edge. Every 1 degree. DRAFT."""
import numpy as np, sys, json
sys.path.insert(0, "../rewind-dial"); sys.path.insert(0, "../knauf-field")
from knauf import chain
from oddpath import setup_odd
from oddrewind_dual import setup_frame
from frames import to_bonds
def Fv(S, z, th):
    X = np.cos(th)*S["a"]*S["M"] + np.sin(th)*S["H"]; out = np.empty(z.size, complex)
    for i in range(0, z.size, 2000): out[i:i+2000] = (S["w"][None, :]*np.exp(np.outer(z[i:i+2000], X))).sum(1)
    return out
def counts(S, th, y0=0.3, y1=36.0, x=1.0):
    y = np.arange(y0, y1, 0.002); f = Fv(S, 1j*y, th).real; online = int(np.sum(np.sign(f[1:]) != np.sign(f[:-1])))
    e = np.concatenate([x + 1j*np.linspace(y0, y1, 9000), np.linspace(x, -x, 600) + 1j*y1,
                        -x + 1j*np.linspace(y1, y0, 9000), np.linspace(-x, x, 600) + 1j*y0])
    ph = np.unwrap(np.angle(Fv(S, e, th))); total = int(np.rint((ph[-1] - ph[0])/(2*np.pi)))
    return online, total
if __name__ == "__main__":
    out = {}
    for n in (8, 10, 12):
        frames = {"spin frame": setup_odd(n), "dual frame": setup_frame(to_bonds(np.log(chain(n).astype(float)), n), n)}
        for name, S in frames.items():
            bad = []
            for deg in range(0, 91):
                on, tot = counts(S, np.radians(deg))
                if on != tot: bad.append((deg, on, tot))
            out["%s n=%d" % (name, n)] = bad
            print(name, "n=%d: angles where some zero is off the line:" % n, bad if bad else "none (0 to 90 deg, every degree)", flush=True)
    json.dump(out, open("linecount.json", "w"), indent=1)
