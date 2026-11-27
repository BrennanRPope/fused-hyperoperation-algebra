"""Remaining angles of the n = 12 odd-half rewind, on a 0.02 grid (Newton refines each zero). DRAFT."""
import numpy as np, json, sys, movie
from oddpath import setup_odd
degs = [int(x) for x in sys.argv[1:]]; S = setup_odd(12); out = {}
for deg in degs:
    z = movie.zeros_box(S, np.radians(deg), h=0.02); z = z[(np.abs(z.real) <= 1) & (z.imag > 0.3) & (z.imag <= 36)]
    z = np.unique(np.round(z, 8)); out[deg] = [[float(v.real), float(v.imag)] for v in z]
    print("odd path n=12 %2d deg: %d zeros, max |Re z| = %.1e" % (deg, z.size, np.abs(z.real).max() if z.size else -1), flush=True)
json.dump(out, open("oddpath12_%s.json" % "_".join(sys.argv[1:]), "w"))
