"""Zeros 6551..6870 of zeta (window around Lehmer's pair near t = 7005) -> zeros_lehmer.json"""
import json, mpmath as mp
from multiprocessing import Pool
mp.mp.dps = 25
def g(n): return n, str(mp.zetazero(n).imag)
if __name__ == '__main__':
    with Pool(4) as P: out = P.map(g, range(6551, 6871), chunksize=8)
    json.dump({'n': [a for a, _ in out], 'gamma': [b for _, b in out]}, open('zeros_lehmer.json', 'w'))
    print('done', out[0], out[-1])
