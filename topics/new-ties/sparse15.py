"""Round 15 (DRAFT): a second search for Conjecture 15.3, in the other direction from kappa15.py: few N/K letters but
long runs of l and m between them (the shape of the words in a bad pair, whose l and m runs grow with a and b).

Every word is x1 u1 x2 u2 ... xk with x_i in {N, K} and each u_i = l^i m^j (ml = lm as maps, and after m only m or N may
follow, so j = 0 before a K). Leading letter x1 = N is enough (psi(Kw) = psi(Nw)); trailing l, m only multiply psi by a
unit, so they are left off. psi = num/19^e is computed in exact integers, the fingerprint is num with 2, 5 and 19 removed.
A clash is two words with the same fingerprint but different canonical forms (kappa15.canon).
Usage: python3 sparse15.py   Writes sparse15.json.
"""
import json, itertools, time
from kappa15 import kappa_int, canon
def runs(R):
    for i in range(R + 1):
        for j in range(R + 1):
            yield i, j
def psi_of(word_runs, letters):
    # right recursion: l: (num*200, e+1); m: (num*20, e); N: (num*19 + 19^e, e); K: (num*10 + 19^e, e); start N = (1, 0)
    num, e = 1, 0
    for (i, j), x in zip(word_runs, letters[1:]):
        num *= 200 ** i; e += i; num *= 20 ** j
        z = 19 ** e
        num = num * 19 + z if x == 'N' else num * 10 + z
    return num, e
def word_str(word_runs, letters):
    s = letters[0]
    for (i, j), x in zip(word_runs, letters[1:]): s += 'l' * i + 'm' * j + x
    return s
PLAN = [(1, 0), (2, 300), (3, 35), (4, 9)]   # (number of N/K letters, longest l run and longest m run allowed)
if __name__ == "__main__":
    t0 = time.time(); first = {}; clashes = []; counted = []; units = []
    for k, R in PLAN:
        cnt = 0
        for letters_tail in itertools.product('NK', repeat=k - 1):
            letters = ('N',) + letters_tail
            for wr in itertools.product(list(runs(R)), repeat=k - 1):
                if any(j > 0 and x == 'K' for (i, j), x in zip(wr, letters[1:])): continue
                num, e = psi_of(wr, letters); kap = kappa_int(num); w = word_str(wr, letters); c = canon(w); cnt += 1
                if kap == 1: units.append(w)
                if kap in first:
                    if first[kap] != c: clashes.append((first[kap], c))
                else: first[kap] = c
        counted.append(dict(NK_letters=k, longest_run=R, words=cnt))
        print(k, R, cnt, "words", len(first), "fingerprints", len(clashes), "clashes", round(time.time() - t0, 1), "s", flush=True)
    out = dict(plan=counted, words=sum(c['words'] for c in counted), distinct_fingerprints=len(first),
               clashes_after_canonical_form=len(clashes), clash_examples=clashes[:20],
               unit_psi_words=sorted(set(units), key=len)[:20], unit_psi_count=len(units))
    json.dump(out, open("sparse15.json", "w"), indent=1)
    print(out['words'], "words;", len(first), "fingerprints;", len(clashes), "clashes; unit-psi words:", out['unit_psi_words'])
