"""The periodic law through the circle operation: shells are squares (SO(4)), periods follow the
Madelung diagonals, and the 20 exceptions only swap electrons between neighbouring diagonals.
Run: python3 checks.py"""
import json
ok = []
def check(name, cond, info=""):
    ok.append(bool(cond)); print(("PASS " if cond else "FAIL ") + name + ("  " + info if info else ""))
L = "spdfghik"
# 1. hydrogen shell n: subshells l = 0..n-1, each 2l+1 orbitals: a square of side n (gnomons)
for n in range(1, 30):
    assert sum(2*l+1 for l in range(n)) == n*n
    m = n-1
    assert n*n == 1 + (m + m + m*m)   # 1 + (n-1)∘₁(n-1)
check("shell n holds n² orbitals = 1 + (n−1)∘₁(n−1); its subshells are the square's L-shaped gnomons (n ≤ 29)", True)
# SO(4) = SU(2)xSU(2): two spins j = (n-1)/2; j⊗j decomposes into l = 0..2j
for n in range(1, 20):
    twoj = n-1
    assert [2*l+1 for l in range(0, twoj+1)] == [2*l+1 for l in range(n)]
check("SO(4): shell n = spin j ⊗ spin j with 2j = n−1; the Clebsch–Gordan series is exactly l = 0…n−1", True)

# 2. Madelung order (n+l, then n)
subs = sorted([(n, l) for n in range(1, 9) for l in range(n)], key=lambda x: (x[0]+x[1], x[0]))
order = [f"{n}{L[l]}" for n, l in subs]
check("Madelung order begins 1s 2s 2p 3s 3p 4s 3d 4p 5s 4d 5p 6s 4f 5d 6p 7s 5f 6d 7p",
      order[:19] == "1s 2s 2p 3s 3p 4s 3d 4p 5s 4d 5p 6s 4f 5d 6p 7s 5f 6d 7p".split())
def madelung(Z):
    conf = {}; left = Z
    for n, l in subs:
        if left == 0: break
        k = min(left, 2*(2*l+1)); conf[(n, l)] = k; left -= k
    return conf
# periods: a new period starts at each new s subshell
noble = []; tot = 0
for n, l in subs:
    tot += 2*(2*l+1)
    nxt = subs[subs.index((n, l))+1] if subs.index((n,l))+1 < len(subs) else None
    if nxt and nxt[1] == 0: noble.append(tot)
check("noble gases fall at Z = 2, 10, 18, 36, 54, 86, 118", noble[:7] == [2, 10, 18, 36, 54, 86, 118])
periods = [noble[0]] + [b-a for a, b in zip(noble, noble[1:7])]
check("period lengths 2, 8, 8, 18, 18, 32, 32 = doubled squares 2·1², 2·2², 2·2², 2·3², …", periods == [2, 8, 8, 18, 18, 32, 32], f"{periods}")
diag = [sum(2*(2*l+1) for n, l in subs if n+l == s) for s in range(1, 9)]
check("each Madelung diagonal n+l = s holds 2·⌈s/2⌉²: 2,2,8,8,18,18,32,32 (Janet's left-step periods)", diag == [2, 2, 8, 8, 18, 18, 32, 32], f"{diag}")

# 3. the experimental exceptions (Wikipedia 'Aufbau principle' table: 20 elements)
def parse(s):
    out = {}
    for tok in s.split():
        n = int(tok[0]); l = L.index(tok[1]); out[(n, l)] = int(tok[2:]) if len(tok) > 2 else 1
    return out
core = {"Ar": 18, "Kr": 36, "Xe": 54, "Rn": 86}
exc = {  # Z: (core, valence) — experimental ground states
 24: ("Ar", "3d5 4s1"), 29: ("Ar", "3d10 4s1"), 41: ("Kr", "4d4 5s1"), 42: ("Kr", "4d5 5s1"),
 44: ("Kr", "4d7 5s1"), 45: ("Kr", "4d8 5s1"), 46: ("Kr", "4d10"), 47: ("Kr", "4d10 5s1"),
 57: ("Xe", "5d1 6s2"), 58: ("Xe", "4f1 5d1 6s2"), 64: ("Xe", "4f7 5d1 6s2"),
 78: ("Xe", "4f14 5d9 6s1"), 79: ("Xe", "4f14 5d10 6s1"),
 89: ("Rn", "6d1 7s2"), 90: ("Rn", "6d2 7s2"), 91: ("Rn", "5f2 6d1 7s2"), 92: ("Rn", "5f3 6d1 7s2"),
 93: ("Rn", "5f4 6d1 7s2"), 96: ("Rn", "5f7 6d1 7s2"), 103: ("Rn", "5f14 7s2 7p1"),
}
moves = {}
for Z, (c, val) in exc.items():
    expt = madelung(core[c]); expt = {k: v for k, v in expt.items()}
    for k, v in parse(val).items(): expt[k] = expt.get(k, 0) + v
    assert sum(expt.values()) == Z, Z
    pred = madelung(Z)
    diff = {k: expt.get(k, 0) - pred.get(k, 0) for k in set(expt) | set(pred) if expt.get(k, 0) != pred.get(k, 0)}
    moves[Z] = diff
moved = sum(abs(v) for d in moves.values() for v in d.values())//2
check("all 20 exceptions differ from Madelung", all(moves[Z] for Z in exc))
maxmove = max(sum(v for v in d.values() if v > 0) for d in moves.values())
check("every exception moves only 1 or 2 electrons", maxmove <= 2, f"largest move {maxmove}")
spans = [max(n+l for (n, l) in d) - min(n+l for (n, l) in d) for d in moves.values()]
check("every exception moves electrons between subshells at most one Madelung diagonal apart (|Δ(n+l)| ≤ 1)", max(spans) <= 1, f"spans {sorted(set(spans))}")
check("no exception involves a p subshell except lawrencium's 6d→7p (same diagonal, n+l = 8)",
      all(all(l != 1 for (n, l) in d) for Z, d in moves.items() if Z != 103))
same_diag = [Z for Z, s in zip(moves, spans) if s == 0]
check("the exceptions on one diagonal (d↔f or d↔p swaps) are La, Ce, Gd, Ac, Th, Pa, U, Np, Cm, Lr", same_diag == [57, 58, 64, 89, 90, 91, 92, 93, 96, 103], f"{same_diag}")
# 4. the dial: q-integers are circle powers of 1 with k = q - 1
def wpow(k, a, n):
    x = 0.0
    for _ in range(n): x = x + a + k*x*a
    return x
for q in [0.5, 0.9, 1.0, 1.1, 2.0]:
    for n in range(1, 10):
        qint = n if q == 1 else (q**n - 1)/(q - 1)
        assert abs(wpow(q-1, 1.0, n) - qint) < 1e-9
check("q-integers [n]_q = n copies of 1 under the circle operation with k = q − 1", True)
json.dump({"order": order, "moves": {str(Z): {f"{n}{L[l]}": v for (n, l), v in d.items()} for Z, d in moves.items()}}, open("results.json", "w"), indent=1)
print(f"{sum(ok)}/{len(ok)} checks pass")
