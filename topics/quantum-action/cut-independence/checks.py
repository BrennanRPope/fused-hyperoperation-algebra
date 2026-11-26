"""Cut-independence from the swap rule's q weights  (DRAFT, Nov 26 2026)

Question: 3D quantum gravity (Turaev-Viro) gives the same answer however space is cut into
tetrahedra. Does that rule follow from our swap rule's face weight q (yx = q xy)?

Chain tested here, one link per section:
  Q1  the swap rule's plane: (x+y)^n expands with q-binomials (exact, sympy)
  Q2  the symmetry that respects the swap rule must carry the same q (exact, sympy)
  Q3  that symmetry is U_q(sl2) on every degree (exact, sympy)
  Q4  gluing two and three copies (the twisted product rule) gives recoupling numbers that
      match the q-6j symbol of quantum gravity (numeric, real q)
  Q5  6j orthogonality, generic q and roots of unity
  Q6  the pentagon (Biedenharn-Elliott) identity, generic q and roots of unity
  Q7  Pachner 2-3 move as a local identity (root of unity)
  Q8  Pachner 1-4 move as a local identity, and what fixes the vertex weight
  Q9  global check: the whole state sum on two different triangulations of the 3-sphere
  Q10 why a root of unity: quantum dimensions are the even terms of our dial's own sequence,
      which hits zero at step 2r; at q = 1 the vertex weight diverges (Ponzano-Regge)
  Q11 our dials (face weight q = b/a, a + b = k, a b = d = 1): after the ruler change
      d' = d - k^2/4, and its sign is the sign of the cosmological constant
Run: python3 checks.py   (writes checks_output.txt and results.json)
"""
import itertools, json, math, random, time
import numpy as np
import sympy as sp

OUT = []
RES = {}
def log(*a):
    s = " ".join(str(x) for x in a); print(s); OUT.append(s)
def record(key, ok, detail):
    RES[key] = {"pass": bool(ok), "detail": detail}
    log(f"[{'PASS' if ok else 'FAIL'}] {key}: {detail}")

T0 = time.time()
q = sp.symbols("q", positive=True)

# ---------------------------------------------------------------------------------------------
# Quantum-plane polynomials: dict {(a, b): coeff} for the normal-ordered monomial x^a y^b.
# The swap rule y x = q x y gives  (x^a y^b)(x^c y^d) = q^(b c) x^(a+c) y^(b+d).
def pmul(f, g, Q=q):
    h = {}
    for (a, b), cf in f.items():
        for (c, d), cg in g.items():
            k = (a + c, b + d); h[k] = h.get(k, 0) + cf * cg * Q**(b * c)
    return h
def padd(f, g, s=1):
    h = dict(f)
    for k, v in g.items(): h[k] = h.get(k, 0) + s * v
    return h
def psimp(f):
    out = {}
    for k, v in f.items():
        v = sp.simplify(v)
        if v != 0: out[k] = v
    return out
X = {(1, 0): 1}; Y = {(0, 1): 1}

log("=" * 100)
log("Q1  The swap rule's plane: (x+y)^n with y x = q x y")
def gauss_binom(n, k, Q):
    num = 1; den = 1
    for i in range(k):
        num *= (1 - Q**(n - i)); den *= (1 - Q**(i + 1))
    return sp.cancel(num / den)
ok1 = True
for n in range(1, 9):
    p = {(0, 0): 1}
    for _ in range(n): p = pmul(p, padd(X, Y))
    for kk in range(n + 1):
        c = sp.expand(p.get((kk, n - kk), 0))
        if sp.simplify(c - sp.expand(gauss_binom(n, kk, q))) != 0: ok1 = False
p3 = pmul(pmul(padd(X, Y), padd(X, Y)), padd(X, Y))
log("  (x+y)^3 =", " + ".join(f"({sp.factor(sp.expand(v))}) x^{a} y^{b}" for (a, b), v in sorted(p3.items())))
record("Q1 swap rule gives q-binomials", ok1,
       "coefficient of x^k y^(n-k) in (x+y)^n is the Gauss q-binomial for every n <= 8 (exact)")

# ---------------------------------------------------------------------------------------------
log("=" * 100)
log("Q2  A symmetry that respects the swap rule must carry the swap rule's own q")
p_ = sp.symbols("p", positive=True)   # the symmetry's parameter, left free
def act(op, f, Q, P):
    """Act on a quantum-plane polynomial (plane parameter Q) with E, F or K of a symmetry whose
    parameter is P, extended to products by the twisted product rule (coproduct)
      K(fg) = K f K g,  E(fg) = E(f) K(g) + f E(g),  F(fg) = F(f) g + K^-1(f) F(g),
    starting from K x = P x, K y = y / P, E x = 0, E y = x, F x = y, F y = 0."""
    out = {}
    for (a, b), c in f.items():
        word = "x" * a + "y" * b
        r = act_word(op, word, Q, P)
        for kk, v in r.items(): out[kk] = out.get(kk, 0) + c * v
    return out
from functools import lru_cache
def act_word(op, word, Q, P):
    if word == "":
        return {(0, 0): 1} if op in ("K", "Ki") else {}
    g, rest = word[0], word[1:]
    G = X if g == "x" else Y
    R = {(rest.count("x"), rest.count("y")): 1}
    if op == "K":
        s = P if g == "x" else 1 / P
        return {kk: s * v for kk, v in pmul(G, act_word("K", rest, Q, P), Q).items()}
    if op == "Ki":
        s = 1 / P if g == "x" else P
        return {kk: s * v for kk, v in pmul(G, act_word("Ki", rest, Q, P), Q).items()}
    if op == "E":   # E(g . rest) = E(g) K(rest) + g E(rest)
        t1 = pmul(X, act_word("K", rest, Q, P), Q) if g == "y" else {}
        t2 = pmul(G, act_word("E", rest, Q, P), Q)
        return padd(t1, t2)
    if op == "F":   # F(g . rest) = F(g) rest + K^-1(g) F(rest)
        t1 = pmul(Y, R, Q) if g == "x" else {}
        s = 1 / P if g == "x" else P
        t2 = {kk: s * v for kk, v in pmul(G, act_word("F", rest, Q, P), Q).items()}
        return padd(t1, t2)
rel = padd(pmul(Y, X), pmul(X, Y), -q)      # y x - q x y, which is zero in the plane
# Act on the *word* "yx" and "xy" (before using the relation), then compare.
Eyx = act_word("E", "yx", q, p_); Exy = act_word("E", "xy", q, p_)
Fyx = act_word("F", "yx", q, p_); Fxy = act_word("F", "xy", q, p_)
Erel = psimp(padd(Eyx, Exy, -q)); Frel = psimp(padd(Fyx, Fxy, -q))
log("  E(yx - q xy) =", {k: sp.factor(v) for k, v in Erel.items()})
log("  F(yx - q xy) =", {k: sp.factor(v) for k, v in Frel.items()})
sol = sp.solve([sp.Eq(v, 0) for v in list(Erel.values()) + list(Frel.values())], p_, dict=True)
ok2 = sol == [{p_: q}]
record("Q2 symmetry's q = swap rule's q", ok2,
       f"E and F keep y x - q x y = 0 only when the symmetry's parameter p solves {sol} (exact)")

# ---------------------------------------------------------------------------------------------
log("=" * 100)
log("Q3  That symmetry is U_q(sl2) on every degree n (spin n/2)")
def op_matrix(op, n, Q, P, numeric=False):
    """Matrix of op on degree-n polynomials, basis x^i y^(n-i), i = 0..n."""
    M = sp.zeros(n + 1, n + 1) if not numeric else np.zeros((n + 1, n + 1))
    for i in range(n + 1):
        r = act_word(op, "x" * i + "y" * (n - i), Q, P)
        for (a, b), v in r.items():
            M[a, i] = (sp.simplify(v) if not numeric else float(v))
    return M
ok3 = True; rows3 = []
for n in range(1, 7):
    E = op_matrix("E", n, q, q); F = op_matrix("F", n, q, q); K = op_matrix("K", n, q, q); Ki = op_matrix("Ki", n, q, q)
    c1 = sp.simplify(K * E - q**2 * E * K)
    c2 = sp.simplify(E * F - F * E - (K - Ki) / (q - 1 / q))
    c3 = sp.simplify(K * Ki - sp.eye(n + 1))
    good = c1 == sp.zeros(n + 1, n + 1) and c2 == sp.zeros(n + 1, n + 1) and c3 == sp.eye(n + 1) * 0
    ok3 &= good; rows3.append(good)
log("  E on x^a y^b (n=3):", [sp.factor(v) for v in op_matrix("E", 3, q, q)[1:, :-1].diagonal()])
record("Q3 U_q(sl2) relations", ok3,
       "K E = q^2 E K and E F - F E = (K - 1/K)/(q - 1/q) hold exactly on degrees 1..6 (spins 1/2..3)")

# ---------------------------------------------------------------------------------------------
# numeric q-numbers and the Racah-Wigner q-6j (labels are doubled spins A = 2a)
def qnum(n, qv):
    """[n] = (q^n - q^-n)/(q - q^-1); qv may be real (>0, !=1) or complex on the unit circle."""
    if isinstance(qv, complex) or (isinstance(qv, (float, int)) and qv == 1):
        if qv == 1: return float(n)
        th = np.angle(qv); return math.sin(n * th) / math.sin(th)
    return (qv**n - qv**-n) / (qv - 1 / qv)
def qfact(n, qv):
    out = 1.0
    for i in range(1, n + 1): out *= qnum(i, qv)
    return out
def adm(A, B, C, rr=None):
    if (A + B + C) % 2: return False
    if A > B + C or B > A + C or C > A + B: return False
    if rr is not None and A + B + C > 2 * (rr - 2): return False
    return True
def delta(A, B, C, qv):
    return math.sqrt(abs(qfact((A + B - C) // 2, qv) * qfact((A - B + C) // 2, qv) * qfact((-A + B + C) // 2, qv)
                     / qfact((A + B + C) // 2 + 1, qv)))
SIXJ = {}
TRUNC = True   # switch used once, for the negative control in Q7
def sixj(A1, A2, A3, A4, A5, A6, qv):
    """q-6j {a1 a2 a3; a4 a5 a6} (Racah-Wigner normalisation, tetrahedral symmetry), doubled spins.
    Triads: (1,2,3), (1,5,6), (4,2,6), (4,5,3)."""
    key = (A1, A2, A3, A4, A5, A6, qv)
    if key in SIXJ: return SIXJ[key]
    tri = [(A1, A2, A3), (A1, A5, A6), (A4, A2, A6), (A4, A5, A3)]
    # at q = e^(i pi/r) every face must also obey the truncation 2(j1+j2+j3) <= 2(r-2)
    rr = round(math.pi / np.angle(qv)) if (isinstance(qv, complex) and TRUNC) else None
    if not all(adm(*t, rr) for t in tri): SIXJ[key] = 0.0; return 0.0
    pref = 1.0
    for t in tri: pref *= delta(*t, qv)
    a = [sum(t) // 2 for t in tri]
    b = [(A1 + A2 + A4 + A5) // 2, (A2 + A3 + A5 + A6) // 2, (A3 + A1 + A6 + A4) // 2]
    s = 0.0
    for z in range(max(a), min(b) + 1):
        den = 1.0
        for ai in a: den *= qfact(z - ai, qv)
        for bi in b: den *= qfact(bi - z, qv)
        s += (-1)**z * qfact(z + 1, qv) / den
    SIXJ[key] = pref * s; return pref * s

# ---------------------------------------------------------------------------------------------
log("=" * 100)
log("Q4  Gluing copies with the twisted product rule gives quantum gravity's recoupling numbers")
QV = 1.3
def rep(N):
    """Spin (N/2): E, F, K numeric matrices from the plane (Q3) and the invariant norms."""
    E = op_matrix("E", N, QV, QV, True); F = op_matrix("F", N, QV, QV, True); K = op_matrix("K", N, QV, QV, True)
    # invariant inner product: E^dagger = F K  (compatible with the twisted product rule)
    nrm = np.ones(N + 1)
    for i in range(N):   # E v_i = e v_(i+1);  F v_(i+1) = f v_i
        e = E[i + 1, i]; f = F[i, i + 1]; w = K[i + 1, i + 1]
        nrm[i + 1] = nrm[i] * f * w / e
    return E, F, K, nrm
REPS = {N: rep(N) for N in range(0, 7)}
# check unitarity of each rep in its norms: G^-1 E^T G = F K
uni = max(np.abs(np.diag(1 / n_) @ E.T @ np.diag(n_) - F @ K).max() for (E, F, K, n_) in REPS.values())
def tensor(*Ns):
    """E, F, K on V_N1 (x) V_N2 (x) ... via the twisted product rule, plus Gram diagonal."""
    mats = [REPS[N] for N in Ns]
    def kron(lst):
        out = np.array([[1.0]])
        for m in lst: out = np.kron(out, m)
        return out
    Id = [np.eye(N + 1) for N in Ns]; Ks = [m[2] for m in mats]
    Ki = [np.linalg.inv(k) for k in Ks]
    E = sum(kron(Id[:i] + [mats[i][0]] + Ks[i + 1:]) for i in range(len(Ns)))
    F = sum(kron(Ki[:i] + [mats[i][1]] + Id[i + 1:]) for i in range(len(Ns)))
    K = kron(Ks)
    G = kron([np.diag(m[3]) for m in mats]).diagonal()
    return E, F, K, G
def null_in(Mat, basis, G):
    """Vectors v in span(basis) with Mat v = 0; returns G-orthonormal list."""
    B = np.array(basis).T
    u, s, vt = np.linalg.svd(Mat @ B)
    rank = int((s > 1e-9 * max(1, s.max() if len(s) else 1)).sum())
    ns = vt[rank:].T
    vs = [B @ ns[:, i] for i in range(ns.shape[1])]
    out = []
    for v in vs:
        for w in out: v = v - (w @ (G * v)) * w
        nv = math.sqrt(v @ (G * v))
        if nv > 1e-10: out.append(v / nv)
    return out
def weight_basis(K, val):
    return [np.eye(K.shape[0])[i] for i in range(K.shape[0]) if abs(K[i, i] - val) < 1e-9 * max(1, abs(val))]
def coupled_space(Ns2, J2):
    """Orthonormal basis (columns) of the spin-J2/2 block inside V_N1 (x) V_N2."""
    E, F, K, G = tensor(*Ns2)
    hw = null_in(E, weight_basis(K, QV**J2), G)
    assert len(hw) == 1
    vecs = [hw[0]]
    for _ in range(J2):
        v = F @ vecs[-1]; v = v / math.sqrt(v @ (G * v)); vecs.append(v)
    return np.array(vecs).T, G
def recoupling(N1, N2, N3, J):
    """Overlaps <((N1 N2)N12 N3) J | (N1 (N2 N3)N23) J> of normalised highest weight vectors."""
    E, F, K, G = tensor(N1, N2, N3)
    left = {}; right = {}
    for N12 in range(abs(N1 - N2), N1 + N2 + 1, 2):
        if not adm(N12, N3, J): continue
        C, _ = coupled_space((N1, N2), N12)      # column t has weight N12 - 2t; e_m has weight 2m - N3
        basis = [np.kron(C[:, t], np.eye(N3 + 1)[m]) for t in range(N12 + 1) for m in range(N3 + 1)
                 if (N12 - 2 * t) + (2 * m - N3) == J]
        hw = null_in(E, basis, G); assert len(hw) == 1; left[N12] = hw[0]
    for N23 in range(abs(N2 - N3), N2 + N3 + 1, 2):
        if not adm(N1, N23, J): continue
        C, _ = coupled_space((N2, N3), N23)
        basis = [np.kron(np.eye(N1 + 1)[m], C[:, t]) for m in range(N1 + 1) for t in range(N23 + 1)
                 if (2 * m - N1) + (N23 - 2 * t) == J]
        hw = null_in(E, basis, G); assert len(hw) == 1; right[N23] = hw[0]
    return {(a, b): left[a] @ (G * right[b]) for a in left for b in right}, E, G
worst = 0.0; count4 = 0; Q4PAIRS = []
for N1, N2, N3 in itertools.product(range(1, 4), repeat=3):
    for J in range(0, N1 + N2 + N3 + 1):
        if (N1 + N2 + N3 - J) % 2: continue
        M, _, _ = recoupling(N1, N2, N3, J)
        if not M: continue
        for (N12, N23), ov in M.items():
            pred = math.sqrt(qnum(N12 + 1, QV) * qnum(N23 + 1, QV)) * abs(sixj(N1, N2, N12, N3, J, N23, QV))
            worst = max(worst, abs(abs(ov) - pred)); count4 += 1; Q4PAIRS.append((abs(ov), pred))
log(f"  unitarity of every spin in its invariant norms: max |G^-1 E^T G - F K| = {uni:.1e}")
record("Q4 recoupling = q-6j", worst < 1e-9,
       f"{count4} recoupling numbers built only from the plane and the twisted product rule (q = {QV}) "
       f"match sqrt([2j12+1][2j23+1])|q-6j| to {worst:.1e}")

# ---------------------------------------------------------------------------------------------
log("=" * 100)
log("Q5  Orthogonality of the q-6j")
def labels(rr, maxN):
    top = maxN if rr is None else min(maxN, rr - 2)
    return range(0, top + 1)
def ortho_err(qv, rr, maxN):
    worst = 0.0; n = 0
    L = list(labels(rr, maxN))
    for A, B, C, D in itertools.product(L, repeat=4):
        XL = range(0, A + B + 1) if rr is None else L      # the sum over x is never cut short
        xs = [x for x in XL if adm(A, B, x, rr) and adm(C, D, x, rr)]
        fs = [f for f in L if adm(A, D, f, rr) and adm(C, B, f, rr)]
        for f, g in itertools.product(fs, repeat=2):
            s = sum(qnum(x + 1, qv) * qnum(f + 1, qv) * sixj(A, B, x, C, D, f, qv) * sixj(A, B, x, C, D, g, qv) for x in xs)
            worst = max(worst, abs(s - (1.0 if f == g else 0.0))); n += 1
    return worst, n
o_gen, n_gen = ortho_err(QV, None, 5)
o_roots = {r: ortho_err(complex(math.cos(math.pi / r), math.sin(math.pi / r)), r, 8) for r in (3, 4, 5, 6, 7)}
log(f"  generic q = {QV}: {n_gen} sums, worst error {o_gen:.1e}")
for r, (e, n) in o_roots.items(): log(f"  q = exp(i pi/{r}), labels truncated at 2j <= r-2: {n} sums, worst error {e:.1e}")
record("Q5 orthogonality", o_gen < 1e-9 and all(e < 1e-9 for e, _ in o_roots.values()),
       f"sum_x [2x+1][2f+1]{{a b x; c d f}}{{a b x; c d g}} = delta_fg for q = {QV} and q = e^(i pi/r), r = 3..7")

# ---------------------------------------------------------------------------------------------
log("=" * 100)
log("Q6  The pentagon (Biedenharn-Elliott) identity")
def be_err(qv, rr, maxN, samples=None, seed=1):
    rnd = random.Random(seed); L = list(labels(rr, maxN)); worst = 0.0; n = 0
    combos = itertools.product(L, repeat=9)
    if samples: combos = (tuple(rnd.choice(L) for _ in range(9)) for _ in range(samples))
    for (a, b, c, d, e, f, p_, q_, r_) in combos:
        # RHS factors must be admissible to be a meaningful test
        if not (adm(p_, q_, r_, rr) and adm(e, a, r_, rr) and adm(e, d, q_, rr) and adm(a, d, p_, rr)
                and adm(f, b, r_, rr) and adm(f, c, q_, rr) and adm(b, c, p_, rr)): continue
        S = (a + b + c + d + e + f + p_ + q_ + r_)
        lhs = 0.0
        for x in (range(0, a + b + 1) if rr is None else L):   # the sum over x is never cut short
            t = sixj(a, b, x, c, d, p_, qv) * sixj(c, d, x, e, f, q_, qv) * sixj(e, f, x, b, a, r_, qv)
            if t: lhs += (-1)**((S + x) // 2) * qnum(x + 1, qv) * t
        rhs = sixj(p_, q_, r_, e, a, d, qv) * sixj(p_, q_, r_, f, b, c, qv)
        worst = max(worst, abs(lhs - rhs)); n += 1
    return worst, n
b_gen = be_err(QV, None, 4, samples=60000)
b_roots = {r: be_err(complex(math.cos(math.pi / r), math.sin(math.pi / r)), r, 8, samples=60000, seed=r) for r in (4, 5, 6, 7)}
log(f"  generic q = {QV}: {b_gen[1]} admissible label sets, worst error {b_gen[0]:.1e}")
for r, (e, n) in b_roots.items(): log(f"  q = exp(i pi/{r}): {n} admissible label sets, worst error {e:.1e}")
TRUNC = False; SIXJ.clear()
ctrl6 = {r: be_err(complex(math.cos(math.pi / r), math.sin(math.pi / r)), r, 8, samples=60000, seed=r) for r in (4, 5, 6, 7)}
TRUNC = True; SIXJ.clear()
for r, (e, n) in ctrl6.items(): log(f"  control, q = exp(i pi/{r}) but faces allowed past 2(j1+j2+j3) <= 2(r-2): worst error {e:.2f}")
CONTROL = {"Q6_untruncated": {str(r): e for r, (e, n) in ctrl6.items()}}
record("Q6 pentagon identity", all(e > 0.1 for e, _ in ctrl6.values()) and b_gen[0] < 1e-9 and all(e < 1e-9 for e, _ in b_roots.values()) and b_gen[1] > 100,
       "Biedenharn-Elliott holds for real q and at roots of unity with the root-of-unity cut on faces, and fails without that cut")
# ---------------------------------------------------------------------------------------------
log("=" * 100)
log("Q7/Q8  Pachner moves as local identities (Turaev-Viro weights)")
def qroot(r): return complex(math.cos(math.pi / r), math.sin(math.pi / r))
def dimw(A, qv): return (-1)**A * qnum(A + 1, qv)          # edge weight (-1)^(2j) [2j+1]
def tet(e01, e02, e12, e23, e13, e03, qv):
    """Tetrahedron weight for vertices 0..3 with the edge labels named; faces 012, 013, 023, 123.
    The phase (-1)^(sum of the six spins) is the usual Ponzano-Regge / Turaev-Viro one; with it and the
    edge weight (-1)^(2j)[2j+1] the 2-3 move is exactly the pentagon identity of Q6."""
    v = sixj(e01, e02, e12, e23, e13, e03, qv)
    return 0.0 if v == 0 else (1j)**(e01 + e02 + e12 + e23 + e13 + e03) * v
def W2(r):
    qv = qroot(r); return sum(qnum(A + 1, qv)**2 for A in range(0, r - 1))
def move23(r, samples=4000, seed=7):
    """Two tets (1234),(2345) sharing face 234  ->  three tets around new edge 15."""
    qv = qroot(r); L = list(range(r - 1)); rnd = random.Random(seed); worst = 0.0; n = 0; nz = 0
    K9 = ["12", "13", "14", "23", "24", "34", "25", "35", "45"]
    pool = (dict(zip(K9, v)) for v in itertools.product(L, repeat=9)) if r <= 5 else \
           ({k: rnd.choice(L) for k in K9} for _ in range(samples))      # every labelling for r <= 5
    for e in pool:
        faces = [("12", "13", "23"), ("12", "14", "24"), ("13", "14", "34"), ("23", "24", "34"),
                 ("23", "25", "35"), ("24", "25", "45"), ("34", "35", "45")]
        if not all(adm(e[a], e[b], e[c], r) for a, b, c in faces): continue
        lhs = tet(e["12"], e["13"], e["23"], e["34"], e["24"], e["14"], qv) * \
              tet(e["23"], e["24"], e["34"], e["45"], e["35"], e["25"], qv)
        rhs = 0.0
        for x in L:   # edge 15
            t = tet(e["12"], e["13"], e["23"], e["35"], e["25"], x, qv) * \
                tet(e["12"], e["14"], e["24"], e["45"], e["25"], x, qv) * \
                tet(e["13"], e["14"], e["34"], e["45"], e["35"], x, qv)
            if t: rhs += dimw(x, qv) * t
        worst = max(worst, abs(lhs - rhs)); n += 1; nz += abs(lhs) > 1e-12
    return worst, n, nz
def move14(r, samples=3000, seed=11):
    """One tet (1234)  ->  four tets around a new vertex 5 (new edges 15, 25, 35, 45)."""
    qv = qroot(r); L = list(range(r - 1)); rnd = random.Random(seed); worst = 0.0; n = 0; ratio = []
    K6 = ["12", "13", "14", "23", "24", "34"]
    pool = (dict(zip(K6, v)) for v in itertools.product(L, repeat=6)) if r <= 5 else \
           ({k: rnd.choice(L) for k in K6} for _ in range(samples))      # every labelling for r <= 5
    for e in pool:
        faces = [("12", "13", "23"), ("12", "14", "24"), ("13", "14", "34"), ("23", "24", "34")]
        if not all(adm(e[a], e[b], e[c], r) for a, b, c in faces): continue
        lhs = tet(e["12"], e["13"], e["23"], e["34"], e["24"], e["14"], qv)
        s = 0.0
        for x1, x2, x3, x4 in itertools.product(L, repeat=4):
            t = tet(e["12"], e["13"], e["23"], x3, x2, x1, qv)
            if not t: continue
            t *= tet(e["12"], e["14"], e["24"], x4, x2, x1, qv)
            if not t: continue
            t *= tet(e["13"], e["14"], e["34"], x4, x3, x1, qv) * tet(e["23"], e["24"], e["34"], x4, x3, x2, qv)
            if t: s += dimw(x1, qv) * dimw(x2, qv) * dimw(x3, qv) * dimw(x4, qv) * t
        if abs(lhs) > 1e-9: ratio.append((s / lhs).real)
        worst = max(worst, abs(lhs - s / W2(r))); n += 1
    return worst, n, ratio
m23 = {r: move23(r, samples=40000) for r in (3, 4, 5, 6, 7)}
for r, (w, n, nz) in m23.items(): log(f"  2-3 move, r = {r}: {n} boundary labellings ({'all' if r <= 5 else 'sampled'}; {nz} nonzero), worst error {w:.1e}")
record("Q7 Pachner 2-3 move", all(w < 1e-9 for w, _, _ in m23.values()),
       "two tetrahedra = sum over the new edge of three tetrahedra, at every tested root of unity")
m14 = {r: move14(r, samples=3000) for r in (3, 4, 5, 6)}
for r, (w, n, ratio) in m14.items():
    log(f"  1-4 move, r = {r}: {n} labellings ({'all' if r <= 5 else 'sampled'}), worst error {w:.1e}; sum/tet = {np.mean(ratio):.6f} "
        f"(spread {np.ptp(ratio):.1e}) vs w^2 = sum [2j+1]^2 = {W2(r):.6f} = r/(2 sin^2(pi/r)) = {r/(2*math.sin(math.pi/r)**2):.6f}")
record("Q8 Pachner 1-4 move", all(w < 1e-9 for w, _, _ in m14.values()),
       "one tetrahedron = (1/w^2) sum over four new edges of four tetrahedra; the factor is forced to be "
       "w^2 = sum of [2j+1]^2 = r/(2 sin^2(pi/r))")

# ---------------------------------------------------------------------------------------------
log("=" * 100)
log("Q9  Global check: the whole state sum on four triangulations of the 3-sphere")
def state_sum(tets, r):
    qv = qroot(r); L = list(range(r - 1))
    edges = sorted({tuple(sorted(p)) for t in tets for p in itertools.combinations(t, 2)})
    verts = sorted({v for t in tets for v in t})
    faces = sorted({tuple(sorted(f)) for t in tets for f in itertools.combinations(t, 3)})
    eidx = {e: i for i, e in enumerate(edges)}
    face_edges = [[eidx[tuple(sorted(p))] for p in itertools.combinations(f, 2)] for f in faces]
    # order edges; check a face once its last edge is set
    last = {}
    for fi, fe in enumerate(face_edges): last.setdefault(max(fe), []).append(fi)
    lab = [0] * len(edges); total = [0.0]
    tet_edges = []
    for t in tets:
        a, b, c, d = sorted(t)
        tet_edges.append([eidx[(a, b)], eidx[(a, c)], eidx[(b, c)], eidx[(c, d)], eidx[(b, d)], eidx[(a, d)]])
    def rec(i):
        if i == len(edges):
            w = 1.0
            for A in lab: w *= dimw(A, qv)
            for te in tet_edges:
                w *= tet(*[lab[k] for k in te], qv)
                if w == 0: return
            total[0] += w; return
        for A in L:
            lab[i] = A
            if all(adm(*[lab[k] for k in face_edges[fi]], r) for fi in last.get(i, [])): rec(i + 1)
    rec(0)
    z = total[0] / W2(r)**len(verts)
    assert abs(complex(z).imag) < 1e-12
    return complex(z).real, len(verts), len(edges), len(tets)
two_tets = [(0, 1, 2, 3), (0, 1, 2, 3)]                      # two tetrahedra glued face to face
bd4 = [tuple(v for v in range(5) if v != i) for i in range(5)]  # boundary of the 4-simplex
def do23(tets, t1, t2):
    s = set(t1) & set(t2); a = (set(t1) - s).pop(); b = (set(t2) - s).pop(); s = sorted(s)
    new = [tuple(sorted([a, b] + [s[i], s[j]])) for i, j in ((0, 1), (0, 2), (1, 2))]
    out = list(tets); out.remove(t1); out.remove(t2); return out + new
def do14(tets, t, v):
    out = list(tets); out.remove(t)
    return out + [tuple(sorted(set(t) - {w} | {v})) for w in t]
bd4_14 = do14(bd4, (1, 2, 3, 4), 5)                 # new vertex 5 inside tet 1234: 8 tets
bd4_14_23 = do23(bd4_14, (2, 3, 4, 5), (0, 2, 3, 4))  # 2-3 move across face 234: new edge 0-5, 9 tets
assert (0, 5) not in {tuple(sorted(p)) for t in bd4_14 for p in itertools.combinations(t, 2)}
Z9 = {}
for r in (3, 4, 5, 6):
    Ts = (two_tets, bd4, bd4_14, bd4_14_23)
    t0 = time.time(); Z9[r] = [state_sum(T, r) for T in Ts]
    log(f"  r = {r}: " + "; ".join(f"V={v} E={e} T={t}: Z = {z:.12f}" for z, v, e, t in Z9[r])
        + f";  1/w^2 = {1/W2(r):.12f}   ({time.time()-t0:.0f} s)")
ok9 = all(max(abs(z - 1 / W2(r)) for z, _, _, _ in Z9[r]) < 1e-9 for r in Z9)
record("Q9 same answer on every triangulation", ok9,
       "Z(3-sphere) is the same, 1/w^2, for 2, 5, 8 and 9 tetrahedra (a 1-4 move, another 1-4 move, then a 2-3 move), r = 3, 4, 5, 6")

# ---------------------------------------------------------------------------------------------
log("=" * 100)
log("Q10 Why a root of unity: the dial sequence and the vertex weight")
# Our dial convention (polar-recast genus pass): 1 + kx + dx^2 = (1 + a x)(1 + b x), a + b = k, a b = d, and the
# face weight is q = b/a. With d = 1: a = q^(-1/2), b = q^(1/2), k = q^(1/2) + q^(-1/2).
# The dial's own sequence u(n+1) = k u(n) - d u(n-1), u(0) = 0, u(1) = 1, is u(n) = (b^n - a^n)/(b - a).
def dial_seq(qv, N):
    if isinstance(qv, complex): sq = complex(math.cos(np.angle(qv) / 2), math.sin(np.angle(qv) / 2))
    else: sq = math.sqrt(qv)
    k = (sq + 1 / sq); k = k.real if isinstance(k, complex) else k
    u = [0.0, 1.0]
    for n in range(1, N): u.append(k * u[n] - u[n - 1])
    return u, k
dim_err = 0.0
for qv in (QV, qroot(5), qroot(9)):
    u, k = dial_seq(qv, 60)
    for n in range(1, 25): dim_err = max(dim_err, abs(qnum(n, qv) - u[2 * n] / u[2]))
first_zero_u = {r: next(n for n in range(1, 10 * r) if abs(dial_seq(qroot(r), 10 * r)[0][n]) < 1e-9) for r in (3, 5, 8)}
first_zero_dim = {r: next(n for n in range(1, 10 * r) if abs(qnum(n, qroot(r))) < 1e-12) for r in (3, 5, 8)}
pos = all(qnum(n, qroot(r)) > 0 for r in (3, 5, 8, 20) for n in range(1, r))
w2 = {r: W2(r) for r in (3, 5, 10, 20, 50, 100, 200)}
log(f"  quantum dimension [n] = u(2n)/u(2), u the dial sequence with k = q^(1/2) + q^(-1/2), d = 1: worst error {dim_err:.1e}")
log(f"  first zero of the dial sequence u: n = {first_zero_u};  first zero of [n]: n = {first_zero_dim};  [n] > 0 below r: {pos}")
log("  w^2 = sum [2j+1]^2 grows like r^3/(2 pi^2): " + ", ".join(f"r={r}: {v:.1f} ({v/(r**3/(2*math.pi**2)):.4f})" for r, v in w2.items()))
record("Q10 the root of unity is what makes the sum finite",
       dim_err < 1e-9 and pos and all(first_zero_dim[r] == r and first_zero_u[r] == 2 * r for r in first_zero_u),
       "quantum dimensions are the even terms of our dial's own sequence; it first hits zero at step 2r, so [2j+1] "
       "first vanishes at 2j+1 = r and the labels stop; as q -> 1 (r -> infinity) the vertex weight w^2 diverges like r^3 "
       "(the Ponzano-Regge infinity)")

# ---------------------------------------------------------------------------------------------
log("=" * 100)
log("Q11 Our dials: face weight q = b/a with a + b = k, a b = d = 1, so k = q^(1/2) + q^(-1/2)")
def kd(qv):
    return dial_seq(qv, 2)[1]
def dprime(qv): return 1 - kd(qv)**2 / 4
tab = []
for name, qv, lam in [("q = e^(i pi/5)", qroot(5), "+"), ("q = e^(i pi/50)", qroot(50), "+"), ("q = 1", 1.0, "0"),
                      ("q = 1.02", 1.02, "-"), ("q = 1.3", 1.3, "-")]:
    dp = dprime(qv); tab.append((name, dp, lam))
    log(f"  {name:>16}: k = {kd(qv):.6f}, d' = d - k^2/4 = {dp:+.6f}   (3D quantum gravity with this q: Lambda {lam})")
sin_ok = all(abs(dprime(qroot(r)) - math.sin(math.pi / (2 * r))**2) < 1e-12 for r in range(3, 200))
sign_ok = all((dp > 1e-15 and lam == "+") or (abs(dp) < 1e-15 and lam == "0") or (dp < -1e-15 and lam == "-") for _, dp, lam in tab)
def comp(qv): return np.array([[kd(qv), -1.0], [1.0, 0.0]])     # step matrix of the dial, eigenvalues a and b
close = all(np.abs(np.linalg.matrix_power(comp(qroot(r)), 2 * r) + np.eye(2)).max() < 1e-8 for r in range(2, 40))
grow = np.abs(np.linalg.matrix_power(comp(1.3), 40)).max()
scal = {r: dprime(qroot(r)) * r**2 for r in (20, 160, 400, 4000)}
log(f"  d' = sin^2(pi/2r) exactly at q = e^(i pi/r): {sin_ok};  step matrix [[k,-1],[1,0]]^(2r) = -1 (r = 2..39): {close};"
    f"  at q = 1.3 it grows (entry {grow:.1e} after 40 steps)")
log("  r^2 d' -> pi^2/4, so d' falls like 1/r^2, the way Lambda does: " + ", ".join(f"r={r}: {v:.6f}" for r, v in scal.items())
    + f"  (pi^2/4 = {math.pi**2/4:.6f})")
record("Q11 sign of d - k^2/4 = sign of Lambda", sin_ok and sign_ok and close,
       "after the ruler change d' = d - k^2/4; roots of unity sit on the sphere side (d' = sin^2(pi/2r) > 0, Lambda > 0), "
       "q = 1 on the flat line (Ponzano-Regge, Lambda = 0), real q on the hyperbolic side (Lambda < 0); d' scales as "
       "1/r^2 like Lambda (the constant in front depends on the length unit, so only the scaling is claimed)")

log("=" * 100)
npass = sum(v["pass"] for v in RES.values() if "pass" in v)
log(f"{npass}/{len(RES)} checks pass   ({time.time() - T0:.1f} s)")
RES["_control"] = CONTROL
RES["_meta"] = {"pass": npass, "total": len([k for k in RES if not k.startswith("_")]), "seconds": round(time.time() - T0, 1),
                "Q9": {str(r): [z for z, *_ in v] for r, v in Z9.items()},
                "w2": {str(r): v for r, v in w2.items()},
                "m14_ratio": {str(r): float(np.mean(v[2])) for r, v in m14.items()}, "Q4pairs": Q4PAIRS}
open("checks_output.txt", "w").write("\n".join(OUT) + "\n")
json.dump(RES, open("results.json", "w"), indent=1, default=str)
