"""Are the circle operation's 1,2,4,8 and the superstring's 3,4,6,10 one fact?
Cayley-Dickson algebras R(1), C(2), H(4), O(8), S(16) built exactly.
Identities checked symbolically (sympy) for dims 1-8 and by an exact integer
counterexample for 16. Writes results_div.json and checks_div_output.txt."""
import json, random, itertools
import sympy as sp

def conj(x):
    return [x[0]] + [-v for v in x[1:]]

def mul(x, y):
    n = len(x)
    if n == 1:
        return [x[0]*y[0]]
    h = n//2
    a, b, c, d = x[:h], x[h:], y[:h], y[h:]
    # (a,b)(c,d) = (ac - d* b, d a + b c*)
    left = [p - q for p, q in zip(mul(a, c), mul(conj(d), b))]
    right = [p + q for p, q in zip(mul(d, a), mul(b, conj(c)))]
    return left + right

def add(x, y): return [p+q for p, q in zip(x, y)]
def sub(x, y): return [p-q for p, q in zip(x, y)]
def scal(s, x): return [s*v for v in x]
def nsq(x): return sum(v*v for v in x)
def one(n): return [1] + [0]*(n-1)

def symvec(name, n): return list(sp.symbols(f'{name}0:{n}'))
def zero(expr_list): return all(sp.expand(e) == 0 for e in expr_list)

NAMES = {1: 'R', 2: 'C', 4: 'H (quaternions)', 8: 'O (octonions)', 16: 'S (sedenions)'}
random.seed(7)

def rint(n): return [random.randint(-3, 3) for _ in range(n)]

def props(n):
    if n <= 8:
        x, y, z = symvec('x', n), symvec('y', n), symvec('z', n)
        test = lambda f: zero(f(x, y, z))
    else:
        # exact integer search: identity holds only if every trial gives 0
        def test(f):
            for _ in range(300):
                if any(v != 0 for v in f(rint(n), rint(n), rint(n))):
                    return False
            return True
    res = {}
    res['norm multiplies |xy|=|x||y|'] = test(lambda x, y, z: [nsq(mul(x, y)) - nsq(x)*nsq(y)])
    res['commutative'] = test(lambda x, y, z: sub(mul(x, y), mul(y, x)))
    res['associative'] = test(lambda x, y, z: sub(mul(mul(x, y), z), mul(x, mul(y, z))))
    res['alternative'] = test(lambda x, y, z: sub(mul(mul(x, x), y), mul(x, mul(x, y))) + sub(mul(mul(y, x), x), mul(y, mul(x, x))))
    res['power-associative'] = test(lambda x, y, z: sub(mul(mul(x, x), x), mul(x, mul(x, x))))
    # circle operation a o_k b = a + b + k ab ; lift phi = 1 + k a
    k = sp.Symbol('k') if n <= 8 else 3
    circ = lambda a, b: add(add(a, b), scal(k, mul(a, b)))
    lift = lambda a: add(one(n), scal(k, a))
    res['lift is exact: 1+k(a o b) = (1+ka)(1+kb)'] = test(lambda x, y, z: sub(lift(circ(x, y)), mul(lift(x), lift(y))))
    res['circle op associative'] = test(lambda x, y, z: sub(circ(circ(x, y), z), circ(x, circ(y, z))))
    res['circle op polar split |1+k(a o b)| = |1+ka||1+kb|'] = test(lambda x, y, z: [nsq(lift(circ(x, y))) - nsq(lift(x))*nsq(lift(y))])
    # superstring 3-psi rule: psi=(a,b) in K^2, (psi psi^dagger)~ psi = 0
    def three_psi(a, b, z):
        c1 = add(scal(-nsq(b), a), mul(mul(a, conj(b)), b))
        c2 = sub(mul(mul(b, conj(a)), a), scal(nsq(a), b))
        return c1 + c2
    res['superstring 3-psi rule'] = test(three_psi)
    return res

out = {}
for n in [1, 2, 4, 8, 16]:
    out[n] = props(n)

# sedenion zero divisor (explicit): search basis pairs e_i+e_j, e_k+-e_l
def basis(n, i): v = [0]*n; v[i] = 1; return v
zd = None
for i, j, kk, l in itertools.product(range(1, 16), repeat=4):
    if i < j and kk < l:
        for s in (1, -1):
            a = add(basis(16, i), basis(16, j)); b = add(basis(16, kk), scal(s, basis(16, l)))
            if all(v == 0 for v in mul(a, b)):
                zd = (i, j, kk, l, s); break
    if zd: break

# 3-psi residual size on random unit-ish spinors, for the bar chart
def residual(n, trials=2000):
    random.seed(11); worst = 0.0; tot = 0.0
    for _ in range(trials):
        a = [random.gauss(0, 1) for _ in range(n)]; b = [random.gauss(0, 1) for _ in range(n)]
        s = (nsq(a) + nsq(b))**0.5; a = scal(1/s, a); b = scal(1/s, b)
        c1 = add(scal(-nsq(b), a), mul(mul(a, conj(b)), b))
        c2 = sub(mul(mul(b, conj(a)), a), scal(nsq(a), b))
        r = nsq(c1 + c2)**0.5; tot += r; worst = max(worst, r)
    return tot/trials, worst
resid = {n: residual(n) for n in [1, 2, 4, 8, 16]}

# spacetime = 2x2 hermitian over K: X = [[t+z, xbar],[x, t-z]], det = (t+z)(t-z) - |x|^2
# the diagonal (t+z)(t-z) is the split-complex (d<0, relativity plane) norm
st = {}
for n in [1, 2, 4, 8]:
    t, zz = sp.symbols('t zz'); xv = symvec('w', n)
    det = (t+zz)*(t-zz) - nsq(xv)
    target = t**2 - zz**2 - sum(v**2 for v in xv)
    st[n] = {'spacetime dim': n+2, 'det is Minkowski norm': sp.expand(det - target) == 0,
             'signature': f'(1,{n+1})', 'Lorentz algebra dim (n+2)(n+1)/2': (n+2)*(n+1)//2}

lines = []
ok = True
expect = {
 1: dict(nm=1, co=1, asc=1, alt=1), 2: dict(nm=1, co=1, asc=1, alt=1),
 4: dict(nm=1, co=0, asc=1, alt=1), 8: dict(nm=1, co=0, asc=0, alt=1), 16: dict(nm=0, co=0, asc=0, alt=0)}
for n in out:
    r = out[n]; e = expect[n]
    good = (r['norm multiplies |xy|=|x||y|'] == bool(e['nm']) and r['commutative'] == bool(e['co'])
            and r['associative'] == bool(e['asc']) and r['alternative'] == bool(e['alt'])
            and r['lift is exact: 1+k(a o b) = (1+ka)(1+kb)']
            and r['circle op associative'] == r['associative']
            and r['circle op polar split |1+k(a o b)| = |1+ka||1+kb|'] == r['norm multiplies |xy|=|x||y|']
            and r['superstring 3-psi rule'] == r['alternative'])
    ok &= good
    lines.append(f"dim {n:2d} {NAMES[n]:18s} " + ", ".join(f"{k}: {'yes' if v else 'no'}" for k, v in r.items()) + f"  [{'PASS' if good else 'FAIL'}]")
lines.append(f"sedenion zero divisor: (e{zd[0]}+e{zd[1]})(e{zd[2]}{'+' if zd[4]>0 else '-'}e{zd[3]}) = 0")
for n in resid:
    lines.append(f"3-psi residual dim {n:2d}: mean {resid[n][0]:.3e}, worst {resid[n][1]:.3e}")
for n in st:
    lines.append(f"spacetime over dim {n}: {st[n]}")
    ok &= st[n]['det is Minkowski norm']
lines.append('ALL PASS' if ok and zd else 'SOMETHING FAILED')
txt = "\n".join(lines); print(txt)
open('checks_div_output.txt', 'w').write(txt + "\n")
json.dump({'props': {str(n): out[n] for n in out}, 'zero_divisor': zd,
           'residual': {str(n): resid[n] for n in resid}, 'spacetime': {str(n): st[n] for n in st}},
          open('results_div.json', 'w'), indent=1)
