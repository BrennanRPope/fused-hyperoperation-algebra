"""Round 4: a three-body style topological word for the rewind.
Rewind (t decreasing) the Lehmer window of zeros (6551..6870, z = 2*gamma) with the zero ODE,
fusing each pair at its merger into (s,p) = (z1+z2, z1*z2) so the flow runs straight through
(round-3 fused coordinates). After a merger a pair is complex: s/2 +- i*sqrt(p - s^2/4).
A letter is a 'syzygy': the real part s/2 of a complex pair passes a real zero or another pair's
centre (the zeros become vertically aligned, the analogue of the three-body collinearity).
Background: zeros outside the window frozen, plus mirrors and tail, as in round 3. -> braid.json"""
import json, numpy as np
G = np.array([float(x) for x in json.load(open('../elliptic-object/zeros.json'))['gamma']])
L = json.load(open('../round3/zeros_lehmer.json')); GL = np.array([float(x) for x in L['gamma']]); NL = L['n']
z_all = 2*np.concatenate([G, GL]); zw0 = 2*GL
outz = 2*G   # zeros below the window (frozen); above the window: tail density
zmax = zw0[-1]
wq = np.linspace(zmax, zmax*400, 200001); rho = np.log(wq/(4*np.pi))/(4*np.pi)
def bgF(x):  # background part of (log H)' at real x
    return np.sum(1/(x-outz)) + np.sum(1/(x+z_all)) + np.trapezoid((1/(x-wq)+1/(x+wq))*rho, wq)
# smooth background: tabulate on a grid across the window
xg = np.linspace(zw0[0]-2, zw0[-1]+2, 4001); bg_tab = np.array([bgF(x) for x in xg])
bg = lambda x: np.interp(x, xg, bg_tab)
# state: list of objects; ('r', z) real zero or ('f', s, p) fused pair
def roots(st):
    rs = []
    for o in st:
        if o[0] == 'r': rs.append(complex(o[1]))
        else:
            s, p = o[1], o[2]; D = s*s-4*p; q = np.sqrt(complex(D))
            rs += [(s+q)/2, (s-q)/2]
    return np.array(rs)
def deriv(st):
    R = roots(st); out = []; i = 0
    for o in st:
        if o[0] == 'r':
            z = R[i]; d = z-np.delete(R, i)
            out.append((2*np.sum(1/d)+2*bg(z.real)).real); i += 1
        else:
            z1, z2 = R[i], R[i+1]; others = np.delete(R, [i, i+1])
            F1 = 2*np.sum(1/(z1-others))+2*bg(z1.real); F2 = 2*np.sum(1/(z2-others))+2*bg(z2.real)
            # background enters as a real smooth field; for complex roots use its value at the real part
            ds = F1+F2; dp = -2+z2*F1+z1*F2
            out.append((ds.real, dp.real)); i += 2
    return out
def step(st, dt):
    def add(st, k, h):
        new = []
        for o, d in zip(st, k):
            if o[0] == 'r': new.append(('r', o[1]+h*d))
            else: new.append(('f', o[1]+h*d[0], o[2]+h*d[1]))
        return new
    k1 = deriv(st); k2 = deriv(add(st, k1, dt/2)); k3 = deriv(add(st, k2, dt/2)); k4 = deriv(add(st, k3, dt))
    new = []
    for o, a, b, c, d in zip(st, k1, k2, k3, k4):
        if o[0] == 'r': new.append(('r', o[1]+dt/6*(a+2*b+2*c+d)))
        else: new.append(('f', o[1]+dt/6*(a[0]+2*b[0]+2*c[0]+d[0]), o[2]+dt/6*(a[1]+2*b[1]+2*c[1]+d[1])))
    return new
def centres(st):
    return [o[1] if o[0] == 'r' else o[1]/2 for o in st]
st = [('r', z) for z in zw0]
lab = list(range(len(st)))         # label of each object = index of its first zero in the window
t, T_END = 0.0, -0.06
events = []; mergers = []; inner = (40, len(zw0)-40)
nstep = 0
while t > T_END:
    R = roots(st); Rr = np.sort_complex(R)
    dmin = min(abs(R[i]-R[j]) for i in range(len(R)) for j in range(i+1, min(i+4, len(R))))
    # fuse adjacent real zeros that are about to merge
    c = centres(st); fused = False
    for i in range(len(st)-1):
        if st[i][0] == 'r' and st[i+1][0] == 'r' and st[i+1][1]-st[i][1] < 0.02:
            a, b = st[i][1], st[i+1][1]
            st[i:i+2] = [('f', a+b, a*b)]; mergers.append(dict(t=t, at=float((a+b)/4), label=lab[i]))
            lab[i:i+2] = [lab[i]]; fused = True; break
    if fused: continue
    # smallest distance between neighbouring roots sets the step
    order = np.argsort([z.real for z in R]); Rs = R[order]
    dmin = np.min(np.abs(np.diff(Rs)))
    for o in st:
        if o[0] == 'f':
            D = o[1]**2-4*o[2]
            if D < 0: pass
            else: dmin = min(dmin, np.sqrt(D))
    dt = -min(1e-4, 0.01*max(dmin, 0.02)**2)
    c0 = centres(st)
    st = step(st, dt); t += dt; nstep += 1
    c1 = centres(st)
    # syzygy: a complex pair's centre passes another object's centre (order of centres changes)
    o0 = np.argsort(c0); o1 = np.argsort(c1)
    if not np.array_equal(o0, o1):
        for i in range(len(st)):
            for j in range(i+1, len(st)):
                if (c0[i]-c0[j])*(c1[i]-c1[j]) < 0:
                    ci = st[i][0] == 'f' and st[i][1]**2-4*st[i][2] < 0
                    cj = st[j][0] == 'f' and st[j][1]**2-4*st[j][2] < 0
                    if (ci or cj) and inner[0] <= lab[i] < inner[1]:
                        events.append(dict(t=t, x=float(c1[i]/2), a=lab[i], b=lab[j], kind=('pair-pair' if ci and cj else 'pair-zero')))
    if nstep % 500 == 0:
        print(round(t, 5), len(mergers), len(events), flush=True)
nz = sum(1 for o in st if o[0] == 'f' and o[1]**2-4*o[2] < 0)
res = dict(T_END=T_END, mergers=[m for m in mergers if inner[0] <= m['label'] < inner[1]], events=events,
           complex_pairs_end=nz, n_inner=inner[1]-inner[0], window_n0=NL[0],
           heights=[float(o[2]-o[1]**2/4)**0.5 if (o[0] == 'f' and o[1]**2-4*o[2] < 0) else 0.0 for o in st],
           centres=[float(c) for c in centres(st)])
json.dump(res, open('braid.json', 'w'))
print('done mergers', len(res['mergers']), 'syzygies', len(events), 'complex pairs', nz)
