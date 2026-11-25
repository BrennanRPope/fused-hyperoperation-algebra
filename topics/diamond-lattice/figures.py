"""Pictures for the page: Julia sets (Fisher-zero sets) of the diamond map in the t-plane as Q slides,
and the parameter picture in the Q-plane for the chromatic marked point t_c = 1/(1-Q)."""
import json, base64, io, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from dynamics import *

HOT = np.array([0.93, 0.55, 0.30]); COLD = np.array([0.36, 0.62, 0.86]); AF = np.array([0.35, 0.72, 0.55])
INK = np.array([0.12, 0.12, 0.16])

def render(lab, step, out, px):
    img = np.zeros(lab.shape + (3,))
    shade = np.where(step <= 3, 1.0, np.where(step <= 6, 0.86, 0.74))
    for k, c in ((0, HOT), (1, COLD), (2, AF)):
        m = lab == k
        img[m] = c[None, :]*shade[m][:, None]
    img[lab == 3] = INK
    # boundary: neighbouring labels differ
    b = np.zeros(lab.shape, bool)
    b[1:, :] |= lab[1:, :] != lab[:-1, :]; b[:, 1:] |= lab[:, 1:] != lab[:, :-1]
    img[b] = INK
    fig = plt.figure(figsize=(px/100, px/100), dpi=100)
    ax = fig.add_axes([0, 0, 1, 1]); ax.axis('off'); ax.imshow(img, origin='lower', interpolation='nearest')
    buf = io.BytesIO(); fig.savefig(buf, format='png'); plt.close(fig)
    open(out, 'wb').write(buf.getvalue())
    return base64.b64encode(buf.getvalue()).decode()

B = beraha
panels = [('0', 0.0, 'Beraha n = 2: the walls merge'), ('0.6', 0.6, 'antiferro basin (green)'), ('1', 1.0, 'Beraha n = 3: percolation'),
          ('32/27', 32/27, 'fold: antiferro point dies'), ('1.41', 1.40980, 'real crossing'), ('1.64', 1.6388969, 'real crossing q1'),
          ('2', 2.0, 'Beraha n = 4: Ising'), ('2.618', B(5), 'Beraha n = 5'), ('3', 3.0, 'Beraha n = 6: three colours'),
          ('3.247', B(7), 'Beraha n = 7'), ('4', 4.0, 'Beraha limit'), ('8', 8.0, ''), ('16', 16.0, 't* = 1/3')]
data = {'julia': [], 'lim': 3.0}
for name, Q, note in panels:
    slow = Q < 0.01 or abs(Q - 32/27) < 1e-6
    it = 1500 if slow else 160
    tol = 3e-2 if slow else 1e-3
    xs, ys, lab, st = basin_grid(Q, xlim=(-3,3), ylim=(-3,3), n=340, iters=it, tol=tol)
    b64 = render(lab, st, f'figures/julia_Q{name.replace("/","_")}.png', 340)
    data['julia'].append({'name': name, 'Q': Q, 'note': note, 'png': b64,
                          'frac': {k: float(np.mean(lab == i)) for i, k in enumerate(['hot', 'cold', 'af', 'undecided'])}})

# Q-plane: where does the chromatic marked point t_c = 1/(1-Q) go?
n = 420
xs = np.linspace(-1.2, 4.6, n); ys = np.linspace(-2.9, 2.9, n)
Qg = xs[None, :] + 1j*ys[:, None]
Z = 1/(1 - Qg)
lab = np.full(Qg.shape, 3, np.int8); step = np.zeros(Qg.shape)
with np.errstate(all='ignore'):
    for k in range(200):
        Z = T(Z, Qg); Z = np.where(np.isfinite(Z), Z, (Qg - 2)/(Qg - 1))
        und = lab == 3
        lab[und & (np.abs(Z) < 1e-4)] = 0
        lab[und & (np.abs(Z - 1) < 1e-4)] = 1
        # third attracting fixed point (antiferro) : test by small move
        W = T(Z, Qg)
        lab[und & (lab == 3) & (np.abs(W - Z) < 1e-9) & (np.abs(Z) > 1e-3) & (np.abs(Z - 1) > 1e-3)] = 2
        step[und & (lab != 3)] = k + 1
data['qplane'] = {'png': render(lab, step, 'figures/qplane.png', 420), 'x0': -1.2, 'x1': 4.6, 'y0': -2.9, 'y1': 2.9}

# nu(Q) curve and the cubic Q(y)
Qs = np.concatenate([np.linspace(0.05, 1, 40), np.linspace(1, 20, 200)[1:]])
data['nu'] = [[float(q), float(nu(q))] for q in Qs]
data['beraha'] = [[n_, float(beraha(n_)), float(nu(beraha(n_)))] for n_ in range(3, 13)]
json.dump(data, open('figures/figdata.json', 'w'))
print('done', [ (d['name'], d['frac']) for d in data['julia']])
