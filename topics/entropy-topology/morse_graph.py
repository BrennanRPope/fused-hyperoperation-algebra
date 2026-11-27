"""Conley-Morse graphs of the dial square on RP^1, computed on a grid (outer approximation).

Method of Arai, Kalies, Kokubu, Mischaikow, Oka and Pilarczyk (2009), in its simplest 1D form:
cut RP^1 into cells, send each cell to every cell its image touches, and take the strongly
connected pieces of that graph that contain a cycle. Those are the Morse sets; the graph between
them says which can flow to which. Attractors receive flow and pass none on; repellers the reverse.
"""
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components


def cell_graph(fh, ncell=1200, samples=12):
    edges_r, edges_c = [], []
    edges = np.linspace(0, np.pi, ncell + 1) + 0.01234567   # offset: keep fixed points off cell walls
    for i in range(ncell):
        th = np.linspace(edges[i], edges[i + 1], samples)
        U, V = fh(np.sin(th), np.cos(th))
        phi = np.mod(np.arctan2(U, V), np.pi)        # image angle on RP^1
        lift = np.unwrap(phi * 2) / 2                  # follow it continuously
        step = np.max(np.abs(np.diff(lift)))           # margin: largest gap between samples
        lo, hi = lift.min() - step, lift.max() + step
        off = 0.01234567
        js = np.arange(int(np.floor((lo - off) / np.pi * ncell)), int(np.floor((hi - off) / np.pi * ncell)) + 1)
        for j in np.mod(js, ncell):
            edges_r.append(i); edges_c.append(j)
    G = csr_matrix((np.ones(len(edges_r)), (edges_r, edges_c)), shape=(ncell, ncell))
    return G, edges


def morse_sets(G):
    n, lab = connected_components(G, directed=True, connection='strong')
    sets = []
    Gd = G.tocsr()
    for c in range(n):
        idx = np.nonzero(lab == c)[0]
        if len(idx) > 1 or Gd[idx[0], idx[0]] > 0:
            sets.append(idx)
    # reachability between Morse sets through the condensation
    from scipy.sparse.csgraph import breadth_first_order
    reach = np.zeros((len(sets), len(sets)), bool)
    owner = -np.ones(G.shape[0], int)
    for a, s in enumerate(sets):
        owner[s] = a
    for a, s in enumerate(sets):
        order = breadth_first_order(Gd, s[0], directed=True, return_predecessors=False)
        for b in set(owner[order]) - {-1}:
            if b != a:
                reach[a, b] = True
    return sets, reach


def classify(sets, reach):
    kinds = []
    for a in range(len(sets)):
        out_, in_ = reach[a].any(), reach[:, a].any()
        kinds.append('attractor' if not out_ else ('repeller' if not in_ else 'saddle'))
    return kinds


def set_centre(edges, s):
    """Centre of a Morse set as a point x of RP^1 (circular mean of its cells' angles)."""
    th = (edges[s] + edges[s + 1]) / 2
    a = np.angle(np.mean(np.exp(2j * th))) / 2
    return np.tan(a)


def merge_adjacent(sets, reach, ncell):
    """Join Morse sets whose cells touch (resolution artefacts next to a fast repeller)."""
    n = len(sets)
    parent = list(range(n))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a
    owner = -np.ones(ncell, int)
    for a, s in enumerate(sets):
        owner[s] = a
    for a, s in enumerate(sets):
        for c in s:
            for nb in ((c - 1) % ncell, (c + 1) % ncell):
                b = owner[nb]
                if b >= 0 and find(a) != find(b):
                    parent[find(a)] = find(b)
    roots = sorted(set(find(a) for a in range(n)))
    newsets = [np.concatenate([sets[a] for a in range(n) if find(a) == r]) for r in roots]
    m = len(roots)
    R = np.zeros((m, m), bool)
    for a in range(n):
        for b in range(n):
            if reach[a, b] and find(a) != find(b):
                R[roots.index(find(a)), roots.index(find(b))] = True
    return newsets, R


def morse_graph(fh, ncell=2000, samples=16):
    G, edges = cell_graph(fh, ncell, samples)
    S, R = morse_sets(G)
    S, R = merge_adjacent(S, R, ncell)
    return S, R, classify(S, R), edges
