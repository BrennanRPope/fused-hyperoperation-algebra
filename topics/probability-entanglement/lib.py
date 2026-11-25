"""Shared helpers for the probability / time / entanglement checks (DRAFT, Nov 25 2026)."""
import numpy as np

def circ(a, b, k):            # a o_k b
    return a + b + k * a * b

def haar_state(rng, dim):
    v = rng.normal(size=dim) + 1j * rng.normal(size=dim)
    return v / np.linalg.norm(v)

def rho_of(psi):
    return np.outer(psi, psi.conj())

def ptrace(rho, dims, keep):
    """partial trace of density matrix rho on subsystems with dims, keeping indices in keep"""
    n = len(dims)
    r = rho.reshape(dims + dims)
    idx = list(range(n))
    trace_out = [i for i in idx if i not in keep]
    # move to einsum
    letters = 'abcdefghijklmnopqrstuvwxyz'
    inp_l = [letters[i] for i in range(n)]
    inp_r = [letters[i + n] for i in range(n)]
    for i in trace_out:
        inp_r[i] = inp_l[i]
    out = ''.join(inp_l[i] for i in keep) + ''.join(inp_r[i] for i in keep)
    res = np.einsum(''.join(inp_l) + ''.join(inp_r) + '->' + out, r)
    d = int(np.prod([dims[i] for i in keep]))
    return res.reshape(d, d)

def vn(rho, base=np.e):
    w = np.linalg.eigvalsh((rho + rho.conj().T) / 2)
    w = w[w > 1e-14]
    return float(-(w * np.log(w)).sum() / np.log(base))

def mutual_info(rho_ab, da, db):
    ra = ptrace(rho_ab, [da, db], [0]); rb = ptrace(rho_ab, [da, db], [1])
    return vn(ra) + vn(rb) - vn(rho_ab)

def random_projector(rng, d, rank=None):
    rank = rank or rng.integers(1, d)
    M = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
    q, _ = np.linalg.qr(M)
    v = q[:, :rank]
    return v @ v.conj().T

def negativity(rho, da, db):
    r = rho.reshape(da, db, da, db).transpose(0, 3, 2, 1).reshape(da * db, da * db)
    w = np.linalg.eigvalsh(r)
    return float(-w[w < 0].sum())

def chord(u, v):
    return 2 * np.abs(np.sin((v - u) / 2))

def xr_circle(a1, a2, b1, b2):
    """cross-ratio x(A,B) for arcs A=(a1,a2), B=(b1,b2) on the circle (same convention as exploration/checks2.py)"""
    return chord(a1, a2) * chord(b1, b2) / (chord(a1, b1) * chord(a2, b2))

def zeta(u, a1, a2):
    """A's modular (Rindler) coordinate for a point u outside A=(a1,a2) on the unit circle"""
    return np.log(np.abs(np.sin((u - a1) / 2) / np.sin((u - a2) / 2)))

def ff_entropy(sites):
    """entanglement entropy (nats) of a set of sites in the half-filled infinite free-fermion chain"""
    s = np.asarray(sites)
    d = s[:, None] - s[None, :]
    with np.errstate(divide='ignore', invalid='ignore'):
        C = np.where(d == 0, 0.5, np.sin(np.pi * d / 2) / (np.pi * d))
    w = np.linalg.eigvalsh(C)
    w = np.clip(w, 1e-15, 1 - 1e-15)
    return float(-(w * np.log(w) + (1 - w) * np.log(1 - w)).sum())
