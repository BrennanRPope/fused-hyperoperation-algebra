"""DRAFT shared helpers for checks.py and figures.py."""
import os, re
import numpy as np
from tbpy import *
HERE = os.path.dirname(os.path.abspath(__file__))

def homographic(e, euler=False):
    """Lagrange (equilateral) or Euler (collinear, body 2 in the middle) homographic solution, released at its far point."""
    if euler:
        c = np.array([[-1, 0], [0, 0], [1, 0]], float); w0 = np.sqrt(5 / 4); K = 1.25
    else:
        s = 1 / np.sqrt(3); c = np.array([[s * np.cos(a), s * np.sin(a)] for a in np.pi / 2 + 2 * np.pi / 3 * np.arange(3)]); w0 = np.sqrt(3.0); K = 1 / np.sqrt(3)
    v = np.sqrt(1 - e) * w0 * np.stack([-c[:, 1], c[:, 0]], 1)
    X0 = np.concatenate([c.ravel(), v.ravel()])
    r = np.linalg.norm(c[0]); E1 = 0.5 * (v[0] ** 2).sum() - K / r; a = -K / (2 * E1)
    return X0, 2 * np.pi * np.sqrt(a ** 3 / K), np.sqrt(-2 * E1)

def circfit(V):
    A = np.column_stack([2 * V[:, 0], 2 * V[:, 1], np.ones(len(V))]); b = (V ** 2).sum(1)
    sol = np.linalg.lstsq(A, b, rcond=None)[0]; c = sol[:2]; rho = np.sqrt(sol[2] + c @ c)
    return c, rho, np.abs(np.linalg.norm(V - c, axis=1) - rho).max() / rho

def fock(V, p0):
    """Fock/Moser stereographic map of planar momenta onto the unit 3-sphere (its z-slice is S^2)."""
    p2 = (V ** 2).sum(1)
    return np.column_stack([2 * p0 * V[:, 0], 2 * p0 * V[:, 1], np.zeros(len(V)), p2 - p0 ** 2]) / (p2 + p0 ** 2)[:, None]

def solid_angle(W, n=np.array([0.3, 0.5, 0.81]) / np.linalg.norm([0.3, 0.5, 0.81])):
    A = W; B = np.roll(W, -1, 0)
    num = np.cross(A, B) @ n; den = 1 + A @ n + (A * B).sum(1) + B @ n
    return 2 * np.arctan2(num, den).sum()

SUV = {"figure-eight": (0.347111, 0.532728, 6.324449), "butterfly I": (0.306893, 0.125507, 6.235641),
       "butterfly II": (0.392955, 0.097579, 7.003707), "moth I": (0.464445, 0.396060, 14.893671),
       "butterfly III": (0.405916, 0.230163, 13.865888), "moth II": (0.439166, 0.452968, 28.670972),
       "moth III": (0.383444, 0.377364, 25.840631), "goggles": (0.083300, 0.127889, 10.466818),
       "dragonfly": (0.080584, 0.588836, 21.270975), "yin-yang Ia": (0.513938, 0.304736, 17.328370),
       "yin-yang Ib": (0.282699, 0.327209, 10.962563), "bumblebee": (0.184279, 0.587188, 63.534541),
       "yarn": (0.559064, 0.349192, 55.501762)}

# catalog of free-group words (Li & Liao 2017, 695 equal-mass families), canonicalized
_wf = os.path.join(HERE, "data", "li-liao-2017-words.md")
if not os.path.exists(_wf):  # Li & Liao's public list (github.com/sjtu-liao/three-body, three-body-free-group-word.md)
    import urllib.request; os.makedirs(os.path.dirname(_wf), exist_ok=True)
    urllib.request.urlretrieve("https://raw.githubusercontent.com/sjtu-liao/three-body/main/three-body-free-group-word.md", _wf)
s = open(_wf).read()
CAT = {}
for th, w in re.findall(r"<th>(.*?)</th>\s*<td>(.*?)</td>", s, re.S):
    CAT.setdefault(canon(w.strip()), []).append(re.sub(r"<[^>]+>", "", th).replace("i.c.", "").strip())
def primitive(w):
    for k in range(1, len(w) + 1):
        if len(w) % k == 0 and w[:k] * (len(w) // k) == w: return w[:k], len(w) // k
    return w, 1
def lookup(w):
    if canon(w) in CAT: return CAT[canon(w)][0], 1
    u, k = primitive(w)
    if k > 1 and canon(u) in CAT: return CAT[canon(u)][0], k
    return None, 1

