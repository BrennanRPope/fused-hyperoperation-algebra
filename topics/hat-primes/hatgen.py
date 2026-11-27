"""Hat tiling generator: a Python port of the metatile construction in
Smith, Myers, Kaplan, Goodman-Strauss, "An aperiodic monotile" (2023), as laid out
in Kaplan's hatviz demo. Pure numpy. Each hat is returned as (affine transform,
kind) where kind in {'H','H1','T','P','F'}; 'H1' is the one reflected hat in each H.

Coordinates: the hat outline is built on the triangular (Eisenstein) lattice
hexPt(x, y) = x + y*e^{i pi/3}; every placed hat is a rigid motion (rotation by a
multiple of 60 degrees, maybe a reflection) of the half-size outline.
"""
import numpy as np
from math import pi, cos, sin, sqrt

hr3 = sqrt(3) / 2

def pt(x, y): return np.array([x, y], dtype=float)
def hexPt(x, y): return pt(x + 0.5 * y, hr3 * y)

HAT = [hexPt(0, 0), hexPt(-1, -1), hexPt(0, -2), hexPt(2, -2), hexPt(2, -1), hexPt(4, -2),
       hexPt(5, -1), hexPt(4, 0), hexPt(3, 0), hexPt(2, 2), hexPt(0, 3), hexPt(0, 2), hexPt(-1, 2)]
HAT_HEX = [(0, 0), (-1, -1), (0, -2), (2, -2), (2, -1), (4, -2), (5, -1), (4, 0), (3, 0), (2, 2),
           (0, 3), (0, 2), (-1, 2)]

IDENT = np.array([1, 0, 0, 0, 1, 0], dtype=float)

def mul(A, B):
    return np.array([A[0]*B[0] + A[1]*B[3], A[0]*B[1] + A[1]*B[4], A[0]*B[2] + A[1]*B[5] + A[2],
                     A[3]*B[0] + A[4]*B[3], A[3]*B[1] + A[4]*B[4], A[3]*B[2] + A[4]*B[5] + A[5]])

def inv(T):
    det = T[0]*T[4] - T[1]*T[3]
    return np.array([T[4]/det, -T[1]/det, (T[1]*T[5] - T[2]*T[4])/det,
                     -T[3]/det, T[0]/det, (T[2]*T[3] - T[0]*T[5])/det])

def trot(ang):
    c, s = cos(ang), sin(ang)
    return np.array([c, -s, 0, s, c, 0])

def ttrans(tx, ty): return np.array([1, 0, tx, 0, 1, ty], dtype=float)

def rotAbout(p, ang): return mul(ttrans(p[0], p[1]), mul(trot(ang), ttrans(-p[0], -p[1])))

def transPt(M, P): return pt(M[0]*P[0] + M[1]*P[1] + M[2], M[3]*P[0] + M[4]*P[1] + M[5])

def matchSeg(p, q): return np.array([q[0]-p[0], p[1]-q[1], p[0], q[1]-p[1], q[0]-p[0], p[1]])

def matchTwo(p1, q1, p2, q2): return mul(matchSeg(p2, q2), inv(matchSeg(p1, q1)))

def intersect(p1, q1, p2, q2):
    d = (q2[1]-p2[1])*(q1[0]-p1[0]) - (q2[0]-p2[0])*(q1[1]-p1[1])
    uA = ((q2[0]-p2[0])*(p1[1]-p2[1]) - (q2[1]-p2[1])*(p1[0]-p2[0])) / d
    return pt(p1[0] + uA*(q1[0]-p1[0]), p1[1] + uA*(q1[1]-p1[1]))


class HatTile:
    def __init__(self, label): self.label = label
    def hats(self, T): yield (T, self.label)


class MetaTile:
    def __init__(self, shape, width):
        self.shape = [np.array(p) for p in shape]
        self.width = width
        self.children = []  # list of (T, geom)

    def addChild(self, T, geom): self.children.append([np.array(T, dtype=float), geom])

    def evalChild(self, n, i):
        T, g = self.children[n]
        return transPt(T, g.shape[i])

    def recentre(self):
        c = sum(self.shape) / len(self.shape)
        self.shape = [p - c for p in self.shape]
        M = ttrans(-c[0], -c[1])
        for ch in self.children:
            ch[0] = mul(M, ch[0])

    def hats(self, T):
        for S, g in self.children:
            yield from g.hats(mul(T, S))


H1_HAT, H_HAT, T_HAT, P_HAT, F_HAT = (HatTile(s) for s in ('H1', 'H', 'T', 'P', 'F'))


def initH():
    o = [pt(0, 0), pt(4, 0), pt(4.5, hr3), pt(2.5, 5*hr3), pt(1.5, 5*hr3), pt(-0.5, hr3)]
    m = MetaTile(o, 2)
    m.addChild(matchTwo(HAT[5], HAT[7], o[5], o[0]), H_HAT)
    m.addChild(matchTwo(HAT[9], HAT[11], o[1], o[2]), H_HAT)
    m.addChild(matchTwo(HAT[5], HAT[7], o[3], o[4]), H_HAT)
    m.addChild(mul(ttrans(2.5, hr3), mul(np.array([-0.5, -hr3, 0, hr3, -0.5, 0]),
                                          np.array([0.5, 0, 0, 0, -0.5, 0]))), H1_HAT)
    return m

def initT():
    o = [pt(0, 0), pt(3, 0), pt(1.5, 3*hr3)]
    m = MetaTile(o, 2)
    m.addChild(np.array([0.5, 0, 0.5, 0, 0.5, hr3]), T_HAT)
    return m

def initP():
    o = [pt(0, 0), pt(4, 0), pt(3, 2*hr3), pt(-1, 2*hr3)]
    m = MetaTile(o, 2)
    m.addChild(np.array([0.5, 0, 1.5, 0, 0.5, hr3]), P_HAT)
    m.addChild(mul(ttrans(0, 2*hr3), mul(np.array([0.5, hr3, 0, -hr3, 0.5, 0]),
                                          np.array([0.5, 0, 0, 0, 0.5, 0]))), P_HAT)
    return m

def initF():
    o = [pt(0, 0), pt(3, 0), pt(3.5, hr3), pt(3, 2*hr3), pt(-1, 2*hr3)]
    m = MetaTile(o, 2)
    m.addChild(np.array([0.5, 0, 1.5, 0, 0.5, hr3]), F_HAT)
    m.addChild(mul(ttrans(0, 2*hr3), mul(np.array([0.5, hr3, 0, -hr3, 0.5, 0]),
                                          np.array([0.5, 0, 0, 0, 0.5, 0]))), F_HAT)
    return m


RULES = [['H'], [0, 0, 'P', 2], [1, 0, 'H', 2], [2, 0, 'P', 2], [3, 0, 'H', 2], [4, 4, 'P', 2],
         [0, 4, 'F', 3], [2, 4, 'F', 3], [4, 1, 3, 2, 'F', 0], [8, 3, 'H', 0], [9, 2, 'P', 0],
         [10, 2, 'H', 0], [11, 4, 'P', 2], [12, 0, 'H', 2], [13, 0, 'F', 3], [14, 2, 'F', 1],
         [15, 3, 'H', 4], [8, 2, 'F', 1], [17, 3, 'H', 0], [18, 2, 'P', 0], [19, 2, 'H', 2],
         [20, 4, 'F', 3], [20, 0, 'P', 2], [22, 0, 'H', 2], [23, 4, 'F', 3], [23, 0, 'F', 3],
         [16, 0, 'P', 2], [9, 4, 0, 2, 'T', 2], [4, 0, 'F', 3]]
PATCH_KINDS = [r[-2] if len(r) > 1 else r[0] for r in RULES]

def constructPatch(H, T, P, F):
    shapes = {'H': H, 'T': T, 'P': P, 'F': F}
    ret = MetaTile([], H.width)
    for r in RULES:
        if len(r) == 1:
            ret.addChild(IDENT, shapes[r[0]])
        elif len(r) == 4:
            Tm, g = ret.children[r[0]]
            poly = g.shape
            P_ = transPt(Tm, poly[(r[1]+1) % len(poly)])
            Q_ = transPt(Tm, poly[r[1]])
            ns = shapes[r[2]]; npoly = ns.shape
            ret.addChild(matchTwo(npoly[r[3]], npoly[(r[3]+1) % len(npoly)], P_, Q_), ns)
        else:
            TP, gP = ret.children[r[0]]
            TQ, gQ = ret.children[r[2]]
            P_ = transPt(TQ, gQ.shape[r[3]])
            Q_ = transPt(TP, gP.shape[r[1]])
            ns = shapes[r[4]]; npoly = ns.shape
            ret.addChild(matchTwo(npoly[r[5]], npoly[(r[5]+1) % len(npoly)], P_, Q_), ns)
    return ret

def padd(a, b): return a + b
def psub(a, b): return a - b

NEW_H_CH = [0, 9, 16, 27, 26, 6, 1, 8, 10, 15]
NEW_P_CH = [7, 2, 3, 4, 28]
NEW_F_CH = [21, 20, 22, 23, 24, 25]
NEW_T_CH = [11]

def constructMetatiles(patch):
    bps1 = patch.evalChild(8, 2)
    bps2 = patch.evalChild(21, 2)
    rbps = transPt(rotAbout(bps1, -2.0*pi/3.0), bps2)
    p72 = patch.evalChild(7, 2)
    p252 = patch.evalChild(25, 2)
    llc = intersect(bps1, rbps, patch.evalChild(6, 2), p72)
    w = psub(patch.evalChild(6, 2), llc)
    Ho = [llc, bps1]
    w = transPt(rotAbout(pt(0, 0), -pi/3), w)
    Ho.append(padd(Ho[1], w))
    Ho.append(patch.evalChild(14, 2))
    w = transPt(rotAbout(pt(0, 0), -pi/3), w)
    Ho.append(psub(Ho[3], w))
    Ho.append(patch.evalChild(6, 2))
    nH = MetaTile(Ho, patch.width*2)
    for c in NEW_H_CH: nH.addChild(*patch.children[c])
    Po = [p72, padd(p72, psub(bps1, llc)), bps1, llc]
    nP = MetaTile(Po, patch.width*2)
    for c in NEW_P_CH: nP.addChild(*patch.children[c])
    Fo = [bps2, patch.evalChild(24, 2), patch.evalChild(25, 0), p252, padd(p252, psub(llc, bps1))]
    nF = MetaTile(Fo, patch.width*2)
    for c in NEW_F_CH: nF.addChild(*patch.children[c])
    AAA = Ho[2]
    BBB = padd(Ho[1], psub(Ho[4], Ho[5]))
    CCC = transPt(rotAbout(BBB, -pi/3), AAA)
    nT = MetaTile([BBB, CCC, AAA], patch.width*2)
    nT.addChild(*patch.children[11])
    for m in (nH, nP, nF, nT): m.recentre()
    return [nH, nT, nP, nF]


def levels(n):
    """Return list of (H,T,P,F) metatile quadruples for levels 0..n."""
    tiles = [initH(), initT(), initP(), initF()]
    out = [tiles]
    for _ in range(n):
        tiles = constructMetatiles(constructPatch(*tiles))
        out.append(tiles)
    return out


def hats_of(geom):
    """List of (T, kind) for every hat in a metatile or patch."""
    return list(geom.hats(IDENT))


def hat_polygon(T):
    return np.array([transPt(T, p) for p in HAT])


def to_hex(P, scale=0.5):
    """Exact triangular-lattice coordinates (x, y) with P = scale*hexPt(x, y)."""
    y = P[1] / (hr3 * scale)
    x = P[0] / scale - 0.5 * y
    return x, y
