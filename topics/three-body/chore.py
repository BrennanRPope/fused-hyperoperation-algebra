"""DRAFT: action minimizer for three-body choreographies (one curve, three bodies a third of a period apart),
in the inertial frame (omega = 0) or a frame turning at rate omega (relative choreographies).
Period fixed at 2*pi. Harmonics divisible by 3 are left out so the centre of mass stays put."""
import numpy as np
from scipy.optimize import minimize

class Chore:
    def __init__(self, N=24, M=720, omega=0.0):
        self.ns = np.array([n for n in range(1, N + 1) if n % 3]); self.M = M; self.omega = omega
        t = np.arange(M) * 2 * np.pi / M; self.t = t; self.dt = 2 * np.pi / M
        n = self.ns[None, :]
        self.C = np.hstack([np.cos(n * t[:, None]), np.sin(n * t[:, None])])
        self.Cd = np.hstack([-n * np.sin(n * t[:, None]), n * np.cos(n * t[:, None])])
        t2 = t + 2 * np.pi / 3
        self.C2 = np.hstack([np.cos(n * t2[:, None]), np.sin(n * t2[:, None])])
        self.nc = self.C.shape[1]

    def unpack(self, x): return x.reshape(self.nc, 2)

    def action(self, x):
        c = self.unpack(x); w = self.omega
        Q = self.C @ c; Qd = self.Cd @ c
        V = Qd + w * np.column_stack([-Q[:, 1], Q[:, 0]])
        D = Q - self.C2 @ c; r = np.sqrt((D ** 2).sum(1))
        A = 3 * self.dt * (0.5 * (V ** 2).sum() + (1 / r).sum())
        gV = V; gQ = w * np.column_stack([V[:, 1], -V[:, 0]])
        gD = -D / r[:, None] ** 3
        g = 3 * self.dt * (self.Cd.T @ gV + self.C.T @ gQ + (self.C - self.C2).T @ gD)
        return A, g.ravel()

    def minimize(self, x0):
        r = minimize(self.action, x0, jac=True, method="BFGS", options={"gtol": 1e-10, "maxiter": 20000})
        return r.x, r.fun, np.linalg.norm(r.jac)

    def initial_state(self, x):
        c = self.unpack(x); w = self.omega
        X = np.zeros(12)
        for k in range(3):
            s = 2 * np.pi * k / 3; n = self.ns
            b = np.concatenate([np.cos(n * s), np.sin(n * s)]); bd = np.concatenate([-n * np.sin(n * s), n * np.cos(n * s)])
            q = b @ c; qd = bd @ c
            X[2 * k:2 * k + 2] = q; X[6 + 2 * k:8 + 2 * k] = qd + w * np.array([-q[1], q[0]])
        return X

    def curve(self, x, M=1000):
        c = self.unpack(x); t = np.linspace(0, 2 * np.pi, M); n = self.ns[None, :]
        B = np.hstack([np.cos(n * t[:, None]), np.sin(n * t[:, None])]); return B @ c
