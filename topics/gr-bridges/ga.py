"""Tiny geometric (Clifford) algebra with a diagonal metric, symbolic coefficients (sympy).

A multivector is a dict {bitmask: coefficient}; bit i set means basis vector e_i is present.
e_i e_i = k_i (the 'dial' along direction i); e_i e_j = -e_j e_i for i != j.
"""
import sympy as sp


class GA:
    def __init__(self, squares):
        self.k = list(squares)
        self.n = len(self.k)

    def _blade_mul(self, A, B):
        # sign from reordering, metric factor from repeated vectors
        s = 1
        a = A >> 1
        while a:
            s *= (-1) ** bin(a & B).count("1")
            a >>= 1
        f = 1
        common = A & B
        for i in range(self.n):
            if common >> i & 1:
                f *= self.k[i]
        return s * f, A ^ B

    def mul(self, X, Y):
        out = {}
        for A, a in X.items():
            for B, b in Y.items():
                c, C = self._blade_mul(A, B)
                out[C] = out.get(C, 0) + c * a * b
        return self.clean(out)

    @staticmethod
    def clean(X):
        return {A: sp.expand(c) for A, c in X.items() if sp.expand(c) != 0}

    def add(self, *Xs):
        out = {}
        for X in Xs:
            for A, a in X.items():
                out[A] = out.get(A, 0) + a
        return self.clean(out)

    def scale(self, s, X):
        return self.clean({A: s * a for A, a in X.items()})

    def sub(self, X, Y):
        return self.add(X, self.scale(-1, Y))

    def e(self, *idx):
        X = {0: 1}
        for i in idx:
            X = self.mul(X, {1 << i: 1})
        return X

    def scalar(self, s):
        return {0: s} if s != 0 else {}

    def vector(self, comps):
        return self.clean({1 << i: c for i, c in enumerate(comps)})

    def general(self, name):
        """A general multivector with symbolic coefficients name_0 ... name_{2^n-1}."""
        return {A: sp.Symbol(f"{name}_{A}") for A in range(2 ** self.n)}

    def comm(self, X, Y):
        return self.sub(self.mul(X, Y), self.mul(Y, X))

    def equal(self, X, Y):
        return self.sub(X, Y) == {}

    def grade(self, X, g):
        return {A: c for A, c in X.items() if bin(A).count("1") == g}
