"""Tests for tesseract.py. Run: python3 -m unittest test_tesseract -v"""

import random
import unittest
from fractions import Fraction
from math import comb, exp

import tesseract as t


def rand_frac(rng):
    return Fraction(rng.randint(-9, 9), rng.randint(1, 5))


class Numbers(unittest.TestCase):
    def test_multadd_is_box_minus_one(self):
        rng = random.Random(1)
        for _ in range(200):
            x, y = rand_frac(rng), rand_frac(rng)
            self.assertEqual(1 + t.multadd(x, y), (1 + x) * (1 + y))

    def test_associative_commutative_identity(self):
        rng = random.Random(2)
        for _ in range(200):
            x, y, z = (rand_frac(rng) for _ in range(3))
            self.assertEqual(t.multadd(t.multadd(x, y), z), t.multadd(x, t.multadd(y, z)))
            self.assertEqual(t.multadd(x, y), t.multadd(y, x))
            self.assertEqual(t.multadd(x, 0), x)

    def test_inverse(self):
        for x in [Fraction(1, 2), Fraction(3), Fraction(-1, 2), Fraction(7, 3)]:
            self.assertEqual(t.multadd(x, t.multadd_inverse(x)), 0)
        with self.assertRaises(ZeroDivisionError):
            t.multadd_inverse(-1)


class Corners(unittest.TestCase):
    def test_tesseract_layers_are_pascal(self):
        layers = t.corners_by_grade(4)
        self.assertEqual([len(L) for L in layers], [1, 4, 6, 4, 1])
        for n in range(1, 9):
            self.assertEqual([len(L) for L in t.corners_by_grade(n)],
                             [comb(n, k) for k in range(n + 1)])

    def test_labels(self):
        self.assertEqual(t.label(0b0101, 4), "1010")
        self.assertEqual(t.name(0b0101, 4), "ac")
        self.assertEqual(t.name(0, 4), "1")
        self.assertEqual(sorted(t.subsets(0b101)), [0, 1, 4, 5])


class Mapping(unittest.TestCase):
    def test_corner_sum_is_one_plus_multadd(self):
        rng = random.Random(3)
        for n in range(1, 7):
            a = [rand_frac(rng) for _ in range(n)]
            self.assertEqual(t.total(t.from_point(a)), 1 + t.multadd_all(a))

    def test_round_trip(self):
        a = [Fraction(1, 2), 2, 3, Fraction(-1, 3)]
        c = t.from_point(a)
        self.assertEqual(t.to_point(c), a)
        self.assertTrue(t.is_pure(c))
        c[0b1111] += 1  # nudge the top corner: no longer a product
        self.assertFalse(t.is_pure(c))

    def test_evaluate_all_ones_is_box(self):
        rng = random.Random(4)
        x = [rand_frac(rng) for _ in range(4)]
        ones = [1] * 16
        box = 1
        for xi in x:
            box *= 1 + xi
        self.assertEqual(t.evaluate(ones, x), box)

    def test_grade_totals_are_elementary_symmetric(self):
        a = [1, 2, 3, 4]
        self.assertEqual(t.grade_totals(t.from_point(a)), [1, 10, 35, 50, 24])


class Transforms(unittest.TestCase):
    def test_zeta_gives_box_sides(self):
        a = [Fraction(1, 2), 2, 3, Fraction(-1, 3)]
        b = t.zeta(t.from_point(a))
        for mask in range(16):
            expect = 1
            for i in range(4):
                if mask >> i & 1:
                    expect *= 1 + a[i]
            self.assertEqual(b[mask], expect)

    def test_mobius_inverts_zeta(self):
        rng = random.Random(5)
        for n in range(0, 6):
            c = [rand_frac(rng) for _ in range(1 << n)]
            self.assertEqual(t.mobius(t.zeta(c)), c)
            self.assertEqual(t.zeta(t.mobius(c)), c)

    def test_compose_is_axiswise_multadd_for_pure_points(self):
        rng = random.Random(6)
        for n in range(1, 6):
            a = [rand_frac(rng) for _ in range(n)]
            b = [rand_frac(rng) for _ in range(n)]
            ab = [t.multadd(x, y) for x, y in zip(a, b)]
            self.assertEqual(t.compose(t.from_point(a), t.from_point(b)), t.from_point(ab))

    def test_compose_general_points(self):
        rng = random.Random(7)
        for n in range(0, 5):
            v = [rand_frac(rng) for _ in range(1 << n)]
            w = [rand_frac(rng) for _ in range(1 << n)]
            u = [rand_frac(rng) for _ in range(1 << n)]
            self.assertEqual(t.compose(v, w), t.compose_fast(v, w))
            self.assertEqual(t.compose(v, w), t.compose(w, v))
            self.assertEqual(t.compose(t.compose(v, w), u), t.compose(v, t.compose(w, u)))
            unit = [1] + [0] * ((1 << n) - 1)
            self.assertEqual(t.compose(v, unit), v)


class Soliton(unittest.TestCase):
    def test_two_soliton_tau(self):
        A12 = Fraction(1, 9)
        w = t.soliton_tau_weights([[1, A12], [A12, 1]])
        self.assertEqual(w, [1, 1, 1, A12])  # 1 + a + b + A ab
        e1, e2 = 0.3, -1.1
        tau = t.evaluate(w, [exp(e1), exp(e2)])
        self.assertAlmostEqual(tau, 1 + exp(e1) + exp(e2) + float(A12) * exp(e1 + e2))

    def test_no_interaction_is_pure(self):
        ones = [[1] * 4 for _ in range(4)]
        self.assertTrue(t.is_pure(t.soliton_tau_weights(ones)))
        A = [[1, 2, 3, 5], [2, 1, 7, 11], [3, 7, 1, 13], [5, 11, 13, 1]]
        w = t.soliton_tau_weights(A)
        self.assertEqual(w[0b1111], 2 * 3 * 5 * 7 * 11 * 13)
        self.assertFalse(t.is_pure(w))


if __name__ == "__main__":
    unittest.main()
