import pytest

from snarkpy.field import GF, GFElement
from snarkpy.polynomial import Polynomial, SparsePolynomial


class TestPolynomial:
    def setup_method(self) -> None:
        self.p = 65537
        self.gf = GF(self.p, 2**17)

    def test_init(self) -> None:
        with pytest.raises(ValueError):
            Polynomial(self.gf)

    def test_call(self) -> None:
        poly = Polynomial(self.gf, coeffs=[1, 2, 3])
        assert poly(2) == self.gf(1 * 2**0 + 2 * 2**1 + 3 * 2**2)
        evals = [1, 2, 3, 4]
        poly = Polynomial(self.gf, evals=evals)
        omega = self.gf.nth_root_of_unity(len(evals))
        for i in range(len(evals)):
            assert poly(omega**i) == self.gf(evals[i])

    def test_repr(self) -> None:
        poly = Polynomial(self.gf, coeffs=[1, 2, 3, 4])
        assert repr(poly) == "1 + 2 * x^1 + 3 * x^2 + 4 * x^3 in F_65537"

    def test_str(self) -> None:
        poly = Polynomial(self.gf, coeffs=[1, 2, 3, 4])
        assert str(poly) == "1 + 2 * x^1 + 3 * x^2 + 4 * x^3"

    def test_len(self) -> None:
        poly = Polynomial(self.gf, coeffs=[1, 2, 3])
        assert len(poly) == 4

    def test_add(self) -> None:
        poly1 = Polynomial(self.gf, evals=[1, 2, 3, 4])
        poly2 = Polynomial(self.gf, coeffs=[1, 1, 1, 1])
        actual = poly1 + poly2
        assert actual.evals == [self.gf(5), self.gf(2), self.gf(3), self.gf(4)]

    def test_sub(self) -> None:
        poly1 = Polynomial(self.gf, evals=[1, 2, 3, 4])
        poly2 = Polynomial(self.gf, coeffs=[1, 1, 1, 1])
        actual = poly1 - poly2
        assert actual.evals == [self.gf(-3), self.gf(2), self.gf(3), self.gf(4)]

    def test_mul(self) -> None:
        poly1 = Polynomial(self.gf, coeffs=[1, 2])
        poly2 = Polynomial(self.gf, coeffs=[2, 3])
        actual = poly1 * poly2
        assert actual.coeffs == [
            self.gf(2),
            self.gf(7),
            self.gf(6),
            self.gf(0),
        ]

    def test_get_item(self) -> None:
        poly = Polynomial(self.gf, evals=[4, 0, 0, 0])
        assert poly[0] == self.gf(1)
        assert poly[1] == self.gf(1)
        assert poly[2] == self.gf(1)
        assert poly[3] == self.gf(1)
        assert poly[100] == self.gf(0)

    def test_degree(self) -> None:
        poly = Polynomial(self.gf, coeffs=[1, 2, 3, 0])
        assert poly.degree() == 2

    def test_extend(self) -> None:
        poly = Polynomial(self.gf, coeffs=[1, 2])
        actual = poly.extend(4)
        assert actual.coeffs == [self.gf(1), self.gf(2), self.gf(0), self.gf(0)]
        poly = Polynomial(self.gf, evals=[2, 0])
        actual = poly.extend(4)
        assert actual.coeffs == [self.gf(1), self.gf(1), self.gf(0), self.gf(0)]

    def test_divmod(self) -> None:
        f = Polynomial(self.gf, coeffs=[-1, 0, 0, 1])
        g = SparsePolynomial(self.gf, coeffs=[(0, -1), (1, 1)])
        actual = divmod(f, g)
        assert actual[0].coeffs == [
            self.gf(1),
            self.gf(1),
            self.gf(1),
            self.gf(0),
        ]
        assert actual[1].coeffs == [
            self.gf(0),
        ]
        f = Polynomial(self.gf, coeffs=[-1, 0, 1])
        g = SparsePolynomial(self.gf, coeffs=[(0, 1), (1, 1)])
        actual = divmod(f, g)
        assert actual[0].coeffs == [
            self.gf(-1),
            self.gf(1),
        ]
        assert actual[1].coeffs == [
            self.gf(0),
        ]
        f = Polynomial(self.gf, coeffs=[8, 5, 5, 4])
        g = SparsePolynomial(self.gf, coeffs=[(0, 1), (1, 4)])
        actual = divmod(f, g)
        assert actual[0].coeffs == [
            self.gf(1),
            self.gf(1),
            self.gf(1),
            self.gf(0),
        ]
        assert actual[1].coeffs == [
            self.gf(7),
        ]
        f = Polynomial(self.gf, coeffs=[5, 0, 2, 1])
        g = SparsePolynomial(self.gf, coeffs=[(0, 4), (1, -1), (2, 1)])
        actual = divmod(f, g)
        assert actual[0].coeffs == [
            self.gf(3),
            self.gf(1),
        ]
        assert actual[1].coeffs == [
            self.gf(-7),
            self.gf(-1),
        ]

    def test_get_evals(self) -> None:
        poly = Polynomial(self.gf, coeffs=[1, 1, 1, 1])
        assert poly.evals == [self.gf(4), self.gf(0), self.gf(0), self.gf(0)]

    def test_get_coeffs(self) -> None:
        poly = Polynomial(self.gf, evals=[4, 0, 0, 0])
        assert poly.coeffs == [self.gf(1), self.gf(1), self.gf(1), self.gf(1)]


class TestSparsePolynomial:
    def setup_method(self) -> None:
        self.p = 65537
        self.gf = GF(self.p, 2**17)

    def test_init(self) -> None:
        poly = SparsePolynomial(self.gf, coeffs=[(5, 10), (1, 20), (0, 30)])
        assert poly.coeffs == [(0, self.gf(30)), (1, self.gf(20)), (5, self.gf(10))]

    def test_repr(self) -> None:
        poly = SparsePolynomial(self.gf, coeffs=[(0, 10), (1, 20), (5, 30)])
        assert repr(poly) == "10 + 20 * x^1 + 30 * x^5 in F_65537"

    def test_degree(self) -> None:
        poly = SparsePolynomial(self.gf, coeffs=[(0, 10), (1, 20), (5, 30)])
        assert poly.degree() == 5
