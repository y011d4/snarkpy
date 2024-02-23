import pytest

from snarkpy.field import GF, GFElement
from snarkpy.polynomial import Polynomial


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

    def test_calc_evals_if_necessary(self) -> None:
        poly = Polynomial(self.gf, coeffs=[1, 1, 1, 1])
        assert poly.evals is None
        poly.calc_evals_if_necessary()
        assert poly.evals == [self.gf(4), self.gf(0), self.gf(0), self.gf(0)]

    def test_calc_coeffs_if_necessary(self) -> None:
        poly = Polynomial(self.gf, evals=[4, 0, 0, 0])
        assert poly.coeffs is None
        poly.calc_coeffs_if_necessary()
        assert poly.coeffs == [self.gf(1), self.gf(1), self.gf(1), self.gf(1)]

    def test_repr(self) -> None:
        poly = Polynomial(self.gf, coeffs=[1, 2, 3, 4])
        assert repr(poly) == "1 + 2 * x^1 + 3 * x^2 + 4 * x^3 in F_65537"

    def test_len(self) -> None:
        poly = Polynomial(self.gf, coeffs=[1, 2, 3])
        assert len(poly) == 4

    def test_add(self) -> None:
        poly1 = Polynomial(self.gf, evals=[1, 2, 3, 4])
        poly2 = Polynomial(self.gf, coeffs=[1, 1, 1, 1])
        actual = poly1 + poly2
        assert actual.evals == [self.gf(5), self.gf(2), self.gf(3), self.gf(4)]
        assert actual.coeffs == None

    def test_sub(self) -> None:
        poly1 = Polynomial(self.gf, evals=[1, 2, 3, 4])
        poly2 = Polynomial(self.gf, coeffs=[1, 1, 1, 1])
        actual = poly1 - poly2
        assert actual.evals == [self.gf(-3), self.gf(2), self.gf(3), self.gf(4)]
        assert actual.coeffs == None

    def test_mul(self) -> None:
        poly1 = Polynomial(self.gf, evals=[1, 2, 3, 4])
        poly2 = Polynomial(self.gf, coeffs=[1, 1, 1, 1])
        actual = poly1 * poly2
        assert actual.evals == [self.gf(4), self.gf(0), self.gf(0), self.gf(0)]
        assert actual.coeffs == None

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
