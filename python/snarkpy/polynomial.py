from math import ceil, log2
from typing import Optional, Sequence, TypeVar, Union

from snarkpy.field import GF, GFElement


class Polynomial:
    Number = TypeVar("Number", int, GFElement)

    _evals: Optional[Sequence[GFElement]]
    _coeffs: Optional[Sequence[GFElement]]

    def __init__(
        self,
        p: int,
        coeffs: Optional[Sequence[Number]] = None,
        evals: Optional[Sequence[Number]] = None,
    ) -> None:
        assert coeffs is not None or evals is not None
        self._p = p
        self._Fp = GF(p, 2**256)
        if coeffs is not None:
            l = len(coeffs)
            assert l > 0
            n = ceil(log2(l))
            coeffs_fp = [self._Fp(c) for c in list(coeffs) + [0] * (2**n - l)]
            self._coeffs = coeffs_fp
            # self._evals = self._coeffs_to_evals(coeffs_fp)
            self._evals = None
        else:
            assert evals is not None
            assert len(evals) & (len(evals) - 1) == 0, "n must be power of 2"
            evals_fp = [self._Fp(e) for e in evals]
            self._evals = evals_fp
            # self._coeffs = self._evals_to_coeffs(evals_fp)
            self._coeffs = None

    def __call__(self, x: Number) -> GFElement:
        if self._coeffs is None:
            assert self._evals is not None
            self._coeffs = self._evals_to_coeffs(self._evals)
        Fp = self._Fp
        x_fp = Fp(x)
        res = Fp(0)
        for c in self._coeffs[::-1]:
            res = res * x_fp
            res += c
        return res

    def __repr__(self) -> str:
        if self._coeffs is None:
            assert self._evals is not None
            self._coeffs = self._evals_to_coeffs(self._evals)
        Fp = self._Fp
        res = []
        for i, c in enumerate(self._coeffs):
            if c == Fp(0):
                continue
            if i == 0:
                res.append(str(c))
            elif c == Fp(1):
                res.append(f"x^{i}")
            else:
                res.append(f"{c} * x^{i}")
        ret = " + ".join(res)
        return ret[:1024] + ("" if len(ret) < 1024 else " ...") + f" in {self._Fp}"

    def __len__(self) -> int:
        if self._coeffs is not None:
            return len(self._coeffs)
        if self._evals is not None:
            return len(self._evals)
        raise RuntimeError

    def _coeffs_to_evals(self, coeffs: Sequence[GFElement]) -> Sequence[GFElement]:
        Fp = self._Fp
        n = len(coeffs)
        assert n & (n - 1) == 0, "n must be power of 2"
        omegas = [Fp.nth_root_of_unity(2**i) for i in range(int(log2(n)) + 1)]
        return self._fft(coeffs, omegas)

    def _fft(
        self, coeffs: Sequence[GFElement], omegas: Sequence[GFElement]
    ) -> Sequence[GFElement]:
        Fp = self._Fp
        n = len(coeffs)
        assert n & (n - 1) == 0, "n must be power of 2"

        def fft(coeffs: Sequence[GFElement]) -> Sequence[GFElement]:
            n = len(coeffs)
            assert n & (n - 1) == 0, "n must be power of 2"
            if n == 1:
                return coeffs
            ye = fft(coeffs[0::2])
            yo = fft(coeffs[1::2])
            omega = omegas[int(log2(n))]
            omega_i = Fp(1)
            evals = [Fp(0)] * len(coeffs)
            for i in range(n // 2):
                tmp = omega_i * yo[i]
                evals[i] = ye[i] + tmp
                evals[i + n // 2] = ye[i] - tmp
                omega_i = omega_i * omega
            return evals

        return fft(coeffs)

    def _evals_to_coeffs(self, evals: Sequence[GFElement]) -> Sequence[GFElement]:
        Fp = self._Fp
        n = len(evals)
        assert n & (n - 1) == 0, "n must be power of 2"
        ninv = Fp(n) ** -1
        omega_invs = [Fp.nth_root_of_unity(2**i) ** -1 for i in range(int(log2(n)) + 1)]
        return [c * ninv for c in self._fft(evals, omega_invs)]

    def _calc_evals_if_necessary(self) -> None:
        if self._evals is None:
            assert self._coeffs is not None
            self._evals = self._coeffs_to_evals(self._coeffs)

    def _calc_coeffs_if_necessary(self) -> None:
        if self._coeffs is None:
            assert self._evals is not None
            self._coeffs = self._evals_to_coeffs(self._evals)

    def extend(self, n: int) -> "Polynomial":
        assert n & (n - 1) == 0, "n must be power of 2"
        self._calc_coeffs_if_necessary()
        assert self._coeffs is not None
        assert n > len(self._coeffs)
        coeffs = list(self._coeffs) + [self._Fp(0)] * (n - len(self._coeffs))
        return Polynomial(self._p, coeffs=coeffs)

    def _prepare_operation(
        self, a: "Polynomial", b: "Polynomial"
    ) -> tuple["Polynomial", "Polynomial"]:
        if len(a) < len(b):
            a, b = b, a
        if len(a) > len(b):
            b = b.extend(len(a))
        a._calc_evals_if_necessary()
        b._calc_evals_if_necessary()
        return a, b

    def __add__(self, other: Union["Polynomial", GFElement]) -> "Polynomial":
        if isinstance(other, Polynomial):
            a, b = self._prepare_operation(self, other)
            assert a._evals is not None
            assert b._evals is not None
            return Polynomial(
                self._p,
                evals=[e1 + e2 for e1, e2 in zip(a._evals, b._evals)],
            )
        elif isinstance(other, GFElement):
            assert self._p == other.gf.p
            if self._coeffs is not None:
                return Polynomial(
                    self._p, coeffs=[self._coeffs[0] + other] + list(self._coeffs[1:])
                )
            elif self._evals is not None:
                return self + Polynomial(self._p, coeffs=[other])
            else:
                raise RuntimeError
        else:
            raise RuntimeError

    def __sub__(self, other: Union["Polynomial", GFElement]) -> "Polynomial":
        if isinstance(other, Polynomial):
            a, b = self._prepare_operation(self, other)
            assert a._evals is not None
            assert b._evals is not None
            return Polynomial(
                self._p,
                evals=[e1 - e2 for e1, e2 in zip(a._evals, b._evals)],
            )
        elif isinstance(other, GFElement):
            assert self._p == other.gf.p
            if self._coeffs is not None:
                return Polynomial(
                    self._p, coeffs=[self._coeffs[0] - other] + list(self._coeffs[1:])
                )
            elif self._evals is not None:
                return self - Polynomial(self._p, coeffs=[other])
            else:
                raise RuntimeError
        else:
            raise RuntimeError

    def __mul__(self, other: Union["Polynomial", GFElement]) -> "Polynomial":
        if isinstance(other, Polynomial):
            a, b = self._prepare_operation(self, other)
            assert a._evals is not None
            assert b._evals is not None
            return Polynomial(
                self._p,
                evals=[e1 * e2 for e1, e2 in zip(a._evals, b._evals)],
            )
        elif isinstance(other, GFElement):
            assert self._p == other.gf.p
            if self._coeffs is not None:
                return Polynomial(self._p, coeffs=[c * other for c in self._coeffs])
            elif self._evals is not None:
                return Polynomial(self._p, evals=[e * other for e in self._evals])
            else:
                raise RuntimeError
        else:
            raise RuntimeError

    def __getitem__(self, idx: int) -> GFElement:
        self._calc_coeffs_if_necessary()
        assert self._coeffs is not None
        return self._coeffs[idx]

    def degree(self) -> int:
        self._calc_coeffs_if_necessary()
        assert self._coeffs is not None
        i = len(self._coeffs) - 1
        while True:
            if self._coeffs[i] != self._Fp(0):
                return i
            i -= 1
            if i < 0:
                return 0

    @property
    def coeffs(self) -> Sequence[GFElement]:
        self._calc_coeffs_if_necessary()
        assert self._coeffs is not None
        return self._coeffs

    @property
    def evals(self) -> Sequence[GFElement]:
        self._calc_evals_if_necessary()
        assert self._evals is not None
        return self._evals
