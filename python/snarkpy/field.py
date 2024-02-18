use_rust = True
if use_rust:
    from snarkpy._snarkpy import field

    GF = field.GF
    GFElement = field.GFElement
else:
    from functools import cache
    from math import ceil
    from typing import Union

    class GF:
        def __init__(self, p: int, R: int) -> None:
            assert R & (R - 1) == 0, "R must be power of 2"
            assert R > p, "R must be greater than p"
            self.p = p
            self.R = R
            self.R_bits = R.bit_length() - 1
            self.R_mask = R - 1
            self.R2 = pow(R, 2, p)
            g, _, p_inv = xgcd(R, p)
            assert g == 1
            self.p_prime = -p_inv

        def add(self, a: int, b: int) -> int:
            return (a + b) % self.p

        def sub(self, a: int, b: int) -> int:
            return (a - b) % self.p

        def mul(self, a: int, b: int) -> int:
            return self.from_montgomery(a * b)

        def inv(self, a: int) -> int:
            return self.pow(a, self.p - 2)

        def div(self, a: int, b: int) -> int:
            return self.mul(a, self.inv(b))

        def pow(self, a: int, b: int) -> int:
            b %= self.p - 1
            bits = [int(bb) for bb in bin(b)[2:][::-1]]
            # bits += [0] * (self.p.bit_length() - len(bits))
            a0 = self.to_montgomery(1)
            a1 = a
            # a0 = a
            # a1 = self.mul(a, a)
            for b in bits[::-1]:
                if b == 0:
                    a1 = self.mul(a0, a1)
                    a0 = self.mul(a0, a0)
                else:
                    a0 = self.mul(a0, a1)
                    a1 = self.mul(a1, a1)
            return a0

        def one(self) -> "GFElement":
            return GFElement(1, self, False)

        def __call__(self, value: Union[int, "GFElement"]) -> "GFElement":
            if isinstance(value, int):
                return GFElement(value, self, False)
            elif isinstance(value, GFElement):
                return value
            else:
                raise RuntimeError

        def __repr__(self) -> str:
            return f"F_{self.p}"

        def __str__(self) -> str:
            return f"F_{self.p}"

        def to_montgomery(self, value: int) -> int:
            return self.from_montgomery(value * self.R2)

        def from_montgomery(self, value: int) -> int:
            t = (value + ((value * self.p_prime) & self.R_mask) * self.p) >> self.R_bits
            if t >= self.p:
                return t - self.p
            else:
                return t - 0

        @property
        def n8(self) -> int:
            return ceil(self.p.bit_length() / 8)

        def from_bytes(self, value: bytes, is_montgomery: bool) -> "GFElement":
            return GFElement(int.from_bytes(value, "little"), self, is_montgomery)

        def to_bytes(self, value: "GFElement") -> bytes:
            return self.from_montgomery(value._value).to_bytes(self.n8, "little")

        def __contains__(self, value: "GFElement") -> bool:
            if not isinstance(value, GFElement):
                return False
            return value._gf == self

        def __eq__(self, other: object) -> bool:
            if not isinstance(other, GF):
                return False
            return self.p == other.p and self.R == other.R

    class GFElement:
        def __init__(self, value: int, gf: GF, is_montgomery: bool = True) -> None:
            self._gf = gf
            if is_montgomery:
                self._value = value
            else:
                self._value = gf.to_montgomery(value)

        def __add__(self, other: "GFElement") -> "GFElement":
            return GFElement(self._gf.add(self._value, other._value), self._gf)

        def __sub__(self, other: "GFElement") -> "GFElement":
            return GFElement(self._gf.sub(self._value, other._value), self._gf)

        def __mul__(self, other: "GFElement") -> "GFElement":
            return GFElement(self._gf.mul(self._value, other._value), self._gf)

        def __truediv__(self, other: "GFElement") -> "GFElement":
            return GFElement(self._gf.div(self._value, other._value), self._gf)

        def __pow__(self, other: int) -> "GFElement":
            return GFElement(self._gf.pow(self._value, other), self._gf)

        def __int__(self) -> int:
            return self._gf.from_montgomery(self._value)

        def __eq__(self, other: object) -> bool:
            if not isinstance(other, GFElement):
                return False
            return self._value == other._value

        def __repr__(self) -> str:
            return f"{int(self)} in {self._gf}"

        def __str__(self) -> str:
            return str(int(self))

        def inv(self) -> "GFElement":
            return GFElement(self._gf.inv(self._value), self._gf)

    def xgcd(a: int, b: int) -> tuple[int, int, int]:
        if b == 0:
            return a, 1, 0
        else:
            d, x, y = xgcd(b, a % b)
            return d, y, x - (a // b) * y

    def primitive_root() -> int:
        p = 21888242871839275222246405745257275088548364400416034343698204186575808495617
        factors = [
            2,
            3,
            13,
            29,
            983,
            11003,
            237073,
            405928799,
            1670836401704629,
            13818364434197438864469338081,
        ]
        k = 2
        while True:
            for f in factors:
                if pow(k, (p - 1) // f, p) == 1:
                    break
            else:
                return k
            k += 1

    @cache
    def calc_omega(power: int) -> int:
        p = 21888242871839275222246405745257275088548364400416034343698204186575808495617
        assert (p - 1) % 2**power == 0, "power is invalid (must be power <= 28)"
        # k = primitive_root()
        k = 5
        return pow(k, (p - 1) // 2**power, p)
