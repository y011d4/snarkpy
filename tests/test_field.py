from snarkpy.field import GF, GFElement


class TestGF:
    def setup_method(self) -> None:
        self.p = 65537
        self.gf = GF(self.p, 2**17)

    def test_call(self) -> None:
        assert self.gf(3) == GFElement(3, self.gf, is_montgomery=False)

    def test_n8(self) -> None:
        assert self.gf.n8 == (self.p.bit_length() + 7) // 8

    def test_from_bytes(self) -> None:
        assert self.gf.from_bytes(b"\xff\x01", is_montgomery=False) == self.gf(
            0x01FF % self.p
        )
        assert self.gf.from_bytes(b"\xff\xff\xff\xff", is_montgomery=False) == self.gf(
            0xFFFFFFFF % self.p
        )

    def test_to_bytes(self) -> None:
        assert self.gf.to_bytes(self.gf(11)) == b"\x0b\x00\x00"

    def test_nth_root_of_unity(self) -> None:
        for i in range(16):
            assert self.gf.nth_root_of_unity(2**i) ** (2**i) == self.gf(1)


class TestGFElement:
    def setup_method(self) -> None:
        self.p = 65537
        self.gf = GF(self.p, 2**17)

    def test_add(self) -> None:
        assert self.gf(3) + self.gf(4) == self.gf(7)
        assert self.gf(3) + self.gf(-1) == self.gf(2)
        assert self.gf(3) + self.gf(self.p - 1) == self.gf(2)

    def test_sub(self) -> None:
        assert self.gf(3) - self.gf(4) == self.gf(-1)
        assert self.gf(3) - self.gf(-1) == self.gf(4)
        assert self.gf(3) - self.gf(self.p - 1) == self.gf(4)

    def test_mul(self) -> None:
        assert self.gf(3) * self.gf(4) == self.gf(12)
        assert self.gf(3) * self.gf(-1) == self.gf(-3)
        assert self.gf(3) * self.gf(self.p - 1) == self.gf(-3)

    def test_invert(self) -> None:
        assert ~self.gf(3) * self.gf(3) == self.gf(1)
        assert ~self.gf(2) == self.gf(32769)

    def test_truediv(self) -> None:
        assert self.gf(1) / self.gf(2) == self.gf(32769)
        assert self.gf(2) / self.gf(2) == self.gf(1)

    def test_pow(self) -> None:
        assert self.gf(3) ** 3 == self.gf(27)
        assert self.gf(3) ** 0 == self.gf(1)
        assert self.gf(2) ** -1 == self.gf(32769)
        assert self.gf(2) ** 17 == self.gf(pow(2, 17, self.p))

    def test_int(self) -> None:
        assert int(self.gf(3)) == 3
        assert int(self.gf(-1)) == self.p - 1

    def test_eq(self) -> None:
        self.gf2 = GF(self.p, 2**17)
        assert id(self.gf2) != id(self.gf)
        assert self.gf2(3) == self.gf(3)

    def test_repr(self) -> None:
        assert repr(self.gf(3)) == "3 in F_65537"

    def test_str(self) -> None:
        assert str(self.gf(3)) == "3"
