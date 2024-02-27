from snarkpy.elliptic_curve import bn128, BN128Element


class TestBN128:
    def test_bn128(self) -> None:
        bn128.fp


class TestBN128Element:
    def test_add(self) -> None:
        two = BN128Element(
            1368015179489954701390400359078579693043519447331113978918064868415326638035,
            9918110051302171585080402603319702774565515993150576347155970296011118125764,
        )
        three = BN128Element(
            3353031288059533942658390886683067124040920775575537747144343083137631628272,
            19321533766552368860946552437480515441416830039777911637913418824951667761761,
        )
        assert bn128.g + bn128.g == two
        assert bn128.g + two == three
        assert two + bn128.g == three

    def test_mul(self) -> None:
        three = BN128Element(
            3353031288059533942658390886683067124040920775575537747144343083137631628272,
            19321533766552368860946552437480515441416830039777911637913418824951667761761,
        )
        assert bn128.g * 3 == three
        assert bn128.g * bn128.order == BN128Element.zero()

    def test_eq(self) -> None:
        assert BN128Element(1, 2) == BN128Element(3, 6, 3)

    def test_repr(self) -> None:
        assert repr(bn128.g) == "BN128Element(1, 2, 1)"

    def test_str(self) -> None:
        assert str(bn128.g) == "BN128Element(1, 2)"
