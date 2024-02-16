import json
from dataclasses import dataclass
from math import log2
from pathlib import Path
from typing import Sequence, Union

from py_ecc import optimized_bn128 as bn128
from py_ecc.optimized_bn128 import FQ, FQ2

from snarkpy.field import GF, GFElement
from snarkpy.polynomial import Polynomial
from snarkpy.file import Reader, Sections


@dataclass
class ZKey:
    protocol: str
    n8q: int
    q: int
    n8r: int
    r: int
    n_vars: int
    n_public: int
    domain_size: int
    power: int
    n_additions: int
    n_constraints: int
    k1: GFElement
    k2: GFElement
    qm: tuple[FQ, FQ, FQ]
    ql: tuple[FQ, FQ, FQ]
    qr: tuple[FQ, FQ, FQ]
    qo: tuple[FQ, FQ, FQ]
    qc: tuple[FQ, FQ, FQ]
    s1: tuple[FQ, FQ, FQ]
    s2: tuple[FQ, FQ, FQ]
    s3: tuple[FQ, FQ, FQ]
    x_2: tuple[FQ2, FQ2, FQ2]
    additions_buf: bytes
    a_idx: Sequence[int]
    b_idx: Sequence[int]
    c_idx: Sequence[int]
    poly_qm: Polynomial
    poly_qm_4: Polynomial
    poly_ql: Polynomial
    poly_ql_4: Polynomial
    poly_qr: Polynomial
    poly_qr_4: Polynomial
    poly_qo: Polynomial
    poly_qo_4: Polynomial
    poly_qc: Polynomial
    poly_qc_4: Polynomial
    poly_s1: Polynomial
    poly_s1_4: Polynomial
    poly_s2: Polynomial
    poly_s2_4: Polynomial
    poly_s3: Polynomial
    poly_s3_4: Polynomial
    poly_ls: Sequence[Polynomial]
    poly_l_4s: Sequence[Polynomial]
    taus: Sequence[tuple[FQ, FQ, FQ]]

    def __post_init__(self) -> None:
        assert bn128.is_on_curve(self.qm, bn128.b)
        assert bn128.is_on_curve(self.ql, bn128.b)
        assert bn128.is_on_curve(self.qr, bn128.b)
        assert bn128.is_on_curve(self.qo, bn128.b)
        assert bn128.is_on_curve(self.qc, bn128.b)
        assert bn128.is_on_curve(self.s1, bn128.b)
        assert bn128.is_on_curve(self.s2, bn128.b)
        assert bn128.is_on_curve(self.s3, bn128.b)
        assert bn128.is_on_curve(self.x_2, bn128.b2)


class Wtns:
    _wtns_buf: bytes
    _internal_wtns_buf: bytes
    _zkey: ZKey
    poly_a: Polynomial
    poly_b: Polynomial
    poly_c: Polynomial

    def __init__(
        self,
        wtns_buf: bytes,
        internal_wtns_buf: bytes,
        zkey: ZKey,
        poly_a: Polynomial,
        poly_b: Polynomial,
        poly_c: Polynomial,
    ) -> None:
        self._wtns_buf = wtns_buf
        self._internal_wtns_buf = internal_wtns_buf
        self._zkey = zkey
        self.poly_a = poly_a
        self.poly_b = poly_b
        self.poly_c = poly_c

    def get_witness(self, idx: int) -> GFElement:
        return _get_witness(idx, self._zkey, self._wtns_buf, self._internal_wtns_buf)


def parse_zkey(filename: Path) -> ZKey:
    reader = Reader(filename)
    sections = Sections.from_reader(reader, b"zkey")
    section = sections.get_section(1)
    reader.move_to_section(section)
    protocol_id = reader.read_uint32()
    assert reader.check_section_end(section)
    assert protocol_id == 2, "only support PLONK protocol"

    section = sections.get_section(2)
    reader.move_to_section(section)
    n8q = reader.read_uint32()
    q = reader.read_uintn(n8q)
    n8r = reader.read_uint32()
    r = reader.read_uintn(n8r)
    assert q == bn128.field_modulus
    curve = bn128
    n_vars = reader.read_uint32()
    n_public = reader.read_uint32()
    domain_size = reader.read_uint32()
    power = int(log2(domain_size))
    n_additions = reader.read_uint32()
    n_constraints = reader.read_uint32()
    Fr = GF(r, 2**256)
    Fq = GF(q, 2**256)
    k1 = Fr.from_bytes(reader.read(n8r), is_montgomery=True)
    k2 = Fr.from_bytes(reader.read(n8r), is_montgomery=True)
    # これらは tau で評価したときの値
    qm = (
        x := FQ(int(Fq.from_bytes(reader.read(n8r), is_montgomery=True))),
        y := FQ(int(Fq.from_bytes(reader.read(n8r), is_montgomery=True))),
        FQ(0) if x == FQ.zero() and y == FQ.zero() else FQ(1),
    )
    ql = (
        x := FQ(int(Fq.from_bytes(reader.read(n8r), is_montgomery=True))),
        y := FQ(int(Fq.from_bytes(reader.read(n8r), is_montgomery=True))),
        FQ(0) if x == FQ.zero() and y == FQ.zero() else FQ(1),
    )
    qr = (
        x := FQ(int(Fq.from_bytes(reader.read(n8r), is_montgomery=True))),
        y := FQ(int(Fq.from_bytes(reader.read(n8r), is_montgomery=True))),
        FQ(0) if x == FQ.zero() and y == FQ.zero() else FQ(1),
    )
    qo = (
        x := FQ(int(Fq.from_bytes(reader.read(n8r), is_montgomery=True))),
        y := FQ(int(Fq.from_bytes(reader.read(n8r), is_montgomery=True))),
        FQ(0) if x == FQ.zero() and y == FQ.zero() else FQ(1),
    )
    qc = (
        x := FQ(int(Fq.from_bytes(reader.read(n8r), is_montgomery=True))),
        y := FQ(int(Fq.from_bytes(reader.read(n8r), is_montgomery=True))),
        FQ(0) if x == FQ.zero() and y == FQ.zero() else FQ(1),
    )
    s1 = (
        x := FQ(int(Fq.from_bytes(reader.read(n8r), is_montgomery=True))),
        y := FQ(int(Fq.from_bytes(reader.read(n8r), is_montgomery=True))),
        FQ(0) if x == FQ.zero() and y == FQ.zero() else FQ(1),
    )
    s2 = (
        x := FQ(int(Fq.from_bytes(reader.read(n8r), is_montgomery=True))),
        y := FQ(int(Fq.from_bytes(reader.read(n8r), is_montgomery=True))),
        FQ(0) if x == FQ.zero() and y == FQ.zero() else FQ(1),
    )
    s3 = (
        x := FQ(int(Fq.from_bytes(reader.read(n8r), is_montgomery=True))),
        y := FQ(int(Fq.from_bytes(reader.read(n8r), is_montgomery=True))),
        FQ(0) if x == FQ.zero() and y == FQ.zero() else FQ(1),
    )
    x_2 = (
        FQ2(
            (
                int(Fq.from_bytes(reader.read(n8r), is_montgomery=True)),
                int(Fq.from_bytes(reader.read(n8r), is_montgomery=True)),
            )
        ),
        FQ2(
            (
                int(Fq.from_bytes(reader.read(n8r), is_montgomery=True)),
                int(Fq.from_bytes(reader.read(n8r), is_montgomery=True)),
            )
        ),
        FQ2((1, 0)),
    )

    additions_section = sections.get_section(3)
    reader.move_to_section(additions_section)
    additions_buf = reader.read(additions_section.size)
    assert reader.check_section_end(additions_section)

    def section_to_idx(section_id: int) -> Sequence[int]:
        section = sections.get_section(section_id)
        reader.move_to_section(section)
        idx = []
        for _ in range(n_constraints):
            idx.append(reader.read_uint32())
        assert reader.check_section_end(section)
        return idx

    def section_to_poly(
        section_id: int, continued: bool = False, check: bool = True
    ) -> tuple[Polynomial, Polynomial]:
        section = sections.get_section(section_id)
        if not continued:
            reader.move_to_section(section)
        coeffs = []
        for _ in range(domain_size):
            coeffs.append(Fr.from_bytes(reader.read(n8r), is_montgomery=True))
        poly = Polynomial(r, coeffs=coeffs)
        evals_4 = []
        for _ in range(domain_size * 4):
            evals_4.append(Fr.from_bytes(reader.read(n8r), is_montgomery=True))
        poly_4 = Polynomial(r, evals=evals_4)
        if check:
            assert reader.check_section_end(section)
        return poly, poly_4

    a_idx = section_to_idx(4)
    b_idx = section_to_idx(5)
    c_idx = section_to_idx(6)
    poly_qm, poly_qm_4 = section_to_poly(7)
    poly_ql, poly_ql_4 = section_to_poly(8)
    poly_qr, poly_qr_4 = section_to_poly(9)
    poly_qo, poly_qo_4 = section_to_poly(10)
    poly_qc, poly_qc_4 = section_to_poly(11)
    poly_s1, poly_s1_4 = section_to_poly(12, False, False)
    poly_s2, poly_s2_4 = section_to_poly(12, True, False)
    poly_s3, poly_s3_4 = section_to_poly(12, True, True)
    poly_ls = []
    poly_l_4s = []
    for i in range(n_public):
        if i == 0:
            poly_l, poly_l_4 = section_to_poly(13, False, False)
        elif 0 < i < n_public - 1:
            poly_l, poly_l_4 = section_to_poly(13, True, False)
        else:
            poly_l, poly_l_4 = section_to_poly(13, True, True)
        poly_ls.append(poly_l)
        poly_l_4s.append(poly_l_4)

    tau_section = sections.get_section(14)
    reader.move_to_section(tau_section)
    taus = []
    for _ in range(tau_section.size // 64):
        x = Fq.from_bytes(reader.read(32), is_montgomery=True)
        y = Fq.from_bytes(reader.read(32), is_montgomery=True)
        z = Fq(0) if x == Fq(0) and y == Fq(0) else Fq(1)
        taus.append((FQ(int(x)), FQ(int(y)), FQ(int(z))))
    assert reader.check_section_end(tau_section)

    return ZKey(
        protocol="plonk",
        n8q=n8q,
        q=q,
        n8r=n8r,
        r=r,
        n_vars=n_vars,
        n_public=n_public,
        domain_size=domain_size,
        power=power,
        n_additions=n_additions,
        n_constraints=n_constraints,
        k1=Fr(k1),
        k2=Fr(k2),
        qm=qm,
        ql=ql,
        qr=qr,
        qo=qo,
        qc=qc,
        s1=s1,
        s2=s2,
        s3=s3,
        x_2=x_2,
        additions_buf=additions_buf,
        a_idx=a_idx,
        b_idx=b_idx,
        c_idx=c_idx,
        poly_qm=poly_qm,
        poly_qm_4=poly_qm_4,
        poly_ql=poly_ql,
        poly_ql_4=poly_ql_4,
        poly_qr=poly_qr,
        poly_qr_4=poly_qr_4,
        poly_qo=poly_qo,
        poly_qo_4=poly_qo_4,
        poly_qc=poly_qc,
        poly_qc_4=poly_qc_4,
        poly_s1=poly_s1,
        poly_s1_4=poly_s1_4,
        poly_s2=poly_s2,
        poly_s2_4=poly_s2_4,
        poly_s3=poly_s3,
        poly_s3_4=poly_s3_4,
        poly_ls=poly_ls,
        poly_l_4s=poly_l_4s,
        taus=taus,
    )


def _get_witness(
    idx: int,
    zkey: ZKey,
    wtns_buf: Union[bytes, Sequence[int]],
    internal_wtns_buf: Union[bytes, Sequence[int]],
) -> GFElement:
    Fr = GF(zkey.r, 2**256)
    if idx < zkey.n_vars - zkey.n_additions:
        return Fr.from_bytes(
            bytes(wtns_buf[idx * zkey.n8r : (idx + 1) * zkey.n8r]),
            is_montgomery=False,
        )
    elif idx < zkey.n_vars:
        return Fr.from_bytes(
            bytes(
                internal_wtns_buf[
                    (idx - (zkey.n_vars - zkey.n_additions))
                    * zkey.n8r : (idx - (zkey.n_vars - zkey.n_additions))
                    * zkey.n8r
                    + zkey.n8r
                ]
            ),
            is_montgomery=False,
        )
    else:
        return Fr(0)


def parse_wtns(filename: Path, zkey: ZKey) -> Wtns:
    reader = Reader(filename)
    sections = Sections.from_reader(reader, b"wtns")
    wtns_section = sections.get_section(2)
    reader.move_to_section(wtns_section)
    wtns_buf = reader.read(wtns_section.size)
    assert reader.check_section_end(wtns_section)
    wtns_buf = b"\x00" * 32 + wtns_buf[32:]

    Fr = GF(zkey.r, 2**256)
    s_sum = 8 + Fr.n8 * 2
    internal_wtns_list: list[int] = [-1] * (Fr.n8 * zkey.n_additions)
    for i in range(zkey.n_additions):
        ai = int.from_bytes(zkey.additions_buf[i * s_sum : i * s_sum + 4], "little")
        bi = int.from_bytes(zkey.additions_buf[i * s_sum + 4 : i * s_sum + 8], "little")
        ac = Fr.from_bytes(
            zkey.additions_buf[i * s_sum + 8 : i * s_sum + 8 + Fr.n8],
            is_montgomery=True,
        )
        bc = Fr.from_bytes(
            zkey.additions_buf[i * s_sum + 8 + Fr.n8 : i * s_sum + 8 + Fr.n8 * 2],
            is_montgomery=True,
        )
        aw = _get_witness(ai, zkey, wtns_buf, internal_wtns_list)
        bw = _get_witness(bi, zkey, wtns_buf, internal_wtns_list)
        r = ac * aw + bc * bw
        internal_wtns_list[i * Fr.n8 : (i + 1) * Fr.n8] = Fr.to_bytes(r)
    assert all([n != -1 for n in internal_wtns_list])
    internal_wtns_buf = bytes(internal_wtns_list)

    def idx_list_to_poly(idx_list: Sequence[int]) -> Polynomial:
        evals = []
        for idx in idx_list:
            evals.append(_get_witness(idx, zkey, wtns_buf, internal_wtns_buf))
        evals += [Fr(0)] * (zkey.domain_size - len(evals))
        return Polynomial(zkey.r, evals=evals)

    poly_a = idx_list_to_poly(zkey.a_idx)
    poly_b = idx_list_to_poly(zkey.b_idx)
    poly_c = idx_list_to_poly(zkey.c_idx)
    return Wtns(
        wtns_buf=wtns_buf,
        internal_wtns_buf=bytes(internal_wtns_buf),
        zkey=zkey,
        poly_a=poly_a,
        poly_b=poly_b,
        poly_c=poly_c,
    )


@dataclass
class VKey:
    protocol: str
    curve: str
    n_public: int
    power: int
    k1: GFElement
    k2: GFElement
    qm: tuple[FQ, FQ, FQ]
    ql: tuple[FQ, FQ, FQ]
    qr: tuple[FQ, FQ, FQ]
    qo: tuple[FQ, FQ, FQ]
    qc: tuple[FQ, FQ, FQ]
    s1: tuple[FQ, FQ, FQ]
    s2: tuple[FQ, FQ, FQ]
    s3: tuple[FQ, FQ, FQ]
    x_2: tuple[FQ2, FQ2, FQ2]
    omega: GFElement

    def __post_init__(self) -> None:
        assert self.protocol == "plonk", "only support PLONK protocol"
        assert self.curve == "bn128", "only support BN128 curve"
        assert bn128.is_on_curve(self.qm, bn128.b)
        assert bn128.is_on_curve(self.ql, bn128.b)
        assert bn128.is_on_curve(self.qr, bn128.b)
        assert bn128.is_on_curve(self.qo, bn128.b)
        assert bn128.is_on_curve(self.qc, bn128.b)
        assert bn128.is_on_curve(self.s1, bn128.b)
        assert bn128.is_on_curve(self.s2, bn128.b)
        assert bn128.is_on_curve(self.s3, bn128.b)
        assert bn128.is_on_curve(self.x_2, bn128.b2)


def parse_vkey(filename: Path) -> VKey:
    with filename.open() as fp:
        data = json.load(fp)
    protocol = data["protocol"]
    assert protocol == "plonk", "only support PLONK protocol"
    curve = data["curve"]
    assert curve == "bn128", "only support BN128 curve"
    n_public = data["nPublic"]
    poewr = data["power"]
    r = bn128.curve_order
    Fr = GF(r, 2**256)
    k1 = Fr(int(data["k1"]))
    k2 = Fr(int(data["k2"]))
    qm = tuple(FQ(int(d)) for d in data["Qm"])
    ql = tuple(FQ(int(d)) for d in data["Ql"])
    qr = tuple(FQ(int(d)) for d in data["Qr"])
    qo = tuple(FQ(int(d)) for d in data["Qo"])
    qc = tuple(FQ(int(d)) for d in data["Qc"])
    s1 = tuple(FQ(int(d)) for d in data["S1"])
    s2 = tuple(FQ(int(d)) for d in data["S2"])
    s3 = tuple(FQ(int(d)) for d in data["S3"])
    x_2 = tuple(FQ2((int(d[0]), int(d[1]))) for d in data["X_2"])
    assert len(qm) == 3
    assert len(ql) == 3
    assert len(qr) == 3
    assert len(qo) == 3
    assert len(qc) == 3
    assert len(s1) == 3
    assert len(s2) == 3
    assert len(s3) == 3
    assert len(x_2) == 3
    omega = Fr(int(data["w"]))
    return VKey(
        protocol=protocol,
        curve=curve,
        n_public=n_public,
        power=poewr,
        k1=k1,
        k2=k2,
        qm=qm,
        ql=ql,
        qr=qr,
        qo=qo,
        qc=qc,
        s1=s1,
        s2=s2,
        s3=s3,
        x_2=x_2,
        omega=omega,
    )
