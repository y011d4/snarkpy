import json
from dataclasses import asdict, dataclass
from pathlib import Path
from secrets import randbelow
from typing import Any, Mapping, Sequence

from py_ecc import optimized_bn128 as bn128
from py_ecc.optimized_bn128 import FQ

from snarkpy.elliptic import exp_tau
from snarkpy.field import GF, GFElement, calc_omega
from snarkpy.hash import keccak
from snarkpy.polynomial import Polynomial
from snarkpy.plonk.parse import parse_wtns, parse_zkey


@dataclass
class Proof:
    a: tuple[FQ, FQ, FQ]
    b: tuple[FQ, FQ, FQ]
    c: tuple[FQ, FQ, FQ]
    z: tuple[FQ, FQ, FQ]
    t_low: tuple[FQ, FQ, FQ]
    t_mid: tuple[FQ, FQ, FQ]
    t_high: tuple[FQ, FQ, FQ]
    w_zeta: tuple[FQ, FQ, FQ]
    w_zetaw: tuple[FQ, FQ, FQ]
    eval_a: GFElement
    eval_b: GFElement
    eval_c: GFElement
    eval_s1: GFElement
    eval_s2: GFElement
    eval_zw: GFElement
    protocol: str = "plonk"
    curve: str = "bn128"

    def __post_init__(self):
        assert self.protocol == "plonk"
        assert self.curve == "bn128"
        Fr = GF(bn128.curve_order, 2**256)
        assert bn128.is_on_curve(self.a, bn128.b)
        assert bn128.is_on_curve(self.b, bn128.b)
        assert bn128.is_on_curve(self.c, bn128.b)
        assert bn128.is_on_curve(self.z, bn128.b)
        assert bn128.is_on_curve(self.t_low, bn128.b)
        assert bn128.is_on_curve(self.t_mid, bn128.b)
        assert bn128.is_on_curve(self.t_high, bn128.b)
        assert bn128.is_on_curve(self.w_zeta, bn128.b)
        assert bn128.is_on_curve(self.w_zetaw, bn128.b)
        assert self.eval_a in Fr
        assert self.eval_b in Fr
        assert self.eval_c in Fr
        assert self.eval_s1 in Fr
        assert self.eval_s2 in Fr
        assert self.eval_zw in Fr

    def to_json(self) -> str:
        def dict_factory(items: Sequence[tuple[str, Any]]) -> Mapping[str, Any]:
            key_map = {
                "a": "A",
                "b": "B",
                "c": "C",
                "z": "Z",
                "t_low": "T1",
                "t_mid": "T2",
                "t_high": "T3",
                "w_zeta": "Wxi",
                "w_zetaw": "Wxiw",
            }
            d = {}
            for k, v in items:
                if isinstance(v, tuple):
                    v = list(map(str, v))
                elif isinstance(v, GFElement):
                    v = str(v)
                if k in key_map:
                    k = key_map[k]
                d[k] = v
            return d

        return json.dumps(asdict(self, dict_factory=dict_factory))

    @classmethod
    def from_json(cls, res: Mapping[str, Any]) -> "Proof":
        assert res["curve"] == "bn128"
        Fr = GF(bn128.curve_order, 2**256)
        return Proof(
            a=(FQ(int(res["A"][0])), FQ(int(res["A"][1])), FQ(int(res["A"][2]))),
            b=(FQ(int(res["B"][0])), FQ(int(res["B"][1])), FQ(int(res["B"][2]))),
            c=(FQ(int(res["C"][0])), FQ(int(res["C"][1])), FQ(int(res["C"][2]))),
            z=(FQ(int(res["Z"][0])), FQ(int(res["Z"][1])), FQ(int(res["Z"][2]))),
            t_low=(FQ(int(res["T1"][0])), FQ(int(res["T1"][1])), FQ(int(res["T1"][2]))),
            t_mid=(FQ(int(res["T2"][0])), FQ(int(res["T2"][1])), FQ(int(res["T2"][2]))),
            t_high=(
                FQ(int(res["T3"][0])),
                FQ(int(res["T3"][1])),
                FQ(int(res["T3"][2])),
            ),
            w_zeta=(
                FQ(int(res["Wxi"][0])),
                FQ(int(res["Wxi"][1])),
                FQ(int(res["Wxi"][2])),
            ),
            w_zetaw=(
                FQ(int(res["Wxiw"][0])),
                FQ(int(res["Wxiw"][1])),
                FQ(int(res["Wxiw"][2])),
            ),
            eval_a=Fr(int(res["eval_a"])),
            eval_b=Fr(int(res["eval_b"])),
            eval_c=Fr(int(res["eval_c"])),
            eval_s1=Fr(int(res["eval_s1"])),
            eval_s2=Fr(int(res["eval_s2"])),
            eval_zw=Fr(int(res["eval_zw"])),
            protocol=res["protocol"],
            curve=res["curve"],
        )


def prove(
    zkey_path: Path, wtns_path: Path, debug: bool = False
) -> tuple[Proof, Sequence[int]]:
    zkey = parse_zkey(zkey_path)
    wtns = parse_wtns(wtns_path, zkey)

    Fr = GF(zkey.r, 2**256)

    omega = Fr(calc_omega(zkey.power))

    zh_coeffs = [Fr(0)] * zkey.domain_size * 4
    zh_coeffs[0] = Fr(-1)
    zh_coeffs[zkey.domain_size] = Fr(1)
    poly_zh_4 = Polynomial(zkey.r, coeffs=zh_coeffs)

    # round 1
    b = [Fr(randbelow(zkey.r)) for _ in range(11)]
    poly_a = wtns.poly_a
    poly_b = wtns.poly_b
    poly_c = wtns.poly_c
    taus = zkey.taus
    poly_s1 = zkey.poly_s1
    poly_s2 = zkey.poly_s2
    poly_s3 = zkey.poly_s3
    poly_l_4s = zkey.poly_l_4s
    poly_qm_4 = zkey.poly_qm_4
    poly_ql_4 = zkey.poly_ql_4
    poly_qr_4 = zkey.poly_qr_4
    poly_qo_4 = zkey.poly_qo_4
    poly_qc_4 = zkey.poly_qc_4
    poly_s1_4 = zkey.poly_s1_4
    poly_s2_4 = zkey.poly_s2_4
    poly_s3_4 = zkey.poly_s3_4
    poly_qm = zkey.poly_qm
    poly_ql = zkey.poly_ql
    poly_qr = zkey.poly_qr
    poly_qo = zkey.poly_qo
    poly_qc = zkey.poly_qc
    poly_ls = zkey.poly_ls

    poly_a_blinded_4 = Polynomial(zkey.r, coeffs=[b[1], b[0]]) * poly_zh_4 + poly_a
    poly_b_blinded_4 = Polynomial(zkey.r, coeffs=[b[3], b[2]]) * poly_zh_4 + poly_b
    poly_c_blinded_4 = Polynomial(zkey.r, coeffs=[b[5], b[4]]) * poly_zh_4 + poly_c

    proof_a = exp_tau(poly_a_blinded_4, taus)
    proof_b = exp_tau(poly_b_blinded_4, taus)
    proof_c = exp_tau(poly_c_blinded_4, taus)

    # round 2
    transcript1 = b""
    transcript1 += int(bn128.normalize(zkey.qm)[0]).to_bytes(32, "big")
    transcript1 += int(bn128.normalize(zkey.qm)[1]).to_bytes(32, "big")
    transcript1 += int(bn128.normalize(zkey.ql)[0]).to_bytes(32, "big")
    transcript1 += int(bn128.normalize(zkey.ql)[1]).to_bytes(32, "big")
    transcript1 += int(bn128.normalize(zkey.qr)[0]).to_bytes(32, "big")
    transcript1 += int(bn128.normalize(zkey.qr)[1]).to_bytes(32, "big")
    transcript1 += int(bn128.normalize(zkey.qo)[0]).to_bytes(32, "big")
    transcript1 += int(bn128.normalize(zkey.qo)[1]).to_bytes(32, "big")
    transcript1 += int(bn128.normalize(zkey.qc)[0]).to_bytes(32, "big")
    transcript1 += int(bn128.normalize(zkey.qc)[1]).to_bytes(32, "big")
    transcript1 += int(bn128.normalize(zkey.s1)[0]).to_bytes(32, "big")
    transcript1 += int(bn128.normalize(zkey.s1)[1]).to_bytes(32, "big")
    transcript1 += int(bn128.normalize(zkey.s2)[0]).to_bytes(32, "big")
    transcript1 += int(bn128.normalize(zkey.s2)[1]).to_bytes(32, "big")
    transcript1 += int(bn128.normalize(zkey.s3)[0]).to_bytes(32, "big")
    transcript1 += int(bn128.normalize(zkey.s3)[1]).to_bytes(32, "big")
    for i in range(zkey.n_public):
        transcript1 += Fr.to_bytes(poly_a(omega**i))[::-1]
    transcript1 += int(proof_a[0]).to_bytes(32, "big")
    transcript1 += int(proof_a[1]).to_bytes(32, "big")
    transcript1 += int(proof_b[0]).to_bytes(32, "big")
    transcript1 += int(proof_b[1]).to_bytes(32, "big")
    transcript1 += int(proof_c[0]).to_bytes(32, "big")
    transcript1 += int(proof_c[1]).to_bytes(32, "big")
    h = keccak(transcript1)
    beta = Fr(int.from_bytes(h, "big"))
    transcript2 = Fr.to_bytes(beta)[::-1]
    h = keccak(transcript2)
    gamma = Fr(int.from_bytes(h, "big"))

    evals = [Fr(1)]
    w = Fr(1)
    for i in range(zkey.domain_size):
        res = evals[-1]
        res = res * (poly_a.evals[i] + beta * w + gamma)
        res = res * (poly_b.evals[i] + beta * zkey.k1 * w + gamma)
        res = res * (poly_c.evals[i] + beta * zkey.k2 * w + gamma)
        res = res / (poly_a.evals[i] + poly_s1.evals[i] * beta + gamma)
        res = res / (poly_b.evals[i] + poly_s2.evals[i] * beta + gamma)
        res = res / (poly_c.evals[i] + poly_s3.evals[i] * beta + gamma)
        evals.append(res)
        w = w * omega
    assert evals.pop() == Fr(1)
    poly_z = Polynomial(zkey.r, evals=evals)
    poly_z_blinded_4 = (
        Polynomial(zkey.r, coeffs=[b[8], b[7], b[6]]) * poly_zh_4 + poly_z
    )
    proof_z = exp_tau(poly_z_blinded_4, taus)

    # round 3
    transcript3 = b""
    transcript3 += Fr.to_bytes(beta)[::-1]
    transcript3 += Fr.to_bytes(gamma)[::-1]
    transcript3 += int(proof_z[0]).to_bytes(32, "big")
    transcript3 += int(proof_z[1]).to_bytes(32, "big")
    h = keccak(transcript3)
    alpha = Fr(int.from_bytes(h, "big"))

    pi_evals = [Fr(0) for _ in range(zkey.domain_size * 4)]
    for i in range(zkey.domain_size * 4):
        for j in range(zkey.n_public):
            pi_evals[i] = pi_evals[i] - poly_l_4s[j].evals[i] * poly_a.evals[j]
    poly_pi_4 = Polynomial(zkey.r, evals=pi_evals)
    poly_zw_blinded_4 = Polynomial(
        zkey.r, coeffs=[c * omega**i for i, c in enumerate(poly_z_blinded_4.coeffs)]
    )
    poly_z_blinded_8 = poly_z_blinded_4.extend(2 ** (zkey.power + 3))
    poly_zw_blinded_8 = poly_zw_blinded_4.extend(2 ** (zkey.power + 3))

    ta = (
        poly_a_blinded_4 * poly_b_blinded_4 * poly_qm_4
        + poly_a_blinded_4 * poly_ql_4
        + poly_b_blinded_4 * poly_qr_4
        + poly_c_blinded_4 * poly_qo_4
        + poly_pi_4
        + poly_qc_4
    )
    tb = (
        (poly_a_blinded_4 + Polynomial(zkey.r, coeffs=[gamma, beta]))
        * (poly_b_blinded_4 + Polynomial(zkey.r, coeffs=[gamma, beta * zkey.k1]))
        * (poly_c_blinded_4 + Polynomial(zkey.r, coeffs=[gamma, beta * zkey.k2]))
        * poly_z_blinded_8
        * alpha
    )
    tc = (
        (poly_a_blinded_4 + poly_s1_4 * beta + gamma)
        * (poly_b_blinded_4 + poly_s2_4 * beta + gamma)
        * (poly_c_blinded_4 + poly_s3_4 * beta + gamma)
        * poly_zw_blinded_8
        * alpha
    )
    td = (poly_z_blinded_4 + Fr(-1)) * poly_l_4s[0] * alpha * alpha
    tzh = ta + tb - tc + td

    t_coeffs = list(tzh.coeffs)
    for i in range(zkey.domain_size):
        t_coeffs[i] = t_coeffs[i] * Fr(-1)
    for i in range(zkey.domain_size, zkey.domain_size * 4 + 6):
        t_coeffs[i] = t_coeffs[i - zkey.domain_size] - t_coeffs[i]
        if i > zkey.domain_size * 3 + 5 and t_coeffs[i] != Fr(0):
            raise RuntimeError
    assert all(c == Fr(0) for c in t_coeffs[zkey.domain_size * 3 + 6 :])

    poly_t_low = Polynomial(zkey.r, coeffs=t_coeffs[: zkey.domain_size] + [b[9]])
    poly_t_mid = Polynomial(
        zkey.r,
        coeffs=[t_coeffs[zkey.domain_size] - b[9]]
        + t_coeffs[zkey.domain_size + 1 : zkey.domain_size * 2]
        + [b[10]],
    )
    poly_t_high = Polynomial(
        zkey.r,
        coeffs=[t_coeffs[zkey.domain_size * 2] - b[10]]
        + t_coeffs[zkey.domain_size * 2 + 1 : zkey.domain_size * 3 + 6],
    )
    proof_t1 = exp_tau(poly_t_low, taus)
    proof_t2 = exp_tau(poly_t_mid, taus)
    proof_t3 = exp_tau(poly_t_high, taus)

    # round 4
    transcript4 = b""
    transcript4 += Fr.to_bytes(alpha)[::-1]
    transcript4 += int(proof_t1[0]).to_bytes(32, "big")
    transcript4 += int(proof_t1[1]).to_bytes(32, "big")
    transcript4 += int(proof_t2[0]).to_bytes(32, "big")
    transcript4 += int(proof_t2[1]).to_bytes(32, "big")
    transcript4 += int(proof_t3[0]).to_bytes(32, "big")
    transcript4 += int(proof_t3[1]).to_bytes(32, "big")
    h = keccak(transcript4)
    zeta = Fr(int.from_bytes(h, "big"))

    eval_a = poly_a_blinded_4(zeta)
    eval_b = poly_b_blinded_4(zeta)
    eval_c = poly_c_blinded_4(zeta)
    eval_s1 = poly_s1(zeta)
    eval_s2 = poly_s2(zeta)
    eval_zw = poly_zw_blinded_4(zeta)

    # round 5
    transcript5 = b""
    transcript5 += Fr.to_bytes(zeta)[::-1]
    transcript5 += int(eval_a).to_bytes(32, "big")
    transcript5 += int(eval_b).to_bytes(32, "big")
    transcript5 += int(eval_c).to_bytes(32, "big")
    transcript5 += int(eval_s1).to_bytes(32, "big")
    transcript5 += int(eval_s2).to_bytes(32, "big")
    transcript5 += int(eval_zw).to_bytes(32, "big")
    h = keccak(transcript5)
    v = Fr(int.from_bytes(h, "big"))

    ra = (
        poly_qm * eval_a * eval_b
        + poly_ql * eval_a
        + poly_qr * eval_b
        + poly_qo * eval_c
        + poly_pi_4(zeta)
        + poly_qc
    )
    rb = (
        poly_z_blinded_4
        * (eval_a + beta * zeta + gamma)
        * (eval_b + beta * zkey.k1 * zeta + gamma)
        * (eval_c + beta * zkey.k2 * zeta + gamma)
        * alpha
    )
    rc = (
        (poly_s3 * beta + eval_c + gamma)
        * (eval_a + eval_s1 * beta + gamma)
        * (eval_b + eval_s2 * beta + gamma)
        * eval_zw
        * alpha
    )
    rd = (poly_z_blinded_4 + Fr(-1)) * poly_ls[0](zeta) * alpha * alpha
    re = (
        poly_t_low
        + poly_t_mid * zeta**zkey.domain_size
        + poly_t_high * zeta ** (zkey.domain_size * 2)
    ) * poly_zh_4(zeta)
    poly_r = ra + rb - rc + rd - re

    def div_poly(poly: Polynomial, x: GFElement) -> Polynomial:
        n = poly.degree()
        coeffs = poly.coeffs[: n + 1]
        ret = [Fr(0)] * n
        ret[n - 1] = coeffs[n]
        for i in range(n - 2, -1, -1):
            ret[i] = ret[i + 1] * x + coeffs[i + 1]
        assert poly[0] == ret[0] * x * Fr(-1)
        return Polynomial(poly._p, coeffs=ret)

    poly_w_zeta = (
        poly_r
        + (poly_a_blinded_4 - eval_a) * v
        + (poly_b_blinded_4 - eval_b) * v**2
        + (poly_c_blinded_4 - eval_c) * v**3
        + (poly_s1_4 - eval_s1) * v**4
        + (poly_s2_4 - eval_s2) * v**5
    )
    poly_w_zeta = div_poly(poly_w_zeta, zeta)

    poly_w_zetaw = poly_z_blinded_4 - eval_zw
    poly_w_zetaw = div_poly(poly_w_zetaw, zeta * omega)

    proof_w_zeta = exp_tau(poly_w_zeta, taus)
    proof_w_zetaw = exp_tau(poly_w_zetaw, taus)

    proof = Proof(
        a=proof_a + (FQ(1),),
        b=proof_b + (FQ(1),),
        c=proof_c + (FQ(1),),
        z=proof_z + (FQ(1),),
        t_low=proof_t1 + (FQ(1),),
        t_mid=proof_t2 + (FQ(1),),
        t_high=proof_t3 + (FQ(1),),
        w_zeta=proof_w_zeta + (FQ(1),),
        w_zetaw=proof_w_zetaw + (FQ(1),),
        eval_a=eval_a,
        eval_b=eval_b,
        eval_c=eval_c,
        eval_s1=eval_s1,
        eval_s2=eval_s2,
        eval_zw=eval_zw,
    )

    transcript6 = b""
    transcript6 += Fr.to_bytes(v)[::-1]
    transcript6 += int(proof_w_zeta[0]).to_bytes(32, "big")
    transcript6 += int(proof_w_zeta[1]).to_bytes(32, "big")
    transcript6 += int(proof_w_zetaw[0]).to_bytes(32, "big")
    transcript6 += int(proof_w_zetaw[1]).to_bytes(32, "big")
    h = keccak(transcript6)
    u = Fr(int.from_bytes(h, "big"))

    public_inputs = []
    for i in range(zkey.n_public):
        public_inputs.append(int(poly_a(omega**i)))

    if debug:
        print(f"{public_inputs = }")
        print(f"{beta = }")
        print(f"{gamma = }")
        print(f"{alpha = }")
        print(f"{zeta = }")
        print(f"{v = }")
        print(f"{u = }")

    return proof, public_inputs
